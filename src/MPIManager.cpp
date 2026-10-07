#include "MPIManager.h"
#include <iostream>
#include <cmath>
#include <cstddef>

#include "SimulationState.h"

void CreateMPIParticleType()
{
    int blockLengths[2] =
    {
        DIM,   // X
        3      // V
    };

    MPI_Aint displacements[2] =
    {
        offsetof(PCL, X),
        offsetof(PCL, V)
    };

    MPI_Datatype types[2] =
    {
        MPI_DOUBLE,
        MPI_DOUBLE
    };

    MPI_Datatype RawType;

    MPI_Type_create_struct(
        2,
        blockLengths,
        displacements,
        types,
        &RawType
    );

    MPI_Type_create_resized(
        RawType,
        0,
        sizeof(PCL),
        &MPI_PCL
    );

    MPI_Type_commit(&MPI_PCL);
    MPI_Type_free(&RawType);
}


void SplitDimension(
    int GlobalN,
    int Coordinate,
    int NProcesses,
    int &LocalN,
    int &StartCell,
    int &EndCell
)
{
    int BaseCells =
        GlobalN / NProcesses;

    int Remainder =
        GlobalN % NProcesses;

    LocalN =
        BaseCells
        + (Coordinate < Remainder ? 1 : 0);

    StartCell =
        Coordinate * BaseCells
        + (Coordinate < Remainder
            ? Coordinate
            : Remainder);

    EndCell =
        StartCell + LocalN;
}

void SetupMPIDomain()
{
    // ---------------------------------------------
    // Decide MPI grid dimensions
    // ---------------------------------------------

    int dims[2] = {0, 0};

    MPI_Dims_create(
        MPISize,
        2,
        dims
    );

    // For our convention:
    // dimension 0 = x
    // dimension 1 = y

    ProcX = dims[0];
    ProcY = dims[1];


    // For a roughly square domain, prefer more
    // processes in x if MPI happened to reverse them.
    if (ProcX < ProcY)
    {
        int temp = ProcX;
        ProcX = ProcY;
        ProcY = temp;
    }

    dims[0] = ProcX;
    dims[1] = ProcY;


    // ---------------------------------------------
    // Non-periodic domain
    // ---------------------------------------------

    int periods[2] =
    {
        0,  // x not periodic
        0   // y not periodic
    };


    // ---------------------------------------------
    // Create Cartesian MPI communicator
    // ---------------------------------------------

    MPI_Cart_create(
        MPI_COMM_WORLD,
        2,
        dims,
        periods,
        0,                  // don't reorder ranks
        &CartesianComm
    );


    // ---------------------------------------------
    // Get this process's Cartesian rank
    // ---------------------------------------------

    MPI_Comm_rank(
        CartesianComm,
        &MPIRank
    );


    // ---------------------------------------------
    // Find our x/y coordinates
    // ---------------------------------------------

    int coords[2];

    MPI_Cart_coords(
        CartesianComm,
        MPIRank,
        2,
        coords
    );

    RankX = coords[0];
    RankY = coords[1];


    // ---------------------------------------------
    // Find neighbors automatically
    // ---------------------------------------------

    MPI_Cart_shift(
        CartesianComm,
        0,                  // x direction
        1,
        &LeftRank,
        &RightRank
    );

    MPI_Cart_shift(
        CartesianComm,
        1,                  // y direction
        1,
        &BottomRank,
        &TopRank
    );


    // ---------------------------------------------
    // Split global x grid
    // ---------------------------------------------

    SplitDimension(
        NX,
        RankX,
        ProcX,
        LocalNX,
        GlobalStartCellX,
        GlobalEndCellX
    );


    // ---------------------------------------------
    // Split global y grid
    // ---------------------------------------------

    SplitDimension(
        NY,
        RankY,
        ProcY,
        LocalNY,
        GlobalStartCellY,
        GlobalEndCellY
    );


    // ---------------------------------------------
    // Physical boundaries of this rank
    // ---------------------------------------------

    LocalX1 =
        X1 + GlobalStartCellX * DX;

    LocalX2 =
        X1 + GlobalEndCellX * DX;

    LocalY1 =
        Y1 + GlobalStartCellY * DY;

    LocalY2 =
        Y1 + GlobalEndCellY * DY;
}



int OwnerRankFromGlobalCell(int globalK)
{
    int baseCells = NX / MPISize;
    int remainder = NX % MPISize;

    int largeBlockCells =
        (baseCells + 1) * remainder;

    if (globalK < largeBlockCells)
    {
        return globalK / (baseCells + 1);
    }

    return remainder
        + (globalK - largeBlockCells) / baseCells;
}


void CreateMPICellType()
{
    const int count = 5;

    int blockLengths[count] =
    {
        1,  // CountNP
        3,  // CountV
        9,  // CountVV
        3,  // CountV2V
        1   // CountNC
    };

    MPI_Aint displacements[count] =
    {
        offsetof(CELL, CountNP),
        offsetof(CELL, CountV),
        offsetof(CELL, CountVV),
        offsetof(CELL, CountV2V),
        offsetof(CELL, CountNC)
    };

    MPI_Datatype types[count] =
    {
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE
    };

    MPI_Datatype RawCellType;

    MPI_Type_create_struct(
        count,
        blockLengths,
        displacements,
        types,
        &RawCellType
    );

    // Ensure arrays of CELL have the correct spacing
    MPI_Type_create_resized(
        RawCellType,
        0,
        sizeof(CELL),
        &MPI_CELL
    );

    MPI_Type_commit(&MPI_CELL);
    MPI_Type_free(&RawCellType);
}

std::vector<CELL> GatherGlobalField()
{
    // =====================================================
    // Pack this rank's local C[][] block
    // =====================================================

    int LocalCellCount =
        LocalNX * LocalNY;

    std::vector<CELL> LocalField(
        LocalCellCount
    );

    for (int k = 0; k < LocalNX; k++)
    {
        for (int l = 0; l < LocalNY; l++)
        {
            LocalField[k * LocalNY + l] =
                C[k][l];
        }
    }


    // =====================================================
    // Tell rank 0 where this block belongs
    // =====================================================

    int LocalInfo[4] =
    {
        GlobalStartCellX,
        GlobalStartCellY,
        LocalNX,
        LocalNY
    };

    std::vector<int> AllInfo;

    if (MPIRank == 0)
    {
        AllInfo.resize(
            MPISize * 4
        );
    }

    MPI_Gather(
        LocalInfo,
        4,
        MPI_INT,

        MPIRank == 0
            ? AllInfo.data()
            : nullptr,
        4,
        MPI_INT,

        0,
        CartesianComm
    );


    // =====================================================
    // Rank 0 calculates receive counts and displacements
    // =====================================================

    std::vector<int> ReceiveCounts;
    std::vector<int> Displacements;

    int TotalCells = 0;

    if (MPIRank == 0)
    {
        ReceiveCounts.resize(MPISize);
        Displacements.resize(MPISize);

        for (int r = 0; r < MPISize; r++)
        {
            int nx =
                AllInfo[4 * r + 2];

            int ny =
                AllInfo[4 * r + 3];

            ReceiveCounts[r] =
                nx * ny;

            Displacements[r] =
                TotalCells;

            TotalCells +=
                ReceiveCounts[r];
        }
    }


    // =====================================================
    // Gather all packed blocks
    // =====================================================

    std::vector<CELL> PackedGlobalField;

    if (MPIRank == 0)
    {
        PackedGlobalField.resize(
            TotalCells
        );
    }

    MPI_Gatherv(
        LocalField.data(),
        LocalCellCount,
        MPI_CELL,

        MPIRank == 0
            ? PackedGlobalField.data()
            : nullptr,

        MPIRank == 0
            ? ReceiveCounts.data()
            : nullptr,

        MPIRank == 0
            ? Displacements.data()
            : nullptr,

        MPI_CELL,
        0,
        CartesianComm
    );


    // =====================================================
    // Non-root ranks are finished
    // =====================================================

    if (MPIRank != 0)
    {
        return {};
    }


    // =====================================================
    // Reconstruct actual NX × NY global field
    // =====================================================

    std::vector<CELL> GlobalField(
        NX * NY
    );


    for (int r = 0; r < MPISize; r++)
    {
        int StartX =
            AllInfo[4 * r + 0];

        int StartY =
            AllInfo[4 * r + 1];

        int nx =
            AllInfo[4 * r + 2];

        int ny =
            AllInfo[4 * r + 3];

        int offset =
            Displacements[r];


        for (int k = 0; k < nx; k++)
        {
            for (int l = 0; l < ny; l++)
            {
                int GlobalK =
                    StartX + k;

                int GlobalL =
                    StartY + l;

                int LocalIndex =
                    offset
                    + k * ny
                    + l;

                int GlobalIndex =
                    GlobalK * NY
                    + GlobalL;

                GlobalField[GlobalIndex] =
                    PackedGlobalField[LocalIndex];
            }
        }
    }


    return GlobalField;
}


int OwnerCoordinate(
    int GlobalCell,
    int GlobalN,
    int NProcesses
)
{
    int BaseCells =
        GlobalN / NProcesses;

    int Remainder =
        GlobalN % NProcesses;

    int LargeRegion =
        (BaseCells + 1) * Remainder;

    if (GlobalCell < LargeRegion)
    {
        return
            GlobalCell / (BaseCells + 1);
    }

    return
        Remainder
        +
        (GlobalCell - LargeRegion)
        / BaseCells;
}


int ParticleOwnerRank(
    double x,
    double y
)
{
    int globalK =
        static_cast<int>(
            floor((x - X1) / DX)
        );

    int globalL =
        static_cast<int>(
            floor((y - Y1) / DY)
        );


    // Particle must already be inside the
    // physical simulation domain.
    if (
        globalK < 0 ||
        globalK >= NX ||
        globalL < 0 ||
        globalL >= NY
    )
    {
        return MPI_PROC_NULL;
    }


    int ownerX =
        OwnerCoordinate(
            globalK,
            NX,
            ProcX
        );

    int ownerY =
        OwnerCoordinate(
            globalL,
            NY,
            ProcY
        );


    int coords[2] =
    {
        ownerX,
        ownerY
    };

    int ownerRank;

    MPI_Cart_rank(
        CartesianComm,
        coords,
        &ownerRank
    );

    return ownerRank;
}

void ExchangeParticles()
{
    // One outgoing bucket for every possible MPI rank
    std::vector<std::vector<PCL>>
        SendBuckets(MPISize);

    // Particles which still belong to this rank
    std::vector<PCL> KeepParticles;

    KeepParticles.reserve(NP);


    // =================================================
    // Decide destination of every local particle
    // =================================================

    for (int i = 0; i < NP; i++)
    {
        int owner =
            ParticleOwnerRank(
                P[i].X[0],
                P[i].X[1]
            );


        if (owner == MPI_PROC_NULL)
        {
            std::cerr
                << "Rank " << MPIRank
                << ": ExchangeParticles received "
                << "particle outside physical domain\n"
                << "x = " << P[i].X[0]
                << '\n'
                << "y = " << P[i].X[1]
                << '\n';

            MPI_Abort(
                CartesianComm,
                30
            );
        }


        if (owner == MPIRank)
        {
            KeepParticles.push_back(P[i]);
        }
        else
        {
            SendBuckets[owner].push_back(P[i]);
        }
    }


    // =================================================
    // Number of particles sent to every rank
    // =================================================

    std::vector<int>
        SendCounts(MPISize, 0);

    std::vector<int>
        ReceiveCounts(MPISize, 0);


    for (int r = 0; r < MPISize; r++)
    {
        SendCounts[r] =
            static_cast<int>(
                SendBuckets[r].size()
            );
    }


    // Everyone tells everyone how many particles
    // they are about to send.
    MPI_Alltoall(
        SendCounts.data(),
        1,
        MPI_INT,

        ReceiveCounts.data(),
        1,
        MPI_INT,

        CartesianComm
    );


    // =================================================
    // Build displacements
    // =================================================

    std::vector<int>
        SendDisplacements(MPISize, 0);

    std::vector<int>
        ReceiveDisplacements(MPISize, 0);


    int TotalSend = 0;
    int TotalReceive = 0;


    for (int r = 0; r < MPISize; r++)
    {
        SendDisplacements[r] =
            TotalSend;

        TotalSend +=
            SendCounts[r];


        ReceiveDisplacements[r] =
            TotalReceive;

        TotalReceive +=
            ReceiveCounts[r];
    }


    // =================================================
    // Flatten outgoing particle buffers
    // =================================================

    std::vector<PCL>
        SendBuffer(TotalSend);

    for (int r = 0; r < MPISize; r++)
    {
        int offset =
            SendDisplacements[r];

        for (
            int i = 0;
            i < SendCounts[r];
            i++
        )
        {
            SendBuffer[offset + i] =
                SendBuckets[r][i];
        }
    }


    // =================================================
    // Allocate incoming buffer
    // =================================================

    std::vector<PCL>
        ReceiveBuffer(TotalReceive);


    // =================================================
    // Exchange actual particles
    // =================================================

    MPI_Alltoallv(
        SendBuffer.empty()
            ? nullptr
            : SendBuffer.data(),

        SendCounts.data(),
        SendDisplacements.data(),
        MPI_PCL,

        ReceiveBuffer.empty()
            ? nullptr
            : ReceiveBuffer.data(),

        ReceiveCounts.data(),
        ReceiveDisplacements.data(),
        MPI_PCL,

        CartesianComm
    );


    // =================================================
    // Check local particle capacity
    // =================================================

    int NewNP =
        static_cast<int>(
            KeepParticles.size()
        )
        +
        TotalReceive;


    if (NewNP > MAX_PCL)
    {
        std::cerr
            << "Rank " << MPIRank
            << ": MAX_PCL exceeded "
            << "during 2D MPI exchange\n"
            << "NewNP = "
            << NewNP
            << '\n';

        MPI_Abort(
            CartesianComm,
            31
        );
    }


    // =================================================
    // Rebuild local P[]
    // =================================================

    NP = 0;


    for (
        const PCL &particle
        : KeepParticles
    )
    {
        P[NP++] =
            particle;
    }


    for (
        const PCL &particle
        : ReceiveBuffer
    )
    {
        P[NP++] =
            particle;
    }
}
