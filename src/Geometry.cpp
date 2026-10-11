#include "Geometry.h"
#include <iostream>
#include <cstring>
#include <cmath>

#include "SimulationState.h"

// Domain bounds and Cartesian grid settings.
double X1 = -20.0;
double X2 = 20.0;
double Y1 = -20.0;
double Y2 = 20.0;
double DZ = 0.1;
int NX = 100;
int NY = 100;
double DL = 0.1;

// Circular obstacle geometry.
double SphereX = 0.0;
double SphereY = 0.0;
double SphereRadius = 3;

// Derived cell dimensions and volume, calculated during setup.
double DX, DY;
double DV;

void SetupGeometry()
{
	DX = (X2-X1)/NX;
	DY = (Y2-Y1)/NY;
	DV=DX*DY*DZ;
}

void Indexing()
{
    memset(NPC, 0, sizeof NPC);

    for (int i = 0; i < NP; i++)
    {
        // -----------------------------------------
        // Global cell coordinates
        // -----------------------------------------

        int globalK =
            static_cast<int>(
                floor(
                    (P[i].X[0] - X1) / DX
                )
            );

        int globalL =
            static_cast<int>(
                floor(
                    (P[i].X[1] - Y1) / DY
                )
            );


        // -----------------------------------------
        // Check global physical domain
        // -----------------------------------------

        if (
            globalK < 0 ||
            globalK >= NX ||
            globalL < 0 ||
            globalL >= NY
        )
        {
            std::cerr
                << "Rank " << MPIRank
                << ": particle outside global domain\n"
                << "particle = " << i << '\n'
                << "x = " << P[i].X[0] << '\n'
                << "y = " << P[i].X[1] << '\n'
                << "globalK = " << globalK << '\n'
                << "globalL = " << globalL << '\n';

            MPI_Abort(
                MPI_COMM_WORLD,
                20
            );
        }


        // -----------------------------------------
        // Check MPI ownership
        // -----------------------------------------

        if (
            globalK < GlobalStartCellX ||
            globalK >= GlobalEndCellX ||
            globalL < GlobalStartCellY ||
            globalL >= GlobalEndCellY
        )
        {
            std::cerr
                << "Rank " << MPIRank
                << ": particle outside local MPI domain\n"
                << "particle = " << i << '\n'
                << "x = " << P[i].X[0] << '\n'
                << "y = " << P[i].X[1] << '\n'
                << "global cell = ("
                << globalK << ", "
                << globalL << ")\n"
                << "owned X cells = "
                << GlobalStartCellX
                << " ... "
                << GlobalEndCellX - 1
                << '\n'
                << "owned Y cells = "
                << GlobalStartCellY
                << " ... "
                << GlobalEndCellY - 1
                << '\n';

            MPI_Abort(
                MPI_COMM_WORLD,
                21
            );
        }


        // -----------------------------------------
        // Convert global → local cell coordinates
        // -----------------------------------------

        int k =
            globalK - GlobalStartCellX;

        int l =
            globalL - GlobalStartCellY;


        // -----------------------------------------
        // Cell capacity check
        // -----------------------------------------

        if (NPC[k][l] >= MAX_PCL_IN_CELL)
        {
            std::cerr
                << "Rank " << MPIRank
                << ": too many particles in cell\n"
                << "local cell = ("
                << k << ", "
                << l << ")\n"
                << "global cell = ("
                << globalK << ", "
                << globalL << ")\n"
                << "NPC = "
                << NPC[k][l]
                << '\n';

            MPI_Abort(
                MPI_COMM_WORLD,
                22
            );
        }


        // -----------------------------------------
        // Store particle index
        // -----------------------------------------

        IPC[k][l][NPC[k][l]++] = i;
    }
}


bool InsideSphere(double x, double y)
{
    double dx = x - SphereX;
    double dy = y - SphereY;

    return dx*dx + dy*dy
           < SphereRadius*SphereRadius;
}

void RemoveParticlesInsideSphere()
{
    for (int i = 0; i < NP;)
    {
        if (InsideSphere(P[i].X[0], P[i].X[1]))
        {
            if (i < NP - 1)
                P[i] = P[NP - 1];

            NP--;
        }
        else
        {
            i++;
        }
    }
}
