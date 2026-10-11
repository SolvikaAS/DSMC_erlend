#include "Geometry.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <limits>

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

// Axis-aligned rectangular obstacle geometry.
double BoxX = 0.0;
double BoxY = 0.0;
double BoxWidth = 2.0;
double BoxHeight = 2.0;

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


bool InsideBox(double x, double y)
{
    return std::abs(x - BoxX) < BoxWidth / 2.0
        && std::abs(y - BoxY) < BoxHeight / 2.0;
}

bool IntersectBox(double x, double y, double dx, double dy,
                  double &tHit, double &normalX, double &normalY)
{
    const double position[2] = {x, y};
    const double displacement[2] = {dx, dy};
    const double lower[2] = {BoxX - BoxWidth / 2.0, BoxY - BoxHeight / 2.0};
    const double upper[2] = {BoxX + BoxWidth / 2.0, BoxY + BoxHeight / 2.0};

    double entry = -std::numeric_limits<double>::infinity();
    double exit = std::numeric_limits<double>::infinity();
    int faceAxis = -1;
    double faceNormal = 0.0;

    // Intersect the time intervals for the x and y slabs. A segment that
    // only touches an edge or corner does not enter the solid interior.
    for (int axis = 0; axis < 2; ++axis)
    {
        if (displacement[axis] == 0.0)
        {
            if (position[axis] <= lower[axis] || position[axis] >= upper[axis])
                return false;
            continue;
        }

        double near = (lower[axis] - position[axis]) / displacement[axis];
        double far = (upper[axis] - position[axis]) / displacement[axis];
        double normal = -1.0;
        if (near > far)
        {
            std::swap(near, far);
            normal = 1.0;
        }

        // Resolve an exact corner entry consistently in favour of the x face.
        if (near > entry)
        {
            entry = near;
            faceAxis = axis;
            faceNormal = normal;
        }
        exit = std::min(exit, far);
    }

    if (faceAxis < 0 || entry < 0.0 || entry > 1.0 || entry >= exit)
        return false;

    tHit = entry;
    normalX = faceAxis == 0 ? faceNormal : 0.0;
    normalY = faceAxis == 1 ? faceNormal : 0.0;
    return true;
}

void RemoveParticlesInsideBox()
{
    for (int i = 0; i < NP;)
    {
        if (InsideBox(P[i].X[0], P[i].X[1]))
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
