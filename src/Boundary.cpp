#include "Boundary.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "SimulationState.h"
#include "RNG.h"

void v2rand_uniform_rect (double &X, double &Y, double a, double b, double c,double d){
	X = a+(b-a)*brng();
	Y = c+(d-c)*brng();
}

double frand_Gaussian (double E, double V){
	return E+sqrt(-2.0*V*log(brng()))*cos(2.0*M_PI*brng());}

void vrand_MB (double *v, double m,double *u,double T){
	double RT = BOLTZMANN_CONSTANT*T/m;
	v[0] = frand_Gaussian (u[0],RT);
	v[1] = frand_Gaussian (u[1],RT);
	v[2] = frand_Gaussian (u[2],RT);
}


void GenerateNewParticles (double X1, double Y1, double X2, double Y2,double Ngas,double Ugasx,
			   double Ugasy,double Tgas, int MoveFlag){
	double V = (X2-X1)*(Y2-Y1)*DZ;
	double Nnew_avg = Ngas*V/Weight;
	int Nnew;
	if (Nnew_avg>20.0){
		Nnew = int (Nnew_avg);
		if (brng()<Nnew_avg-Nnew) Nnew++;
	}
	else {
		Nnew = irand_Poisson (Nnew_avg);
	}

	if (NP + Nnew > MAX_PCL)
    {
        std::cerr << "ERROR: Particle array capacity exceeded\n";
        std::cerr << "NP = " << NP << '\n';
        std::cerr << "Nnew = " << Nnew << '\n';
        std::cerr << "MAX_PCL = " << MAX_PCL << '\n';
        std::exit(1);
    }

	double Ugas[3] = {Ugasx,Ugasy,0.0};
	for (int i =NP; i<NP+Nnew;i++){
		v2rand_uniform_rect(P[i].X[0],P[i].X[1],X1,X2,Y1,Y2);
		vrand_MB (P[i].V,MoleculeMass,Ugas,Tgas);
		if (MoveFlag ==1){
			P[i].X[0] += Dt*P[i].V[0];
			P[i].X[1] += Dt*P[i].V[1];
		}
	}
	NP +=Nnew;
}

void BoundaryConditions()
{
    // =============================================
    // LEFT physical boundary
    // =============================================

    if (LeftRank == MPI_PROC_NULL)
    {
        GenerateNewParticles(
            X1 - DL,
            LocalY1,
            X1,
            LocalY2,
            NFree,
            UxFree,
            UyFree,
            TFree,
            1
        );
    }


    // =============================================
    // RIGHT physical boundary
    // =============================================

    if (RightRank == MPI_PROC_NULL)
    {
        GenerateNewParticles(
            X2,
            LocalY1,
            X2 + DL,
            LocalY2,
            NFree,
            UxFree,
            UyFree,
            TFree,
            1
        );
    }


    // =============================================
    // BOTTOM physical boundary
    // =============================================

    if (BottomRank == MPI_PROC_NULL)
    {
        GenerateNewParticles(
            LocalX1,
            Y1 - DL,
            LocalX2,
            Y1,
            NFree,
            UxFree,
            UyFree,
            TFree,
            1
        );
    }


    // =============================================
    // TOP physical boundary
    // =============================================

    if (TopRank == MPI_PROC_NULL)
    {
        GenerateNewParticles(
            LocalX1,
            Y2,
            LocalX2,
            Y2 + DL,
            NFree,
            UxFree,
            UyFree,
            TFree,
            1
        );
    }


    // =============================================
    // Remove particles which left the PHYSICAL
    // simulation domain
    // =============================================

    for (int i = 0; i < NP;)
    {
        if (
            P[i].X[0] <= X1 ||
            P[i].X[0] >= X2 ||
            P[i].X[1] <= Y1 ||
            P[i].X[1] >= Y2
        )
        {
            if (i < NP - 1)
            {
                P[i] =
                    P[NP - 1];
            }

            NP--;
        }
        else
        {
            i++;
        }
    }
}

void DiffuseScattering(
    double *v,
    double m,
    double Tw,
    double Nx,
    double Ny
)
{
    double RTw = BOLTZMANN_CONSTANT * Tw / m;

    // Normal velocity: molecules leave the wall
    double Vn = frand_Rayleigh(sqrt(RTw));

    // Tangential velocity
    double Vt = frand_Gaussian(0.0, RTw);

    // z component
    double Vz = frand_Gaussian(0.0, RTw);

    // Unit tangent perpendicular to the normal
    double Tx = -Ny;
    double Ty =  Nx;

    // Transform local normal/tangential velocities
    // back into global x-y coordinates
    v[0] = Vn * Nx + Vt * Tx;
    v[1] = Vn * Ny + Vt * Ty;
    v[2] = Vz;
}

void MoveParticles()
{
    for (int i = 0; i < NP; ++i)
    {
        const double x0 = P[i].X[0];
        const double y0 = P[i].X[1];
        const double dx = Dt * P[i].V[0];
        const double dy = Dt * P[i].V[1];
        double tHit, normalX, normalY;

        if (!IntersectBox(x0, y0, dx, dy, tHit, normalX, normalY))
        {
            P[i].X[0] = x0 + dx;
            P[i].X[1] = y0 + dy;
            continue;
        }

        const double xHit = x0 + tHit * dx;
        const double yHit = y0 + tHit * dy;
        const double vxIn = P[i].V[0];

        DiffuseScattering(P[i].V, MoleculeMass, Tw, normalX, normalY);
        DragImpulseTotal += Weight * MoleculeMass * (vxIn - P[i].V[0]);
        ++DragCollisionCount;

        // The outgoing normal velocity keeps the rest of the trajectory
        // outside this convex box, including after a corner entry.
        const double remainingTime = Dt * (1.0 - tHit);
        constexpr double epsilon = 1e-12;
        P[i].X[0] = xHit + remainingTime * P[i].V[0] + epsilon * normalX;
        P[i].X[1] = yHit + remainingTime * P[i].V[1] + epsilon * normalY;
    }
}
