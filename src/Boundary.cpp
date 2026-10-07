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
    for (int i = 0; i < NP; i++)
    {
        // --------------------------------------------
        // Old particle position
        // --------------------------------------------

        double X0 = P[i].X[0];
        double Y0 = P[i].X[1];

        // --------------------------------------------
        // Proposed new particle position
        // --------------------------------------------

        double Xnew = X0 + Dt * P[i].V[0];
        double Ynew = Y0 + Dt * P[i].V[1];

        // Particle displacement during timestep
        double dx = Xnew - X0;
        double dy = Ynew - Y0;

        // Vector from sphere center to old position
        double fx = X0 - SphereX;
        double fy = Y0 - SphereY;

        // --------------------------------------------
        // Line-circle intersection
        //
        // |P0 + t*d - C|^2 = R^2
        //
        // 0 <= t <= 1 means collision occurred
        // during this timestep.
        // --------------------------------------------

        double a = dx*dx + dy*dy;
        double b = 2.0 * (fx*dx + fy*dy);
        double c = fx*fx + fy*fy
                 - SphereRadius*SphereRadius;

        double discriminant = b*b - 4.0*a*c;

        bool collision = false;
        double tHit = 1.0;

        if (discriminant >= 0.0 && a > 0.0)
        {
            double sqrtD = sqrt(discriminant);

            double t1 = (-b - sqrtD) / (2.0*a);
            double t2 = (-b + sqrtD) / (2.0*a);

            // First intersection encountered
            if (t1 >= 0.0 && t1 <= 1.0)
            {
                tHit = t1;
                collision = true;
            }
            else if (t2 >= 0.0 && t2 <= 1.0)
            {
                tHit = t2;
                collision = true;
            }
        }

        // --------------------------------------------
        // No collision
        // --------------------------------------------

        if (!collision)
        {
            P[i].X[0] = Xnew;
            P[i].X[1] = Ynew;
            continue;
        }

        // --------------------------------------------
        // Collision position
        // --------------------------------------------

        double Xhit = X0 + tHit * dx;
        double Yhit = Y0 + tHit * dy;

        // --------------------------------------------
        // Surface normal
        // --------------------------------------------

        double Nx = (Xhit - SphereX) / SphereRadius;
        double Ny = (Yhit - SphereY) / SphereRadius;

        // Numerical normalization
        double Nmag = sqrt(Nx*Nx + Ny*Ny);

        Nx /= Nmag;
        Ny /= Nmag;

        // --------------------------------------------
        // Diffuse reflection
        // --------------------------------------------

	double VxIn = P[i].V[0];
	double VyIn = P[i].V[1];
	double VzIn = P[i].V[2];

        DiffuseScattering(
            P[i].V,
            MoleculeMass,
            Tw,
            Nx,
            Ny
        );

	double VxOut = P[i].V[0];
	DragImpulseTotal += Weight*MoleculeMass*(VxIn-VxOut);
	DragCollisionCount++;

        // --------------------------------------------
        // Remaining fraction of timestep
        // --------------------------------------------

        double DtRemaining = Dt * (1.0 - tHit);

        P[i].X[0] =
            Xhit + DtRemaining * P[i].V[0];

        P[i].X[1] =
            Yhit + DtRemaining * P[i].V[1];

        // Tiny displacement away from surface to avoid
        // numerical re-collision at exactly R
        constexpr double epsilon = 1e-12;

        P[i].X[0] += epsilon * Nx;
        P[i].X[1] += epsilon * Ny;
    }
}
