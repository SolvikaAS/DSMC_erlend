#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <vector>
#include <cstddef>


#include "BRNG.cxx"
#include "DiscreteUniformRNG.cxx"
#include "PoissonRNG.cxx"
#include "RayleighRNG.cxx"
#include "IsotropicVectorRNG.cxx"
#include "Gamma.cxx"

#define MAX_PCL_IN_CELL 1000
#define MAX_CELL_X 200
#define MAX_CELL_Y 200
#define DIM 2
#define MAX_PCL 1000000

#define M_PI                    3.14159265358979323846
#define BOLTZMANN_CONSTANT      1.380662E-23                    // Boltzmann constant (J/K) 
#define AVOGADRO_CONSTANT       6.022045E+23                    // Avogadro constant (1/mole)
#define UNIT_DALTON             1.66056e-27                     // Atomic mass unit, Dalton (kg)

double Dt = 1.0e-05;
int NStep = 100000;
int FirstSamplingStep = 10000;
int SamplingPeriod = 1;
int PrintPeriod = 10000;
double X1 = -20.0;
double X2 = 20.0;
double Y1 = -20.0;
double Y2 = 20.0;
double DZ = 0.1;
int NX = 100;
int NY = 100;
double NPCFree = 10;
double DL = 0.1;
double MolarMass = 0.040;
double MaFree = 4;
double NFree = 1e18;
double TFree = 200;
double SphereX = 0.0;
double SphereY = 0.0;
double SphereRadius = 3;
double AttackAngle = 0;
double Tw = 300.0;
double Omega;
double SigmaRef;
double CrRef;
double MFTFree;
double MFPFree;
double KnFree;
double DRef = 4.17e-10;
double ViscosityIndex = 0.81;
double Tref = 273.0;
double DragForce = 0;
double DragForceSum = 0;
int DragSamples =0;
double DragImpulse = 0;
long long DragCollisionCount = 0;
double DragImpulseTotal = 0.0;

MPI_Datatype MPI_CELL;
MPI_Datatype MPI_PCL;

int MPIRank;
int MPISize;

// Cartesian MPI communicator
MPI_Comm CartesianComm = MPI_COMM_NULL;

// Number of MPI domains in each direction
int ProcX;
int ProcY;

// This rank's position in the MPI grid
int RankX;
int RankY;

// Neighbor ranks
int LeftRank;
int RightRank;
int BottomRank;
int TopRank;

// Local number of DSMC cells
int LocalNX;
int LocalNY;

// Global cell range owned by this rank
int GlobalStartCellX;
int GlobalEndCellX;

int GlobalStartCellY;
int GlobalEndCellY;

// Physical boundaries of this MPI subdomain
double LocalX1;
double LocalX2;

double LocalY1;
double LocalY2;

typedef struct pcl {
	double X[DIM];
	double V[3];
} PCL;


typedef struct cell { //////////////////////////////////////////////////////////////////////
                double  CountNP;                // Particle number -> Number density
                double  CountV[3];              // Velocity vector V -> Macroscopic gas velocity
                double  CountVV[3][3];          // Velocity tensor V x V -> Stress tensor
                double  CountV2V[3];            // Vector of V^2 * V -> Heat flux vector
                double  CountNC;                // Collsion number -> Collsion frequency, M.F.T. and M.F.P.
} CELL;





int SampleStep;
CELL C[MAX_CELL_X][MAX_CELL_Y];
#define sqr3( V ) ( V[0] * V[0] + V[1] * V[1] + V[2] * V[2] )


int NP;
PCL P[MAX_PCL];



double DX, DY;
double DV;
double Weight;
double MoleculeMass;

double PFree;
double UFree;
double UxFree, UyFree;

int NPC[MAX_CELL_X][MAX_CELL_Y];
int IPC[MAX_CELL_X][MAX_CELL_Y][MAX_PCL_IN_CELL];

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


void Setup(){
	MoleculeMass =MolarMass/AVOGADRO_CONSTANT;
	PFree= NFree*(BOLTZMANN_CONSTANT*TFree);
	UFree = MaFree*sqrt((5.0/3.0)*BOLTZMANN_CONSTANT*TFree/MoleculeMass);
	UxFree =UFree*cos(M_PI*AttackAngle/180.0);
	UyFree =-UFree*sin(M_PI*AttackAngle/180.0);
	DX = (X2-X1)/NX;
	DY = (Y2-Y1)/NY;
	DV=DX*DY*DZ;
	Weight = DV*NFree/NPCFree;
	SigmaRef= M_PI*DRef*DRef;
	Omega = 2.0*ViscosityIndex-1.0;
	CrRef = sqrt(4.0*BOLTZMANN_CONSTANT*Tref/MoleculeMass)/pow(Gamma(2.5-ViscosityIndex),1.0/Omega);
	double q = 1.0 - Omega;

double MeanRelativeMoment =
    pow(
        4.0 * BOLTZMANN_CONSTANT * TFree / MoleculeMass,
        q / 2.0
    )
    *
    Gamma((q + 3.0) / 2.0)
    /
    Gamma(1.5);

double CollisionRateCoefficient =
    SigmaRef
    *
    pow(CrRef, Omega)
    *
    MeanRelativeMoment;

double CollisionFrequency =
    NFree * CollisionRateCoefficient;

MFTFree = 1.0 / CollisionFrequency;

double MeanMolecularSpeed =
    sqrt(
        8.0 * BOLTZMANN_CONSTANT * TFree
        /
        (M_PI * MoleculeMass)
    );

MFPFree = MeanMolecularSpeed * MFTFree;

KnFree =
    MFPFree / (2.0 * SphereRadius);

/*
std::cout << "\n";
std::cout << "Freestream properties\n";
std::cout << "-----------------------------\n";
std::cout << "P_inf       = " << PFree << " 1/m^3\n";
std::cout << "U_inf       = " << UFree << " m/s\n";
std::cout << "MFT_inf     = " << MFTFree << " s\n";
std::cout << "MFP_inf     = " << MFPFree << " m\n";
std::cout << "Kn_inf      = " << KnFree << "\n";
std::cout << "DX / lambda = " << DX / MFPFree << "\n";
std::cout << "DY / lambda = " << DY / MFPFree << "\n";
std::cout << "Dt / MFT    = " << Dt / MFTFree << "\n";
*/

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


double SampleTime;

void Sampling () ////////////////////////////////////////////////////////////////////////////
        {
                
                SampleStep++;
                SampleTime += Dt;
                for ( int k = 0; k < LocalNX; k++ )
                        for ( int l = 0; l < LocalNY; l++ ) 
                                for ( int i = 0; i < NPC[k][l]; i++ ) {
                                        C[k][l].CountNP += 1.0;
                                        double V2 = sqr3 ( P[IPC[k][l][i]].V );
                                        for ( int m = 0; m < 3; m++ ) {
                                                C[k][l].CountV[m] += P[IPC[k][l][i]].V[m];
                                                C[k][l].CountV2V[m] += V2 * P[IPC[k][l][i]].V[m];
                                                for ( int n = 0; n < 3; n++ )
                                                        C[k][l].CountVV[m][n] += P[IPC[k][l][i]].V[m] * P[IPC[k][l][i]].V[n];
                                        }
                        }
        }




void Printing(
    const std::vector<CELL> &GlobalField
)
{
    if (MPIRank != 0)
    {
        return;
    }


    if (
        GlobalField.size()
        != static_cast<size_t>(NX * NY)
    )
    {
        std::cerr
            << "ERROR: GlobalField has wrong size\n"
            << "Expected = "
            << NX * NY
            << '\n'
            << "Received = "
            << GlobalField.size()
            << '\n';

        return;
    }


    FILE *F =
        fopen(
            "DSMC2D_Field.plt",
            "wt"
        );


    if (F == nullptr)
    {
        std::cerr
            << "ERROR: Could not open "
            << "DSMC2D_Field.plt\n";

        return;
    }


    fprintf(
        F,
        "VARIABLES = "
        "\"X\" \"Y\" \"N\" "
        "\"Ux\" \"Uy\" "
        "\"Tx\" \"Ty\" \"T\"\n"
    );

    fprintf(
        F,
        "ZONE I=%d J=%d\n",
        NX,
        NY
    );


    // =====================================================
    // Write complete global field
    // =====================================================

    for (int k = 0; k < NX; k++)
    {
        for (int l = 0; l < NY; l++)
        {
            const CELL &cell =
                GlobalField[
                    k * NY + l
                ];


            double X =
                X1 + (k + 0.5) * DX;

            double Y =
                Y1 + (l + 0.5) * DY;


            double N = 0.0;

            double U[3] =
            {
                0.0,
                0.0,
                0.0
            };

            double Txyz[3] =
            {
                0.0,
                0.0,
                0.0
            };


            // =================================================
            // Number density
            // =================================================

            if (SampleStep > 0)
            {
                N =
                    Weight
                    * cell.CountNP
                    / SampleStep
                    / DV;
            }


            // =================================================
            // Velocity and temperature
            // =================================================

            if (cell.CountNP > 0.0)
            {
                for (int m = 0; m < 3; m++)
                {
                    U[m] =
                        cell.CountV[m]
                        / cell.CountNP;


                    Txyz[m] =
                        MoleculeMass
                        *
                        (
                            cell.CountVV[m][m]
                            / cell.CountNP

                            - U[m] * U[m]
                        )
                        /
                        BOLTZMANN_CONSTANT;
                }
            }


            double T =
                (
                    Txyz[0]
                    + Txyz[1]
                    + Txyz[2]
                )
                / 3.0;


            fprintf(
                F,

                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e\n",

                X,
                Y,
                N,
                U[0],
                U[1],
                Txyz[0],
                Txyz[1],
                T
            );
        }
    }


    fclose(F);

    std::cout
        << "Global field written to "
        << "DSMC2D_Field.plt\n";
}


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

int Step;
double Time;
//double SampleTime;



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

void ElasticCollision ( double *V1, double *V2, double Cr){
	double N[3];
	v3rand_isotropic (N);
	Cr *=0.5;
	double VC[3] = {0.5*(V2[0]+V1[0]),0.5*(V2[1]+V1[1]),0.5*(V2[2]+V1[2])};
	double VCr[3] = {Cr*N[0],Cr*N[1],Cr*N[2]};
	V1[0]=VC[0]+VCr[0];
	V1[1]=VC[1]+VCr[1];
	V1[2]=VC[2]+VCr[2];
	V2[0]=VC[0]-VCr[0];
	V2[1]=VC[1]-VCr[1];
	V2[2]=VC[2]-VCr[2];
}

int BinaryCollisionsInCell (double CellVolume,int NPC,int *IPC){
	if (NPC<2)return 0;
	double SCMax = 9.0*SigmaRef*sqrt(BOLTZMANN_CONSTANT*TFree/MoleculeMass);
	double Npair = 0.5*NPC*(NPC-1)*Weight*SCMax*Dt/CellVolume;
	int NN=int (Npair);
	double N1= Npair-NN;
	NN =(brng()<N1)?NN+1:NN;
	int i,j;
	int NC=0;
	for (int k=0; k<NN;k++){
		i = irand_uniform(NPC);
		do{j = irand_uniform(NPC);}while (j==i);
		i = IPC[i];
		j = IPC[j];
		double VCr[3] = {P[j].V[0]-P[i].V[0],P[j].V[1]-P[i].V[1],P[j].V[2]-P[i].V[2]};
		double Cr = sqrt(VCr[0]*VCr[0]+VCr[1]*VCr[1]+VCr[2]*VCr[2]);
		double Sigma = SigmaRef*pow(CrRef/Cr, Omega);
		if ( brng()<Sigma*Cr/SCMax){
			ElasticCollision (P[i].V,P[j].V,Cr);
			NC++;
		}
	}
	return NC;
}
int SamplingMode;

int CollideParticles () ////////////////////////////////////////////////////////////////////
        {
        int     NC = 0; // Counter of collisions
                for ( int k = 0; k < LocalNX; k++ )
                        for ( int l = 0; l < LocalNY; l++ ) { // Sampling of collisions in particle (i,j)
                                int NCC = BinaryCollisionsInCell ( DV, NPC[k][l], IPC[k][l] ); // NCC is the number of collisions sampled in the current cell
                                NC += NCC;
                                // Here we update CountNC in the cell for fure calculations of the mean free time and path, if necessary
                                if ( SamplingMode == 1 ) C[k][l].CountNC += NCC;
                        }
                return NC; // NC is the total number of real collsions sampled during the time step
}

void InitialConditions()
{
    Step = 0;
    Time = 0.0;
    NP = 0;

    SampleStep = 0;
    SampleTime = 0.0;

    memset(C, 0, sizeof C);

    // Generate particles only inside this MPI rank's domain
    GenerateNewParticles(
        LocalX1, LocalY1,
        LocalX2, LocalY2,
        NFree,
        UxFree,
        UyFree,
        TFree,
        0
    );

    RemoveParticlesInsideSphere();
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

double DragSamplingStart = 0.0;
double DragSamplingTime = 0.0;

int main(int argc, char **argv)
{

        MPI_Init(&argc,&argv);
        
        CreateMPICellType();
        CreateMPIParticleType();
        
        MPI_Comm_rank(MPI_COMM_WORLD, &MPIRank);
        MPI_Comm_size(MPI_COMM_WORLD, &MPISize);

    // Give each MPI process its own RNG seed
    time_t t;
    unsigned baseSeed = static_cast<unsigned>(time(&t));

    SetSeed(
        baseSeed + 7919u * MPIRank,
        362436069u + 104729u * MPIRank
    );

	Setup();
	std::cout<<"Setup \n";
	SetupMPIDomain();
	
	/*
	std::cout
        << "Rank " << MPIRank
        << " / " << MPISize
        << " : cells "
        << GlobalStartCell
        << " to "
        << GlobalEndCell - 1
        << ", LocalNX = "
        << LocalNX
        << ", x = "
        << LocalX1
        << " to "
        << LocalX2
        << '\n';
	*/
	
	InitialConditions();
	std::cout
    << "Rank " << MPIRank
    << " initial particles = "
    << NP
    << '\n';
    
    int GlobalNP = 0;

    MPI_Reduce(
        &NP,
        &GlobalNP,
        1,
        MPI_INT,
        MPI_SUM,
        0,
        MPI_COMM_WORLD
    );

    if (MPIRank == 0)
    {
        std::cout
            << "Total particles across all ranks = "
            << GlobalNP
            << '\n';
    }
    
	
	std::cout<<"InitialConditions \n";
	do {
		if (Step > FirstSamplingStep && Step % SamplingPeriod ==0)
			SamplingMode =1;
		else
			SamplingMode =0;
		MoveParticles();
		if(Step>FirstSamplingStep){
			DragSamplingTime +=Dt;
		}
		//std::cout<<"moved \n";
		BoundaryConditions();
		//std::cout<<"boundary \n";
		
		ExchangeParticles();
		
		Indexing();
		//std::cout<<"indexed \n";
		int NC = CollideParticles();
		if (SamplingMode ==1)
		{
		 Sampling();
		}
		Step++;
		Time += Dt;
		if (Step % 100 == 0)
                {
                    int GlobalNP = 0;
                    int GlobalNC = 0;

                    MPI_Reduce(
                        &NP,
                        &GlobalNP,
                        1,
                        MPI_INT,
                        MPI_SUM,
                        0,
                        MPI_COMM_WORLD
                    );

                    MPI_Reduce(
                        &NC,
                        &GlobalNC,
                        1,
                        MPI_INT,
                        MPI_SUM,
                        0,
                        MPI_COMM_WORLD
                    );

                    if (MPIRank == 0)
                    {
                        std::cout
                            << "Time step " << Step
                            << ", Number of particles "
                            << GlobalNP
                            << ", Number of collisions "
                            << GlobalNC
                            << '\n';
                    }
                }
		//if (Step > FirstSamplingStep && Step % PrintPeriod ==0) Printing();
	} while (Step < NStep);
	
	std::vector<CELL> FinalField =
        GatherGlobalField();
	
	
	if (MPIRank== 0)
      {
        Printing(FinalField);
      }      

/*

	double MeanDrag =
    DragImpulseTotal / DragSamplingTime;

double RhoFree = NFree * MoleculeMass;
double AreaRef = 2.0 * SphereRadius * DZ;

double Cd =
    MeanDrag /
    (0.5 * RhoFree * UFree * UFree * AreaRef);

std::cout << "\n";
std::cout << "====================================\n";
std::cout << "Sphere collisions = "
          << DragCollisionCount << "\n";

std::cout << "Sampling time = "
          << DragSamplingTime << " s\n";

std::cout << "Mean drag force = "
          << MeanDrag << " N\n";

std::cout << "Drag coefficient Cd = "
          << Cd << "\n";

std::cout << "====================================\n";
*/

MPI_Type_free(&MPI_PCL);
MPI_Type_free(&MPI_CELL);

MPI_Comm_free(&CartesianComm);

MPI_Finalize();

return 0;
}
