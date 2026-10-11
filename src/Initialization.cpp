#include "Initialization.h"
#include <iostream>
#include <cstring>
#include <cmath>

#include "SimulationState.h"
#include "Boundary.h"
#include "Geometry.h"
#include "Gamma.h"

double Dt = 1.0e-05;
int NStep = 100000;
int FirstSamplingStep = 10000;
int SamplingPeriod = 1;
int PrintPeriod = 10000;
double NPCFree = 10;
double MolarMass = 0.040;
double MaFree = 4;
double NFree = 1e18;
double TFree = 200;
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
int SampleStep;
CELL C[MAX_CELL_X][MAX_CELL_Y];
int NP;
PCL P[MAX_PCL];



double Weight;
double MoleculeMass;

double PFree;
double UFree;
double UxFree, UyFree;

int NPC[MAX_CELL_X][MAX_CELL_Y];
int IPC[MAX_CELL_X][MAX_CELL_Y][MAX_PCL_IN_CELL];
double SampleTime;
int Step;
double Time;
int SamplingMode;
double DragSamplingStart = 0.0;
double DragSamplingTime = 0.0;

void Setup(){
	MoleculeMass =MolarMass/AVOGADRO_CONSTANT;
	PFree= NFree*(BOLTZMANN_CONSTANT*TFree);
	UFree = MaFree*sqrt((5.0/3.0)*BOLTZMANN_CONSTANT*TFree/MoleculeMass);
	UxFree =UFree*cos(M_PI*AttackAngle/180.0);
	UyFree =-UFree*sin(M_PI*AttackAngle/180.0);
	SetupGeometry();
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
