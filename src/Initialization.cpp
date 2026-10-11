#include "Initialization.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include <cfloat>

#include "SimulationState.h"
#include "Boundary.h"
#include "Geometry.h"

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

namespace {
////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculation of the Gamma function.
//--------------------------------------------------------------------------------------------------
// This code is based on the code downloaded from http://www.johndcook.com/cpp_gamma.html.
//--------------------------------------------------------------------------------------------------
// Alexey Volkov, University of Alabama, avolkov1@ua.edu, 2016, Version 12.01
////////////////////////////////////////////////////////////////////////////////////////////////////

double Gamma ( double x ) //////////////////////////////////////////////////////////////////
// x must be positive
{ //////////////////////////////////////////////////////////////////////////////////////////
// Split the function domain into three intervals:
// (0, 0.001), [0.001, 12), and (12, infinity)
////////////////////////////////////////////////////////////////////////////////////////////
// First interval: (0, 0.001)
//
// For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
// So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
// The relative error over this interval is less than 6e-7.
const   double gamma = 0.577215664901532860606512090; // Euler's gamma constant
        if ( x < 0.001 ) return 1.0 / ( x * ( 1.0 + gamma * x ) );

////////////////////////////////////////////////////////////////////////////////////////////
// Second interval: [0.001, 12)
        if ( x < 12.0 ) {
        // The algorithm directly approximates gamma over (1,2) and uses
        // reduction identities to reduce other arguments to this interval.
                double y = x;
                int n = 0;
                bool arg_was_less_than_one = ( y < 1.0 );
                // Add or subtract integers as necessary to bring y into (1,2)
                // Will correct for this below
                if ( arg_was_less_than_one ) {
                        y += 1.0;
                } else {
                        n = static_cast<int> ( floor ( y ) ) - 1;  // will use n later
                        y -= n;
                }
                // numerator coefficients for approximation over the interval (1,2)
                static const double p[] = {
                        -1.71618513886549492533811E+0,
                        2.47656508055759199108314E+1,
                        -3.79804256470945635097577E+2,
                        6.29331155312818442661052E+2,
                        8.66966202790413211295064E+2,
                        -3.14512729688483675254357E+4,
                        -3.61444134186911729807069E+4,
                        6.64561438202405440627855E+4
                };
                // denominator coefficients for approximation over the interval (1,2)
                static const double q[] = {
                        -3.08402300119738975254353E+1,
                        3.15350626979604161529144E+2,
                        -1.01515636749021914166146E+3,
                        -3.10777167157231109440444E+3,
                        2.25381184209801510330112E+4,
                        4.75584627752788110767815E+3,
                        -1.34659959864969306392456E+5,
                        -1.15132259675553483497211E+5
                };
                double num = 0.0;
                double den = 1.0;
                int i;
                double z = y - 1;
                for ( i = 0; i < 8; i++ ) {
                        num = ( num + p[i] ) * z;
                        den = den*z + q[i];
                }
                double result = num/den + 1.0;
                // Apply correction if argument was not initially in (1,2)
                if ( arg_was_less_than_one ) {
                        // Use identity gamma(z) = gamma(z+1)/z
                        // The variable "result" now holds gamma of the original y + 1
                        // Thus we use y-1 to get back the orginal y.
                        result /= ( y - 1.0 );
                } else {
                        // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
                        for ( i = 0; i < n; i++ ) result *= y++;
                }
                return result;
        }

////////////////////////////////////////////////////////////////////////////////////////////
// Third interval: [12, infinity)
        if ( x > 171.624 ) {
                // Correct answer too large to display. Force +infinity.
                double temp = DBL_MAX;
                return temp * 2.0;
        }
        // Abramowitz and Stegun 6.1.41
        // Asymptotic series should be good to at least 11 or 12 figures
        // For error analysis, see Whittiker and Watson
        // A Course in Modern Analysis (1927), page 252
        static const double c[8] = { 1.0/12.0, -1.0/360.0, 1.0/1260.0, -1.0/1680.0, 1.0/1188.0, -691.0/360360.0, 1.0/156.0, -3617.0/122400.0 };
        double z = 1.0 / ( x * x );
        double sum = c[7];
        for ( int i = 6; i >= 0; i-- ) {
                sum *= z;
                sum += c[i];
        }
        double series = sum/x;
        static const double halfLogTwoPi = 0.91893853320467274178032973640562;
        double logGamma = ( x - 0.5 ) * log ( x ) - x + halfLogTwoPi + series;
        return exp ( logGamma );
} //////////////////////////////////////////////////////////////////////////////////////////
} // namespace

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
    MFPFree / BoxHeight;


std::cout << "\n";
std::cout << "Freestream properties\n";
std::cout << "-----------------------------\n";
std::cout << "Box size    = " << BoxWidth << " x " << BoxHeight << " m\n";
std::cout << "Angle       = " << AttackAngle << " degrees\n";
std::cout << "P_inf       = " << PFree << " Pa\n";
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

    RemoveParticlesInsideBox();
}
