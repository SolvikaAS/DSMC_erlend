#pragma once

#include <mpi.h>
#include "Particle.h"
#include "Cell.h"
#include "Geometry.h"

// Shared state for the existing procedural solver. Geometry is declared in Geometry.h;
// the remaining definitions live in Initialization.cpp.
extern double Dt;
extern int NStep;
extern int FirstSamplingStep;
extern int SamplingPeriod;
extern int PrintPeriod;
extern double NPCFree;
extern double MolarMass;
extern double MaFree;
extern double NFree;
extern double TFree;
extern double AttackAngle;
extern double Tw;
extern double Omega;
extern double SigmaRef;
extern double CrRef;
extern double MFTFree;
extern double MFPFree;
extern double KnFree;
extern double DRef;
extern double ViscosityIndex;
extern double Tref;
extern double DragForce;
extern double DragForceSum;
extern int DragSamples;
extern double DragImpulse;
extern long long DragCollisionCount;
extern double DragImpulseTotal;
extern MPI_Datatype MPI_CELL;
extern MPI_Datatype MPI_PCL;
extern int MPIRank;
extern int MPISize;
extern MPI_Comm CartesianComm;
extern int ProcX;
extern int ProcY;
extern int RankX;
extern int RankY;
extern int LeftRank;
extern int RightRank;
extern int BottomRank;
extern int TopRank;
extern int LocalNX;
extern int LocalNY;
extern int GlobalStartCellX;
extern int GlobalEndCellX;
extern int GlobalStartCellY;
extern int GlobalEndCellY;
extern double LocalX1;
extern double LocalX2;
extern double LocalY1;
extern double LocalY2;
extern int SampleStep;
extern CELL C[MAX_CELL_X][MAX_CELL_Y];
extern int NP;
extern PCL P[MAX_PCL];
extern double Weight;
extern double MoleculeMass;
extern double PFree;
extern double UFree;
extern double UxFree, UyFree;
extern int NPC[MAX_CELL_X][MAX_CELL_Y];
extern int IPC[MAX_CELL_X][MAX_CELL_Y][MAX_PCL_IN_CELL];
extern double SampleTime;
extern int Step;
extern double Time;
extern int SamplingMode;
extern double DragSamplingStart;
extern double DragSamplingTime;
