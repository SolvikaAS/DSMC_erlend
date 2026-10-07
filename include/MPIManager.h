#pragma once

#include <vector>
#include "Cell.h"

// Domain decomposition, MPI data types, particle exchange, and field gathering.
void CreateMPIParticleType();
void CreateMPICellType();
void SplitDimension(int GlobalN, int Coordinate, int NProcesses,
                    int &LocalN, int &StartCell, int &EndCell);
void SetupMPIDomain();
int OwnerRankFromGlobalCell(int globalK);
int OwnerCoordinate(int GlobalCell, int GlobalN, int NProcesses);
int ParticleOwnerRank(double x, double y);
void ExchangeParticles();
std::vector<CELL> GatherGlobalField();
