#pragma once

// Binary collision selection and elastic scattering.
void ElasticCollision(double *V1, double *V2, double Cr);
int BinaryCollisionsInCell(double CellVolume, int NPC, int *IPC);
int CollideParticles();
