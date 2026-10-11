#pragma once

// Domain and obstacle settings. Definitions live in Geometry.cpp.
extern double X1, X2;
extern double Y1, Y2;
extern double DZ;
extern int NX, NY;
extern double DL;
extern double SphereX, SphereY, SphereRadius;
extern double DX, DY;
extern double DV;

// Calculate cell dimensions, index particles, and exclude the obstacle interior.
void SetupGeometry();
void Indexing();
bool InsideSphere(double x, double y);
void RemoveParticlesInsideSphere();
