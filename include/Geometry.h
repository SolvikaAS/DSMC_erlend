#pragma once

// Domain and obstacle settings. Definitions live in Geometry.cpp.
extern double X1, X2;
extern double Y1, Y2;
extern double DZ;
extern int NX, NY;
extern double DL;
extern double BoxX, BoxY, BoxWidth, BoxHeight;
extern double DX, DY;
extern double DV;

// Calculate cell dimensions, index particles, and exclude the obstacle interior.
void SetupGeometry();
void Indexing();
bool InsideBox(double x, double y);
void RemoveParticlesInsideBox();

// First entry into the box along a segment; tHit is in [0, 1].
// The returned unit normal points out of the impacted face.
bool IntersectBox(double x, double y, double dx, double dy,
                  double &tHit, double &normalX, double &normalY);
