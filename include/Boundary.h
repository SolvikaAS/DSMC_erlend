#pragma once

enum class WallModel { Diffuse, Specular };
extern WallModel WallInteraction;
const char *WallModelName();

// Particle injection, physical boundaries, and motion with wall scattering.
void GenerateNewParticles(double X1, double Y1, double X2, double Y2,
                          double Ngas, double Ugasx, double Ugasy,
                          double Tgas, int MoveFlag);
void BoundaryConditions();
void DiffuseScattering(double *v, double m, double Tw, double Nx, double Ny);
// Stationary wall with an outward unit normal: reverse normal velocity only.
void SpecularScattering(double *v, double Nx, double Ny);
void MoveParticles();
