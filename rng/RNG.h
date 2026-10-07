#pragma once

void SetSeed(unsigned int u, unsigned int v);
unsigned int IntUniform(unsigned int &u, unsigned int &v);
double brng();
int irand_uniform(int N);
int irand_Poisson(double E);
double frand_Rayleigh(double S);
void v3rand_isotropic(double *R);
