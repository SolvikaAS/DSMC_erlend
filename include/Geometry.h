#pragma once

// Cell indexing and exclusion of particles inside the obstacle.
void Indexing();
bool InsideSphere(double x, double y);
void RemoveParticlesInsideSphere();
