#pragma once

#include "Constants.h"

// Two position components and three velocity components per particle.
typedef struct pcl {
    double X[DIM];
    double V[3];
} PCL;
