#pragma once

// Accumulated moments used to calculate macroscopic field properties.
typedef struct cell {
    double CountNP;          // Particle count
    double CountV[3];        // Velocity sum
    double CountVV[3][3];    // Velocity tensor sum
    double CountV2V[3];      // Squared speed times velocity sum
    double CountNC;          // Collision count
} CELL;
