#include "Boundary.h"
#include "Geometry.h"
#include "Initialization.h"
#include "RNG.h"
#include "SimulationState.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {
void Check(bool condition, const char *message)
{
    if (!condition)
    {
        std::cerr << message << '\n';
        std::exit(1);
    }
}

void Hit(double x, double y, double dx, double dy,
         double expectedTime, double expectedNX, double expectedNY)
{
    double time, nx, ny;
    Check(IntersectBox(x, y, dx, dy, time, nx, ny), "Expected a box impact");
    Check(std::abs(time - expectedTime) < 1e-12, "Incorrect impact time");
    Check(nx == expectedNX && ny == expectedNY, "Incorrect outward normal");
}

void Miss(double x, double y, double dx, double dy)
{
    double time, nx, ny;
    Check(!IntersectBox(x, y, dx, dy, time, nx, ny), "Unexpected box impact");
}
}

int main()
{
    BoxX = 0.0;
    BoxY = 0.0;
    BoxWidth = 6.0;
    BoxHeight = 6.0;
    Check(InsideBox(0.0, 0.0), "Box centre must be inside");
    Check(InsideBox(2.9, -2.9), "Interior near corner must be inside");
    Check(!InsideBox(3.0, 0.0), "Surface must not count as interior");
    Check(!InsideBox(0.0, 3.1), "Exterior must be outside");

    Hit(-4, 0, 2, 0, 0.5, -1, 0);
    Hit(4, 0, -2, 0, 0.5, 1, 0);
    Hit(0, -4, 0, 2, 0.5, 0, -1);
    Hit(0, 4, 0, -2, 0.5, 0, 1);
    Hit(-10, 0, 20, 0, 0.35, -1, 0); // Ends beyond the far face.
    Hit(-4, 0, 1, 0, 1.0, -1, 0);    // Impact at end of step.
    Hit(-3, 0, 1, 0, 0.0, -1, 0);    // Starts on surface, heads inward.
    Hit(-4, -4, 2, 2, 0.5, -1, 0);   // Exact corner entry.
    Hit(-4, 2, 4, -1, 0.25, -1, 0);  // Oblique impact.
    Miss(-4, 0, -1, 0);
    Miss(-3, 0, -1, 0);
    Miss(-4, 3, 8, 0);               // Along a face.
    Miss(-4, 2, 2, 2);               // Touches a corner without entering.
    Miss(4, 4, 0, 0);
    Miss(-4, 0, 0.5, 0);

    BoxX = 2.0;
    BoxY = -1.0;
    BoxWidth = 4.0;
    BoxHeight = 2.0;
    Hit(-1, -1, 2, 0, 0.5, -1, 0);
    Hit(2, 1, 0, -2, 0.5, 0, 1);
    Check(InsideBox(3.9, -1.9), "Translated rectangular box interior");
    BoxX = BoxY = 0.0;
    BoxWidth = BoxHeight = 6.0;

    AttackAngle = 90.0;
    Setup();
    Check(std::abs(UxFree) < 1e-10, "90 degree flow must have zero x velocity");
    Check(std::abs(UyFree + UFree) < 1e-10, "Positive angle must tilt toward negative y");
    AttackAngle = -45.0;
    Setup();
    Check(UxFree > 0 && UyFree > 0, "Negative angle must tilt toward positive y");

    NP = 3;
    P[0] = {{0, 0}, {0, 0, 0}};
    P[1] = {{4, 0}, {0, 0, 0}};
    P[2] = {{2, 2}, {0, 0, 0}};
    RemoveParticlesInsideBox();
    Check(NP == 1 && P[0].X[0] == 4, "Initial interior particles must be removed");

    SetSeed(12345, 67890);
    Dt = 0.002;
    const PCL approaches[] = {
        {{-4, 0}, {1000, 0, 0}}, {{4, 0}, {-1000, 0, 0}},
        {{0, -4}, {0, 1000, 0}}, {{0, 4}, {0, -1000, 0}},
        {{-4, -4}, {1000, 1000, 0}},
    };
    for (const auto &approach : approaches)
    {
        for (int sample = 0; sample < 100; ++sample)
        {
            NP = 1;
            P[0] = approach;
            double time, nx, ny;
            Check(IntersectBox(P[0].X[0], P[0].X[1], Dt * P[0].V[0],
                               Dt * P[0].V[1], time, nx, ny), "Reflection setup");
            DragCollisionCount = 0;
            MoveParticles();
            Check(DragCollisionCount == 1, "Impact must be counted once");
            Check(!InsideBox(P[0].X[0], P[0].X[1]), "Reflected particle entered box");
            Check(P[0].V[0] * nx + P[0].V[1] * ny > 0, "Reflection points inward");
            Check(std::isfinite(P[0].X[0]) && std::isfinite(P[0].X[1]), "Invalid position");
        }
    }
    NP = 1;
    P[0] = {{-10, 10}, {1000, 0, 0}};
    DragCollisionCount = 0;
    MoveParticles();
    Check(P[0].X[0] == -8 && P[0].X[1] == 10 && P[0].V[0] == 1000,
          "Free flight changed unexpectedly");
    Check(DragCollisionCount == 0, "Free flight counted as impact");
    std::cout << "Box intersection, initialization, flow angles, and reflection checks passed.\n";
}
