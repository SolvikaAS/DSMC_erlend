#include "Boundary.h"
#include "Geometry.h"
#include "Initialization.h"
#include "RNG.h"
#include "SimulationState.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void Check(bool condition, const char *message)
{
    if (!condition)
    {
        std::cerr << message << '\n';
        std::exit(1);
    }
}

void CheckMoments(double alphaN, double alphaT)
{
    // Set kB*Tw/m = 1 and sample in a rotated frame. Incident components
    // are (-2, 3, 1.5). Check analytic conditional first/second moments.
    const double n = std::sqrt(0.5);
    const double mass = BOLTZMANN_CONSTANT * 300.0;
    constexpr int samples = 100000;
    double sumNormal2 = 0, sumTangent = 0, sumTangent2 = 0, sumZ = 0, sumZ2 = 0;
    SetSeed(12345, 67890);
    for (int i = 0; i < samples; ++i)
    {
        double v[] = {-5*n, n, 1.5};
        CLLScattering(v, mass, 300, n, n, alphaN, alphaT);
        const double normal = (v[0] + v[1]) * n;
        const double tangent = (-v[0] + v[1]) * n;
        Check(normal > 0 && std::isfinite(normal), "CLL must emit a finite outward velocity");
        sumNormal2 += normal * normal;
        sumTangent += tangent;
        sumTangent2 += tangent * tangent;
        sumZ += v[2];
        sumZ2 += v[2] * v[2];
    }
    const double normalMean2 = (1-alphaN)*4 + 2*alphaN;
    const double normalVariance2 = 4*alphaN*(1-alphaN)*4 + 4*alphaN*alphaN;
    const double tangentMean = std::sqrt(1-alphaT)*3;
    const double zMean = std::sqrt(1-alphaT)*1.5;
    const double meanTolerance = 6*std::sqrt(alphaT/samples) + 1e-10;
    const double tangentSecondTolerance =
        6*std::sqrt((2*alphaT*alphaT + 4*tangentMean*tangentMean*alphaT)/samples) + 1e-10;
    const double zSecondTolerance =
        6*std::sqrt((2*alphaT*alphaT + 4*zMean*zMean*alphaT)/samples) + 1e-10;
    Check(std::abs(sumNormal2/samples - normalMean2)
          < 6*std::sqrt(normalVariance2/samples) + 1e-10, "Incorrect CLL normal energy moment");
    Check(std::abs(sumTangent/samples - tangentMean) < meanTolerance, "Incorrect CLL tangential mean");
    Check(std::abs(sumZ/samples - zMean) < meanTolerance, "Incorrect CLL z mean");
    Check(std::abs(sumTangent2/samples - tangentMean*tangentMean - alphaT)
          < tangentSecondTolerance, "Incorrect CLL tangential energy moment");
    Check(std::abs(sumZ2/samples - zMean*zMean - alphaT)
          < zSecondTolerance, "Incorrect CLL z energy moment");
}
}

int main()
{
    // The zero-accommodation limit is exactly specular and uses no RNG draws.
    double specular[] = {-3, -4, 5};
    double cll[] = {-3, -4, 5};
    const double diagonal = std::sqrt(0.5);
    SpecularScattering(specular, diagonal, diagonal);
    SetSeed(111, 222);
    const double nextRandom = brng();
    SetSeed(111, 222);
    CLLScattering(cll, BOLTZMANN_CONSTANT*300, 300, diagonal, diagonal, 0, 0);
    for (int m = 0; m < 3; ++m)
        Check(cll[m] == specular[m], "CLL (0,0) must be exactly specular");
    Check(brng() == nextRandom, "CLL specular limit should not consume randomness");

    CheckMoments(1, 1);     // Flux-weighted diffuse limit.
    CheckMoments(0.3, 0.7); // Independent partial normal/tangential accommodation.
    CheckMoments(0, 0.6);   // Specular normal, partially thermalized tangential velocities.
    CheckMoments(0.6, 0);   // Incident tangential components retained exactly.

    const double invalid[] = {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()};
    for (double coefficient : invalid)
    {
        for (bool normalInvalid : {false, true})
        {
            bool rejected = false;
            double v[] = {-2, 3, 1.5};
            try
            {
                CLLScattering(v, BOLTZMANN_CONSTANT*300, 300, 1, 0,
                              normalInvalid ? coefficient : 0.8,
                              normalInvalid ? 0.8 : coefficient);
            }
            catch (const std::invalid_argument &) { rejected = true; }
            Check(rejected, "Invalid CLL coefficient was accepted");
        }
    }

    // Exercise actual box impacts on every face and at the corner.
    Setup();
    BoxX = BoxY = 0;
    BoxWidth = BoxHeight = 6;
    WallInteraction = WallModel::CLL;
    Check(std::string(WallModelName()) == "cll", "Wrong CLL result folder name");
    CLLNormalAccommodation = 0.3;
    CLLTangentialAccommodation = 0.7;
    Dt = 0.002;
    const PCL approaches[] = {
        {{-4, 0}, {1000, 100, 50}}, {{4, 0}, {-1000, 100, 50}},
        {{0, -4}, {100, 1000, 50}}, {{0, 4}, {100, -1000, 50}},
        {{-4, -4}, {1000, 1000, 50}},
    };
    for (const auto &approach : approaches)
    {
        for (int sample = 0; sample < 100; ++sample)
        {
            NP = 1;
            P[0] = approach;
            double time, nx, ny;
            Check(IntersectBox(P[0].X[0], P[0].X[1], Dt*P[0].V[0], Dt*P[0].V[1],
                               time, nx, ny), "Test particle must hit box");
            DragCollisionCount = 0;
            MoveParticles();
            Check(DragCollisionCount == 1, "CLL impact count is incorrect");
            Check(!InsideBox(P[0].X[0], P[0].X[1]), "CLL particle entered the box");
            Check(P[0].V[0]*nx + P[0].V[1]*ny > 0, "CLL reflection points inward");
        }
    }
    std::cout << "CLL limits, conditional moments, coefficient validation, and box impacts passed.\n";
}
