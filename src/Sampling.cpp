#include "Sampling.h"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "SimulationState.h"
#include "Boundary.h"

#define sqr3( V ) ( V[0] * V[0] + V[1] * V[1] + V[2] * V[2] )

void Sampling () ////////////////////////////////////////////////////////////////////////////
        {

                SampleStep++;
                SampleTime += Dt;
                for ( int k = 0; k < LocalNX; k++ )
                        for ( int l = 0; l < LocalNY; l++ )
                                for ( int i = 0; i < NPC[k][l]; i++ ) {
                                        C[k][l].CountNP += 1.0;
                                        double V2 = sqr3 ( P[IPC[k][l][i]].V );
                                        for ( int m = 0; m < 3; m++ ) {
                                                C[k][l].CountV[m] += P[IPC[k][l][i]].V[m];
                                                C[k][l].CountV2V[m] += V2 * P[IPC[k][l][i]].V[m];
                                                for ( int n = 0; n < 3; n++ )
                                                        C[k][l].CountVV[m][n] += P[IPC[k][l][i]].V[m] * P[IPC[k][l][i]].V[n];
                                        }
                        }
        }




void Printing(
    const std::vector<CELL> &GlobalField
)
{
    if (MPIRank != 0)
    {
        return;
    }


    if (
        GlobalField.size()
        != static_cast<size_t>(NX * NY)
    )
    {
        std::cerr
            << "ERROR: GlobalField has wrong size\n"
            << "Expected = "
            << NX * NY
            << '\n'
            << "Received = "
            << GlobalField.size()
            << '\n';

        return;
    }


    const std::filesystem::path outputDirectory = DSMC_RESULTS_DIRECTORY;
    const std::filesystem::path outputFile = outputDirectory / "DSMC2D_Field.plt";

    std::error_code directoryError;
    std::filesystem::create_directories(outputDirectory, directoryError);
    if (directoryError)
    {
        std::cerr << "ERROR: Could not create results directory: "
                  << directoryError.message() << '\n';
        return;
    }

    FILE *F =
        fopen(
            outputFile.string().c_str(),
            "wt"
        );


    if (F == nullptr)
    {
        std::cerr
            << "ERROR: Could not open "
            << outputFile.string() << '\n';

        return;
    }


    fprintf(
        F,
        "VARIABLES = "
        "\"X\" \"Y\" \"N\" "
        "\"Ux\" \"Uy\" "
        "\"Tx\" \"Ty\" \"T\"\n"
    );

    fprintf(
        F,
        "ZONE I=%d J=%d\n",
        NX,
        NY
    );


    // =====================================================
    // Write complete global field
    // =====================================================

    for (int k = 0; k < NX; k++)
    {
        for (int l = 0; l < NY; l++)
        {
            const CELL &cell =
                GlobalField[
                    k * NY + l
                ];


            double X =
                X1 + (k + 0.5) * DX;

            double Y =
                Y1 + (l + 0.5) * DY;


            double N = 0.0;

            double U[3] =
            {
                0.0,
                0.0,
                0.0
            };

            double Txyz[3] =
            {
                0.0,
                0.0,
                0.0
            };


            // =================================================
            // Number density
            // =================================================

            if (SampleStep > 0)
            {
                N =
                    Weight
                    * cell.CountNP
                    / SampleStep
                    / DV;
            }


            // =================================================
            // Velocity and temperature
            // =================================================

            if (cell.CountNP > 0.0)
            {
                for (int m = 0; m < 3; m++)
                {
                    U[m] =
                        cell.CountV[m]
                        / cell.CountNP;


                    Txyz[m] =
                        MoleculeMass
                        *
                        (
                            cell.CountVV[m][m]
                            / cell.CountNP

                            - U[m] * U[m]
                        )
                        /
                        BOLTZMANN_CONSTANT;
                }
            }


            double T =
                (
                    Txyz[0]
                    + Txyz[1]
                    + Txyz[2]
                )
                / 3.0;


            fprintf(
                F,

                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e "
                "%18.10e\n",

                X,
                Y,
                N,
                U[0],
                U[1],
                Txyz[0],
                Txyz[1],
                T
            );
        }
    }


    fclose(F);

    // Keep the plotting geometry with the field rather than duplicating it
    // in Python settings. The field's two-line header remains unchanged.
    const auto metadataFile = outputDirectory / "DSMC2D_Field.json";
    std::ofstream metadata(metadataFile);
    metadata << std::setprecision(17)
             << "{\n  \"obstacle\": {\"type\": \"box\", \"center\": ["
             << BoxX << ", " << BoxY << "], \"width\": " << BoxWidth
             << ", \"height\": " << BoxHeight << "},\n"
             << "  \"angle_of_attack_degrees\": " << AttackAngle << ",\n"
             << "  \"wall_model\": \"" << WallModelName() << "\"\n}\n";
    metadata.close();
    if (!metadata)
        std::cerr << "ERROR: Could not write " << metadataFile.string() << '\n';

    std::cout
        << "Global field written to "
        << outputFile.string() << '\n';
}
