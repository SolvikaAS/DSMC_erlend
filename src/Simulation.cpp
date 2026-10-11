#include "Simulation.h"
#include <iostream>
#include <ctime>

#include "Initialization.h"
#include "SimulationState.h"
#include "Boundary.h"
#include "Geometry.h"
#include "Collision.h"
#include "Sampling.h"
#include "MPIManager.h"
#include "RNG.h"

int RunSimulation(int argc, char **argv)
{

        MPI_Init(&argc,&argv);

        CreateMPICellType();
        CreateMPIParticleType();

        MPI_Comm_rank(MPI_COMM_WORLD, &MPIRank);
        MPI_Comm_size(MPI_COMM_WORLD, &MPISize);

    // Give each MPI process its own RNG seed
    time_t t;
    unsigned baseSeed = static_cast<unsigned>(time(&t));

    SetSeed(
        baseSeed + 7919u * MPIRank,
        362436069u + 104729u * MPIRank
    );

	Setup();
	std::cout<<"Setup \n";
	SetupMPIDomain();

	/*
	std::cout
        << "Rank " << MPIRank
        << " / " << MPISize
        << " : cells "
        << GlobalStartCell
        << " to "
        << GlobalEndCell - 1
        << ", LocalNX = "
        << LocalNX
        << ", x = "
        << LocalX1
        << " to "
        << LocalX2
        << '\n';
	*/

	InitialConditions();
	std::cout
    << "Rank " << MPIRank
    << " initial particles = "
    << NP
    << '\n';

    int GlobalNP = 0;

    MPI_Reduce(
        &NP,
        &GlobalNP,
        1,
        MPI_INT,
        MPI_SUM,
        0,
        MPI_COMM_WORLD
    );

    if (MPIRank == 0)
    {
        std::cout
            << "Total particles across all ranks = "
            << GlobalNP
            << '\n';
    }


	std::cout<<"InitialConditions \n";
	do {
		if (Step > FirstSamplingStep && Step % SamplingPeriod ==0)
			SamplingMode =1;
		else
			SamplingMode =0;
		MoveParticles();
		if(Step>FirstSamplingStep){
			DragSamplingTime +=Dt;
		}
		//std::cout<<"moved \n";
		BoundaryConditions();
		//std::cout<<"boundary \n";

		ExchangeParticles();

		Indexing();
		//std::cout<<"indexed \n";
		int NC = CollideParticles();
		if (SamplingMode ==1)
		{
		 Sampling();
		}
		Step++;
		Time += Dt;
		if (Step % 100 == 0)
                {
                    int GlobalNP = 0;
                    int GlobalNC = 0;

                    MPI_Reduce(
                        &NP,
                        &GlobalNP,
                        1,
                        MPI_INT,
                        MPI_SUM,
                        0,
                        MPI_COMM_WORLD
                    );

                    MPI_Reduce(
                        &NC,
                        &GlobalNC,
                        1,
                        MPI_INT,
                        MPI_SUM,
                        0,
                        MPI_COMM_WORLD
                    );

                    if (MPIRank == 0)
                    {
                        std::cout
                            << "Time step " << Step
                            << ", Number of particles "
                            << GlobalNP
                            << ", Number of collisions "
                            << GlobalNC
                            << '\n';
                    }
                }
		//if (Step > FirstSamplingStep && Step % PrintPeriod ==0) Printing();
	} while (Step < NStep);

	std::vector<CELL> FinalField =
        GatherGlobalField();


	if (MPIRank== 0)
      {
        Printing(FinalField);
      }

/*

	double MeanDrag =
    DragImpulseTotal / DragSamplingTime;

double RhoFree = NFree * MoleculeMass;
double AreaRef = 2.0 * SphereRadius * DZ;

double Cd =
    MeanDrag /
    (0.5 * RhoFree * UFree * UFree * AreaRef);

std::cout << "\n";
std::cout << "====================================\n";
std::cout << "Sphere collisions = "
          << DragCollisionCount << "\n";

std::cout << "Sampling time = "
          << DragSamplingTime << " s\n";

std::cout << "Mean drag force = "
          << MeanDrag << " N\n";

std::cout << "Drag coefficient Cd = "
          << Cd << "\n";

std::cout << "====================================\n";
*/

MPI_Type_free(&MPI_PCL);
MPI_Type_free(&MPI_CELL);

MPI_Comm_free(&CartesianComm);

MPI_Finalize();

return 0;
}
