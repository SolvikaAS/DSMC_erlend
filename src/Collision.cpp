#include "Collision.h"
#include <cmath>

#include "SimulationState.h"
#include "RNG.h"

void ElasticCollision ( double *V1, double *V2, double Cr){
	double N[3];
	v3rand_isotropic (N);
	Cr *=0.5;
	double VC[3] = {0.5*(V2[0]+V1[0]),0.5*(V2[1]+V1[1]),0.5*(V2[2]+V1[2])};
	double VCr[3] = {Cr*N[0],Cr*N[1],Cr*N[2]};
	V1[0]=VC[0]+VCr[0];
	V1[1]=VC[1]+VCr[1];
	V1[2]=VC[2]+VCr[2];
	V2[0]=VC[0]-VCr[0];
	V2[1]=VC[1]-VCr[1];
	V2[2]=VC[2]-VCr[2];
}

int BinaryCollisionsInCell (double CellVolume,int NPC,int *IPC){
	if (NPC<2)return 0;
	double SCMax = 9.0*SigmaRef*sqrt(BOLTZMANN_CONSTANT*TFree/MoleculeMass);
	double Npair = 0.5*NPC*(NPC-1)*Weight*SCMax*Dt/CellVolume;
	int NN=int (Npair);
	double N1= Npair-NN;
	NN =(brng()<N1)?NN+1:NN;
	int i,j;
	int NC=0;
	for (int k=0; k<NN;k++){
		i = irand_uniform(NPC);
		do{j = irand_uniform(NPC);}while (j==i);
		i = IPC[i];
		j = IPC[j];
		double VCr[3] = {P[j].V[0]-P[i].V[0],P[j].V[1]-P[i].V[1],P[j].V[2]-P[i].V[2]};
		double Cr = sqrt(VCr[0]*VCr[0]+VCr[1]*VCr[1]+VCr[2]*VCr[2]);
		double Sigma = SigmaRef*pow(CrRef/Cr, Omega);
		if ( brng()<Sigma*Cr/SCMax){
			ElasticCollision (P[i].V,P[j].V,Cr);
			NC++;
		}
	}
	return NC;
}
int CollideParticles () ////////////////////////////////////////////////////////////////////
        {
        int     NC = 0; // Counter of collisions
                for ( int k = 0; k < LocalNX; k++ )
                        for ( int l = 0; l < LocalNY; l++ ) { // Sampling of collisions in particle (i,j)
                                int NCC = BinaryCollisionsInCell ( DV, NPC[k][l], IPC[k][l] ); // NCC is the number of collisions sampled in the current cell
                                NC += NCC;
                                // Here we update CountNC in the cell for fure calculations of the mean free time and path, if necessary
                                if ( SamplingMode == 1 ) C[k][l].CountNC += NCC;
                        }
                return NC; // NC is the total number of real collsions sampled during the time step
}
