// Form streamlines from saddle points and seed points
// --- Input: 1. 3D vector field
//            2. Seed points
// --- Output: skelelton file
// --- Author: Xiaosong Yuan, Bala, Vizlab, Rutgers University
// --- Date: 9/28/2002
// --- Original file name: skel_streamline_float.cpp
//
// --- Changed to StreamLn.cpp by Nicu D. Cornea, on Tuesday, July 22, 2003
//	What it does:
//		Is finds the critical points by dividing each cell into a 5x5x5 grid and tries to find
//			force vectors that are smaller than a certain threashold.
//			The problem is how we chose that threshold.
//

#define TRACE

#include "StreamLn.h"

#define MAX_NUM_SKEL	30000
#define MAX_NUM_SEEDS	50000
#define MAX_NUM_CRITPTS	100

inline double veclength(Vector vecin) {
    return sqrt(vecin.xd * vecin.xd + vecin.yd * vecin.yd +vecin.zd * vecin.zd);
}
inline Vector interpolation(float x, float y, float z, int sizx, int sizy, int sizz, Vector *forcevec);
inline void rk2(float x, float y, float z, int sizx, int sizy, int sizz, double steps, Vector *Force_ini, VoxelPositionDouble *nextPos);


bool GetStreamLines(
	Vector* ForceField, int L, int M, int N,
	unsigned char *flags,
	VoxelPosition *BdSeeds, int numBoundSeeds,
	VoxelPositionDouble **Skel, int *numSkel)
{

#ifdef TRACE
	printf("TRACE: Starting GetStreamLines function...\n");
#endif

	printf("***************!!!!!!! StreamLn version 1 - don't use this !!!!!!!!!!!");

	// Vector *force;
	Vector *Normforce;
	long idx, iidx, slsz;

	VoxelPositionDouble Startpos, Nextpos;
	VoxelPositionDouble seeds[MAX_NUM_SEEDS];
	VoxelPositionDouble critPts[MAX_NUM_CRITPTS];
	int idxSeeds, nrCritPts;

	Vector vecin;

	int cc;
	int ii,jj,kk;

	int *FlagOnSkeleton;
	float vecLength;
	int NumCritPoints=0;
	float Ngrid;

	Vector OutForce;
	int streamSteps=0;

	int numSeeds;
	int i,j,k;

	double thCritPt = 0.00;
	bool critPtFound;


	slsz = L*M;		// slice size


	// force = new Vector[L*M*N];
	Normforce = new Vector[L*M*N];
	FlagOnSkeleton = new int[L*M*N];
	if(	(FlagOnSkeleton == NULL)	||
		(Normforce == NULL))	{
		printf("UPS! - Error allocating memory for working data structures. Abort.\n");
		exit(1);
	}


	(*Skel) = NULL;
	(*numSkel) = 0;

	if(((*Skel) = new VoxelPositionDouble[MAX_NUM_SKEL]) == NULL) {
		printf("UPS! - Error allocating memory for the output array. Abort.\n");
		exit(1);
	}


#ifdef _DEBUG
	PrintElapsedTime("Phase4.1: Allocating memory");
#endif

	// analyze vector field and find threshold value for vector length
	for (k = 0; k < N; k++) {
		for (j = 0; j < M; j++) {
			for (i = 0; i < L; i++) {
				idx = k*slsz + j*L +i;

				vecLength = sqrt(ForceField[idx].xd * ForceField[idx].xd +
							ForceField[idx].yd * ForceField[idx].yd +
							ForceField[idx].zd * ForceField[idx].zd);
				if(vecLength == 0.00) {
					vecLength = 1.00;
				}
				// normalize the vectors
				Normforce[idx].xd = ForceField[idx].xd / vecLength;
				Normforce[idx].yd = ForceField[idx].yd / vecLength;
				Normforce[idx].zd = ForceField[idx].zd / vecLength;

				FlagOnSkeleton[idx] = 0;
			}
		}
	}
/*
	// write out the normalized force vectors
	for (k = 0; k < N; k++) {
		for (j = 0; j < M; j++) {
			for (i = 0; i < L; i++) {
				printf("%f %f %f\n",Normforce[k*slsz+j*L+i].xd,Normforce[k*slsz+j*L+i].yd, Normforce[k*slsz+j*L+i].zd);
			}
		}
	}
*/
#ifdef _DEBUG
	PrintElapsedTime("Phase4.2: Fill in the force and Normforce arrays");
#endif


	// copy boundary seeds to the seeds array
	numSeeds = numBoundSeeds;
	for (i = 0; i< numSeeds; i++)  {
		seeds[i].x = BdSeeds[i].x;
		seeds[i].y = BdSeeds[i].y;
		seeds[i].z = BdSeeds[i].z;
	}


#ifdef _DEBUG
	PrintElapsedTime("Phase4.4: Copy boundary seeds");
#endif

////////////////////////////////
	// find the max and min vector lenghts:
	double maxvl, minvl, vl;

	maxvl = -0.00;
	minvl = 9999.99;
	for(idx=0; idx < L*M*N; idx++) {
		if(flags[idx] == INTERIOR) {
			vl = 	sqrt((ForceField[idx].xd * ForceField[idx].xd) + (ForceField[idx].yd * ForceField[idx].yd) + (ForceField[idx].zd * ForceField[idx].zd));

			if(vl > maxvl) {
				maxvl = vl;
			}

			if(vl < minvl) {
				minvl = vl;
			}
		}
	}

	// calculate threshold as minvalue + 0.0001% of the difference between maxvl and minvl
	// thCritPt = (minvl + ((maxvl - minvl) * 0.00001 / 100.00));
	thCritPt = 0.1;

#ifdef _DEBUG
	printf("Vector field analysis:\n");
	printf("\tmaximum vector length: %f\n", maxvl);
	printf("\tminimum vector length: %f\n", minvl);
	printf("\tCritical point threshold: %f\n", thCritPt);
#endif

/////////////////////////


	// find all critical points -- method 4: use points with small length of vector
	//
	// see other methods in original file
	//
	Ngrid = 5;
	NumCritPoints = 0;

	for (k = 0; k < N-1; k++) {

		printf("\tProcessing plane %d out of %d\r", k, N-1);
		fflush(stdout);

		for (j = 0; j < M-1; j++) {
			for (i = 0; i < L-1; i++) {
#ifdef TRACE
			//printf("\t%d, %d, %d\n", k, j, i);
#endif
				idx = k*slsz + j*L +i;

				// - if this point one of it's neighbors is external, ignore this point,
				//	since the interpolation might give wrong results if one of the neighbors
				//	has a force = 0 because it's an external point not because it is a critical point
				// - face neighbors are covered by the flags.
				//	If a face neighbor is external, then this point is a SURF or BOUNDARY point
				if ((flags[idx] == EXTERIOR)) continue;
				if ((flags[idx] == SURF) || (flags[idx] == BOUNDARY)) continue;
				// - further, we have to check a few edge neighbors
				//		and a vertex neighbor that will be touched by the interpolation
				if(flags[idx + L + 1] == EXTERIOR) continue;		//(i+1, j+1, k)
				if(flags[idx + slsz + 1] == EXTERIOR) continue;		//(i+1, j, k+1)
				if(flags[idx + slsz + L] == EXTERIOR) continue;		//(i, j+1, k+1)
				if(flags[idx + slsz + L + 1] == EXTERIOR) continue;	//(i+1, j+1, k+1)

				critPtFound = false;

				for (kk = 0; kk < Ngrid && !critPtFound; kk++) {
					for (jj = 0; jj < Ngrid && !critPtFound; jj++) {
						for (ii = 0; ii < Ngrid && !critPtFound; ii++) {
							OutForce=interpolation(i+ii/Ngrid, j+jj/Ngrid, k+kk/Ngrid, L, M, N, Normforce);
							if(veclength(OutForce) < thCritPt) {
								//vector length threshold for critical points

								seeds[numSeeds].x = i + ii/Ngrid;
								seeds[numSeeds].y = j + jj/Ngrid;
								seeds[numSeeds].z = k + kk/Ngrid;
								numSeeds++;
								if(numSeeds >= MAX_NUM_SEEDS) {
									printf("UPS! - Too many seed points found! - Abort.\n");
									exit(1);
								}
								NumCritPoints++;

								// once we find a critical point inside a cell, move to the next cell
								critPtFound = true;
							}
						}
					}
				}
			}
		}
	}

#ifdef _DEBUG
	PrintElapsedTime("\nPhase4.5: Finding critical points 1");
#endif

	printf("Number of critical points is: %d, and number of seeds %d\n", NumCritPoints, numSeeds);
	//fprintf(fout,"%d %d %d %f %f\n", 1, 1, 1, -1.0, -1.0);

#ifdef TRACE
	// print the critical points
	printf("Critical points found in reverse order:\n");
	for(i = numSeeds-1; i >= numSeeds - NumCritPoints; i--) {
		printf("\t%f, %f, %f\n", seeds[i].x, seeds[i].y, seeds[i].z);
	}

#endif

	// reduce the number of critical points.

	exit(1);

	// seeds[0..numBoundSeeds-1] = boundary seeds, seeds[numBoundSeeds.. numSeeds-1] = critical points
	idxSeeds = numSeeds-1;
	while (idxSeeds >= 0) {
		Startpos.x=seeds[idxSeeds].x; Startpos.y=seeds[idxSeeds].y; Startpos.z=seeds[idxSeeds].z;
		idx = (int)Startpos.z * slsz + (int)Startpos.y *L + (int)Startpos.x;

#ifdef TRACE
		printf("Point: (%f, %f, %f): idx = %d, force here: (%f, %f, %f)\n",
			Startpos.x, Startpos.y, Startpos.z, idx,
			ForceField[idx].xd, ForceField[idx].yd, ForceField[idx].zd);
#endif
		//Check whether the high curv points are within D-voxel distance to the existing skeleton
		int Dvoxel = 1;
		int FlagWithin = 0;
		if (idxSeeds < numBoundSeeds) {
			for (kk = -Dvoxel+1; kk <= Dvoxel-1; kk++) {
				for (jj = -Dvoxel+1; jj <= Dvoxel-1; jj++) {
					for (ii = -Dvoxel+1; ii <= Dvoxel-1; ii++) {
						if(FlagOnSkeleton[idx + kk*slsz + jj*L + ii] == 1)  FlagWithin = 1;

					}
				}
			}
			if(FlagWithin == 1) {
				idxSeeds--;
				continue;
			}
		}

		// being able to not show critical point in the streamlines
		(*Skel)[(*numSkel)].x = Startpos.x;
		(*Skel)[(*numSkel)].y = Startpos.y;
		(*Skel)[(*numSkel)].z = Startpos.z;
		(*numSkel) = (*numSkel) + 1;

		if((*numSkel) >= MAX_NUM_SKEL ) {
			printf("UPS! - Too many skeleton points to output ! - Abort.\n");
			exit(1);
		}

		FlagOnSkeleton[idx] = 1;

#ifdef TRACE
		// if(idxSeeds < numBoundSeeds) {
			printf("Starting to follow the force field...\n");
		// }
#endif


		streamSteps = 0;
		while(streamSteps<5000)   {
			rk2(Startpos.x, Startpos.y, Startpos.z, L, M, N, 0.2, ForceField, &Nextpos);

#ifdef TRACE
			//if(idxSeeds < numBoundSeeds) {
				printf("Step %d. Next position: (%f, %f, %f) idx = %d: force here: (%f, %f, %f)\n",
					streamSteps, Nextpos.x, Nextpos.y, Nextpos.z, idx,
					ForceField[idx].xd, ForceField[idx].yd, ForceField[idx].zd);
			//}
#endif


			// if Nextpos = Startpos then, we should stop
			if(	(Nextpos.x == Startpos.x)	&&
				(Nextpos.y == Startpos.y)	&&
				(Nextpos.z == Startpos.z))
			{
#ifdef TRACE
				printf("New position is the same as start position. Break\n");
#endif
				break;
			}

			streamSteps++;


			// idx = (int)Startpos.z *slsz + (int)Startpos.y *L + (int)Startpos.x;
			idx = (int)Nextpos.z *slsz + (int)Nextpos.y *L + (int)Nextpos.x;


			if (FlagOnSkeleton[idx] != 1) {

				(*Skel)[(*numSkel)].x = Nextpos.x;
				(*Skel)[(*numSkel)].y = Nextpos.y;
				(*Skel)[(*numSkel)].z = Nextpos.z;
				(*numSkel) = (*numSkel) + 1;

				if((*numSkel) >= MAX_NUM_SKEL ) {
					printf("UPS! - Too many skeleton points to output ! - Abort.\n");
					exit(1);
				}

				//printf("%f %f %f %d\n", Nextpos.x, Nextpos.y, Nextpos.z, 1);
				FlagOnSkeleton[idx] = 1;
			}
			else {
				// FlagOnSkeleton[idx] is 1. It means we've been here before
				// if we moved to a new cell, then if that cell is already part of the skeleton we can stop.
				if(	((int)Startpos.x != (int)Nextpos.x) ||
					((int)Startpos.y != (int)Nextpos.y) ||
					((int)Startpos.z != (int)Nextpos.z))
				{
#ifdef TRACE
				printf("New position is in already visited cell. Break\n");
#endif

					break;
				}
			}

			Startpos.x = Nextpos.x;
			Startpos.y = Nextpos.y;
			Startpos.z = Nextpos.z;
		}

		idxSeeds--;
	}

#ifdef _DEBUG
	PrintElapsedTime("Phase4.6: Finding critical points 2");
#endif

	// free allocated memory
	// delete [] Normforce;
	delete [] FlagOnSkeleton;

	return true;
}


inline Vector interpolation(float x, float y, float z, int sizx, int sizy, int sizz, Vector *forcevec)
    {
	float alpha, beta, gamma;
	Vector forceInt;
	long slsz;

	alpha=x-int(x);
	beta=y-int(y);
	gamma=z-int(z);
	slsz=sizy*sizx;

	forceInt.xd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*alpha*beta*(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*(alpha*beta*gamma);

	forceInt.yd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*gamma;

	forceInt.zd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*gamma;

	return(forceInt);
    }


inline void rk2(float x, float y, float z, int sizx, int sizy, int sizz, double steps, Vector *Force_ini, VoxelPositionDouble *nextPos)
   {
	long slsz;
	Vector OutForce;
	float x1, y1, z1;
	slsz=sizy*sizx;

	OutForce=interpolation(x,y,z,sizx,sizy,sizz,Force_ini);

	x = x + OutForce.xd * steps;
	y = y + OutForce.yd * steps;
	z = z + OutForce.zd * steps;

	nextPos->x = x;
	nextPos->y = y;
	nextPos->z = z;
   }



