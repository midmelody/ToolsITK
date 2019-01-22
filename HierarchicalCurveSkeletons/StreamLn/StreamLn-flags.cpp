// Form streamlines from saddle points and seed points
// --- Input: 1. normalized 3D vector field
//            2. Seed points: 
//			critical points, 
//			high curvature points on the boundary, 
//			and interior high divergence points
// --- Output: skelelton file
// --- Author: Xiaosong Yuan, Bala, Nicu D. Cornea - Vizlab, Rutgers University
// --- Date: 9/28/2002
// --- Original file name: skel_streamline_float.cpp
//
// --- Changed to StreamLn.cpp by Nicu D. Cornea, on Tuesday, July 22, 2003
//
//

// #define TRACE

#include "StreamLn.h"

#define MAX_NUM_SKEL	80000


inline double veclength(Vector vecin) {
    return sqrt(vecin.xd * vecin.xd + vecin.yd * vecin.yd +vecin.zd * vecin.zd);
}
inline Vector interpolation(double x, double y, double z, int sizx, int sizy, int sizz, Vector *forcevec);
inline void rk2(double x, double y, double z, int sizx, int sizy, int sizz, double steps, Vector *Force_ini, VoxelPositionDouble *nextPos);

bool FollowStreamlines(double x, double y, double z,
	int crtFlag, int* FlagOnSkeleton, int* inds,
	int sX, int sY, int sZ, Vector* ForceField,
	VoxelPositionDouble **Skel, int *numSkel);

bool GetStreamLines(
	Vector* ForceField, int L, int M, int N,
	unsigned char *flags,
	CriticalPoint *CritPts, int numCritPts,
	VoxelPosition *BdSeeds, int numBoundSeeds,
	VoxelPositionDouble **Skel, int *numSkel)
{

#ifdef TRACE
	printf("TRACE: Starting GetStreamLines function...\n");
#endif

	// Vector *force;
	// Vector *Normforce;
	long idx, iidx, slsz, prvidx;


	int cc;
	int ii,jj,kk;

	int *FlagOnSkeleton;
	int crtFlagOnSkeleton, maxFlagOnSkeleton;

	float vecLength;

	int i,j,k;

	slsz = L*M;		// slice size

#ifdef TRACE
	printf("Critical points:\n");
	for(i=0; i < numCritPts; i++) {
		printf("%d: (%f %f %f)\n", i, CritPts[i].position.x, CritPts[i].position.y, CritPts[i].position.z);
	}
#endif

	FlagOnSkeleton = new int[L*M*N];
	if((FlagOnSkeleton == NULL)) {
		printf("UPS! - Error allocating memory for working data structures. Abort.\n");
		exit(1);
	}


	(*Skel) = NULL;
	(*numSkel) = 0;

	if(((*Skel) = new VoxelPositionDouble[MAX_NUM_SKEL]) == NULL) {
		printf("UPS! - Error allocating memory for the output array. Abort.\n");
		exit(1);
	}


	// set FlagOnSekeleton to 0
	for (idx=0; idx < L*M*N; idx++) {
		FlagOnSkeleton[idx] = 0;
	}


/////////////////////////////////////
// define neighborhood of a voxel
int inds[30];

inds[0]	=  + slsz + 0 + 0;
inds[1]	=  - slsz + 0 + 0;
inds[2]	=  +    0 + L + 0;
inds[3]	=  +    0 - L + 0;
inds[4]	=  +    0 + 0 + 1;
inds[5]	=  +    0 + 0 - 1;
	// v-neighbors
inds[6]	=  - slsz - L - 1;
inds[7]	=  - slsz - L + 1;
inds[8]	=  - slsz + L - 1;
inds[9]	=  - slsz + L + 1;
inds[10]=  + slsz - L - 1;
inds[11]=  + slsz - L + 1;
inds[12]=  + slsz + L - 1;
inds[13]=  + slsz + L + 1;
	// e-neighbors
inds[14]=  + slsz + L + 0;
inds[15]=  + slsz - L + 0;
inds[16]=  - slsz + L + 0;
inds[17]=  - slsz - L + 0;
inds[18]=  + slsz + 0 + 1;
inds[19]=  + slsz + 0 - 1;
inds[20]=  - slsz + 0 + 1;
inds[21]=  - slsz + 0 - 1;
inds[22]=  +    0 + L + 1;
inds[23]=  +    0 + L - 1;
inds[24]=  +    0 - L + 1;
inds[25]=  +    0 - L - 1;
	// itself
inds[26]=  +    0 + 0 + 0;


	// set the flags for every critical point as being part of the skeleton
	// they have to be set in advance because when we follow the force field,
	//	we will check to see if we are close to one of the critical points by examining
	//	the flags of the neighbors of the current cell
	crtFlagOnSkeleton = 1;
	for(i=0; i < numCritPts; i++) {
		idx = (int)CritPts[i].position.z * slsz + (int)CritPts[i].position.y *L + (int)CritPts[i].position.x;
		FlagOnSkeleton[idx] = crtFlagOnSkeleton;
		crtFlagOnSkeleton++;
	}
	maxFlagOnSkeleton = crtFlagOnSkeleton;

	// follow the streamlines starting at saddles
	// for every critical point that is not an attracting node
	//	follow the streamlines in the direction of the positive eigenvector(s)
	//	until we reach an already visited cell, or one of the critical points



	for(i = 0; i < numCritPts; i++) {

		idx = (int)CritPts[i].position.z * slsz + (int)CritPts[i].position.y *L + (int)CritPts[i].position.x;

#ifdef TRACE
		printf("Point: (%f, %f, %f): idx = %d, force here: (%f, %f, %f)\n",
			CritPts[i].position.x, CritPts[i].position.y, CritPts[i].position.z, idx,
			ForceField[idx].xd, ForceField[idx].yd, ForceField[idx].zd);
#endif



		(*Skel)[(*numSkel)].x = CritPts[i].position.x;
		(*Skel)[(*numSkel)].y = CritPts[i].position.y;
		(*Skel)[(*numSkel)].z = CritPts[i].position.z;
		(*numSkel) = (*numSkel) + 1;

		if((*numSkel) >= MAX_NUM_SKEL ) {
			printf("UPS! - Too many skeleton points to output ! - Abort.\n");
			exit(1);
		}

		crtFlagOnSkeleton = FlagOnSkeleton[idx];


		// start to follow forcefield at saddle points only
		if(CritPts[i].type != CPT_SADDLE) {
			continue;
		}


#ifdef TRACE
		printf("Starting to follow the force field...\n");
#endif

		// the point is a saddle, so we should go in the direction pointed by the pozitive eigenvectors
		for(j=0; j < 3; j++) {
			if(CritPts[i].eval[j] > 0) {
				FollowStreamlines(
							CritPts[i].position.x + (CritPts[i].evect[j].xd * 0.80),
							CritPts[i].position.y + (CritPts[i].evect[j].yd * 0.80),
							CritPts[i].position.z + (CritPts[i].evect[j].zd * 0.80),
							crtFlagOnSkeleton, FlagOnSkeleton, inds,
							L, M, N, ForceField, Skel, numSkel);

				FollowStreamlines(
							CritPts[i].position.x - (CritPts[i].evect[j].xd * 0.80),
							CritPts[i].position.y - (CritPts[i].evect[j].yd * 0.80),
							CritPts[i].position.z - (CritPts[i].evect[j].zd * 0.80),
							crtFlagOnSkeleton, FlagOnSkeleton, inds,
							L, M, N, ForceField, Skel, numSkel);

			}
			//
		}



	}

#ifdef _DEBUG
	PrintElapsedTime("\tSL-2: Following streamlines starting at the critical points.");
#endif

	// free allocated memory
	// delete [] Normforce;
	delete [] FlagOnSkeleton;



	return true;
}







inline Vector interpolation(double x, double y, double z, int sizx, int sizy, int sizz, Vector *forcevec)
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


inline void rk2(double x, double y, double z, int sizx, int sizy, int sizz, double steps, Vector *Force_ini, VoxelPositionDouble *nextPos)
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

#ifdef TRACE
	//if(idxSeeds < numBoundSeeds) {

	printf("Next position: (%f, %f, %f) force here: (%f, %f, %f)\n",
		nextPos->x, nextPos->y, nextPos->z,
		OutForce.xd, OutForce.yd, OutForce.zd);

	//}
#endif
   }







bool FollowStreamlines(double x, double y, double z,
	int crtFlag, int* FlagOnSkeleton, int* inds,
	int sX, int sY, int sZ, Vector* ForceField,
	VoxelPositionDouble **Skel, int *numSkel) {

long idx, prividx, slsz;
int i, streamSteps;
VoxelPositionDouble Startpos, Nextpos;


	slsz = sX*sY;

	Startpos.x = x;
	Startpos.y = y;
	Startpos.z = z;

	idx = (int)Startpos.z *slsz + (int)Startpos.y *sX + (int)Startpos.x;
	// if this point is already marked as part of the skeleton but with a different flag than the current flag,
	// then we should return
	if(	(FlagOnSkeleton[idx] != 0) &&
		(FlagOnSkeleton[idx] != crtFlag))
	{
		return true;
	}

	// if not,
	// add the point to the skeleton
	(*Skel)[(*numSkel)].x = Startpos.x;
	(*Skel)[(*numSkel)].y = Startpos.y;
	(*Skel)[(*numSkel)].z = Startpos.z;
	(*numSkel) = (*numSkel) + 1;

	if((*numSkel) >= MAX_NUM_SKEL ) {
		printf("UPS! - Too many skeleton points to output ! - Abort.\n");
		exit(1);
	}
	FlagOnSkeleton[idx] = crtFlag;


	streamSteps = 0;
	while(streamSteps < 5000)   {

#ifdef TRACE
		//if(idxSeeds < numBoundSeeds) {
			printf("Step %d.", streamSteps);
		//}
#endif

		rk2(Startpos.x, Startpos.y, Startpos.z, sX, sY, sZ, 0.2, ForceField, &Nextpos);

		// if Nextpos = Startpos then, we should stop
		if(	EQUAL(Nextpos.x, Startpos.x)	&&
			EQUAL(Nextpos.y, Startpos.y)	&&
			EQUAL(Nextpos.z, Startpos.z))
		{

#ifdef TRACE
			printf("New position is the same as start position. Break\n");
#endif
			return true;;
		}

		streamSteps++;

		prividx = idx;
		idx = (int)Nextpos.z *slsz + (int)Nextpos.y *sX + (int)Nextpos.x;
#ifdef TRACE
		printf("Index of new position is: %d\n", idx);
#endif


		// if we just moved to a new cell,
		if(	((int)Startpos.x != (int)Nextpos.x) ||
			((int)Startpos.y != (int)Nextpos.y) ||
			((int)Startpos.z != (int)Nextpos.z))
		{
#ifdef TRACE
			printf("New cell.\n");
#endif

			if (FlagOnSkeleton[idx] == 0) {
				// if that cell is not already part of the skeleton
				// add it to the skeleton

				(*Skel)[(*numSkel)].x = Nextpos.x;
				(*Skel)[(*numSkel)].y = Nextpos.y;
				(*Skel)[(*numSkel)].z = Nextpos.z;
				(*numSkel) = (*numSkel) + 1;

				if((*numSkel) >= MAX_NUM_SKEL ) {
					printf("UPS! - Too many skeleton points to output ! - Abort.\n");
					exit(1);
				}

				FlagOnSkeleton[idx] = crtFlag;

				// we should keep going. so we will reset the streamSteps counter
				streamSteps = 0;

				// if this new cell has a neighbor that is already part of the skeleton, stop
				// but not if it's a neighbor from the current line (same flag value )
				for(i=0; i < 26; i++) {
					if(	((idx + inds[i]) >= 0) &&
						((idx + inds[i]) < sX*sY*sZ))

					{

						if(	(FlagOnSkeleton[idx + inds[i]] != 0) &&
							(FlagOnSkeleton[idx + inds[i]] != crtFlag))
						{
							// stop here
#ifdef TRACE
							printf("It has a neighbor that is already in the skeleton. Skip to the next seed.\n");
#endif
							return true;
						}
					}
				}

			}
			else {
				// FlagOnSkeleton[idx] is not 0. It means we've been here before, so we can stop.
#ifdef TRACE
					printf("New position is in already visited cell. Break\n");
#endif
					return true;

			}


		}

		Startpos.x = Nextpos.x;
		Startpos.y = Nextpos.y;
		Startpos.z = Nextpos.z;
	}

	return true;
}

