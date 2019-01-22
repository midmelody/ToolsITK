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

#define SP_SMALL_DISTANCE	0.5		// defines "close to ..."
	// in terms of Manhattan distance |x1-x2| + |y1-y2| + |z1-z2|
	//

#define HD_SMALL_DISTANCE	2.00	// defines how close a high divergence
        // point should be
	// to the skeleton, for it to be ignored

inline Vector interpolation(double x, double y, double z, int sizx, int sizy,
int sizz, Vector *forcevec);
inline void rk2(double x, double y, double z, int sizx, int sizy, int sizz,
double steps,
		Vector *Force_ini, VoxelPositionDouble *nextPos);

inline double Distance(	double x1, double y1, double z1,
						double x2, double y2, double z2);


bool CloseToSkel(	double x, double y, double z,
			VoxelPositionDouble **Skel, int nrAlreadyInSkel,
			double maxDistance);

bool CloseToCP (	double x, double y, double z,
			CriticalPoint *CritPts, int numCritPts, int	originCP,
			double maxDistance);


bool FollowStreamlines(int originCP, double x, double y, double z,
	int sX, int sY, int sZ, Vector* ForceField,
	CriticalPoint *CritPts, int numCritPts,
	VoxelPositionDouble **Skel, int *numSkel);




bool GetStreamLines(
	Vector* ForceField, int L, int M, int N,
	unsigned char *flags,
	CriticalPoint *CritPts, int numCritPts,
	VoxelPosition *BdSeeds, int numBoundSeeds,
	VoxelPositionDouble *HDPoints, int numHDPoints,
	VoxelPositionDouble **Skel, int *numSkel)
{

#ifdef TRACE
	printf("TRACE: Starting GetStreamLines function...\n");
#endif

	long idx, iidx, slsz, prvidx;


	int cc;
	int ii,jj,kk;

	int i,j,k;

	slsz = L*M;		// slice size

#ifdef TRACE
	printf("Critical points:\n");
	for(i=0; i < numCritPts; i++) {
		printf("%d: (%f %f %f)\n", i, CritPts[i].position.x,
		CritPts[i].position.y, CritPts[i].position.z);
	}
#endif


	(*Skel) = NULL;
	(*numSkel) = 0;

	if(((*Skel) = new VoxelPositionDouble[MAX_NUM_SKEL]) == NULL) {
		printf("UPS! - Error allocating memory for the output array. Abort.\n");
		exit(1);
	}


	// follow the streamlines starting at saddles in the direction of the
	// positive eigenvector(s)
	//	until we are close to one of the skeleton points or a critical point,
	// ignoring the points in the
	//	current streamline.

	for(i = 0; i < numCritPts; i++) {

		idx = (int)CritPts[i].position.z * slsz + (int)CritPts[i].position.y *L
		+ (int)CritPts[i].position.x;

#ifdef TRACE
		printf("Point: (%f, %f, %f): idx = %d, force here: (%f, %f, %f)\n",
			CritPts[i].position.x, CritPts[i].position.y, CritPts[i].position.z,
			idx,
			ForceField[idx].xd, ForceField[idx].yd, ForceField[idx].zd);
#endif

		// start to follow forcefield at saddle points only
		if(CritPts[i].type != CPT_SADDLE) {
			continue;
		}


#ifdef TRACE
		printf("Starting to follow the force field...\n");
#endif

		// the point is a saddle, so we should go in the direction pointed by
		// the pozitive eigenvectors
		for(j=0; j < 3; j++) {
			if(CritPts[i].eval[j] > 0) {
				FollowStreamlines(
							i,
							CritPts[i].position.x + (CritPts[i].evect[j].xd * SP_SMALL_DISTANCE),
							CritPts[i].position.y + (CritPts[i].evect[j].yd * SP_SMALL_DISTANCE),
							CritPts[i].position.z + (CritPts[i].evect[j].zd * SP_SMALL_DISTANCE),
							L, M, N, ForceField,
							CritPts, numCritPts,
							Skel, numSkel);

				FollowStreamlines(
							i,
							CritPts[i].position.x - (CritPts[i].evect[j].xd * SP_SMALL_DISTANCE),
							CritPts[i].position.y - (CritPts[i].evect[j].yd * SP_SMALL_DISTANCE),
							CritPts[i].position.z - (CritPts[i].evect[j].zd * SP_SMALL_DISTANCE),
							L, M, N, ForceField,
							CritPts, numCritPts,
							Skel, numSkel);

			}
		}



	}

	// add all the critical points to the skeleton
	for(i=0; i < numCritPts; i++) {
		(*Skel)[(*numSkel)].x = CritPts[i].position.x;
		(*Skel)[(*numSkel)].y = CritPts[i].position.y;
		(*Skel)[(*numSkel)].z = CritPts[i].position.z;
		(*numSkel) = (*numSkel) + 1;

		if((*numSkel) >= MAX_NUM_SKEL ) {
			printf("UPS! - Too many skeleton points to output (%d)! - Abort.\n",
			(*numSkel));
			exit(1);
		}
	}

#ifdef _DEBUG
	PrintElapsedTime("\tSL-2: Following streamlines starting at the critical points.");
#endif

	// start at high divergence points and follow streamlines
	if(HDPoints != NULL) {
		for(i=0; i < numHDPoints; i++) {
			// if this point is close to the existing skeleton, just ignore it
			if(CloseToSkel(HDPoints[i].x, HDPoints[i].y, HDPoints[i].z,
							Skel, *numSkel, HD_SMALL_DISTANCE))
			{
				// ignore this point and move on
				continue;
			}
			
			// it's not close enough
			// follow streamline starting here
			FollowStreamlines(
							-1,
							HDPoints[i].x, HDPoints[i].y, HDPoints[i].z,
							L, M, N, ForceField,
							CritPts, numCritPts,
							Skel, numSkel);
		}
	}

	return true;
}







inline Vector interpolation(double x, double y, double z, int sizx, int sizy,
int sizz, Vector *forcevec)
    {
	float alpha, beta, gamma;
	Vector forceInt;
	long slsz;

	alpha=x-int(x);
	beta=y-int(y);
	gamma=z-int(z);
	slsz=sizy*sizx;

	forceInt.xd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-
	beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*
			(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*alpha*beta*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*
			beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*(alpha*
			beta*gamma);

	forceInt.yd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-
	beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*
			(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*
			beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*
			beta*gamma;

	forceInt.zd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-
	beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*
			(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*
			beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*
			beta*gamma;

	return(forceInt);
    }


inline void rk2(double x, double y, double z, int sizx, int sizy, int sizz,
double steps, Vector *Force_ini, VoxelPositionDouble *nextPos)
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







bool FollowStreamlines(int originCP, double x, double y, double z,
	int sX, int sY, int sZ, Vector* ForceField,
	CriticalPoint *CritPts, int numCritPts,
	VoxelPositionDouble **Skel, int *numSkel)
{

	long idx, prividx, slsz;
	int i, streamSteps;
	VoxelPositionDouble Startpos, Nextpos, LastAddedSkelPoint;
	int nrAlreadyInSkel;
	int step;
	bool stop;

	nrAlreadyInSkel = (*numSkel);

	Startpos.x = x;
	Startpos.y = y;
	Startpos.z = z;

	idx = (int)Startpos.z *slsz + (int)Startpos.y *sX + (int)Startpos.x;

	// add the point to the skeleton
	(*Skel)[(*numSkel)].x = Startpos.x;
	(*Skel)[(*numSkel)].y = Startpos.y;
	(*Skel)[(*numSkel)].z = Startpos.z;
	(*numSkel) = (*numSkel) + 1;

	if((*numSkel) >= MAX_NUM_SKEL ) {
		printf("UPS! - Too many skeleton points to output ! - Abort.\n");
		exit(1);
	}

	LastAddedSkelPoint.x = Startpos.x;
	LastAddedSkelPoint.y = Startpos.y;
	LastAddedSkelPoint.z = Startpos.z;

	step = 0;
	stop = false;
	while(!stop)   {
		step++;
		if(step > 8000) {
			// stop if more than 8000 steps were performed
#ifdef TRACE
		printf("Step counter is %d. Break", step);
#endif
			return true;
		}
#ifdef TRACE
		printf("Step %d.", step);
#endif

		// if this point is close to another point that is already part of the
		// skeleton, or to a critical point
		// we should stop
		if(	CloseToSkel(	Startpos.x, Startpos.y, Startpos.z,
							Skel, nrAlreadyInSkel, SP_SMALL_DISTANCE)
			||
			CloseToCP(	Startpos.x, Startpos.y, Startpos.z,
						CritPts, numCritPts, originCP, SP_SMALL_DISTANCE)) 
		{
#ifdef TRACE
			printf("Position is close to a skeleton point. Break\n");
#endif
			return true;
		}

		rk2(Startpos.x, Startpos.y, Startpos.z, sX, sY, sZ, 0.2,
			ForceField, &Nextpos);

		// if Nextpos = Startpos then, we should stop
		if(	EQUAL(Nextpos.x, Startpos.x)	&&
			EQUAL(Nextpos.y, Startpos.y)	&&
			EQUAL(Nextpos.z, Startpos.z))
		{

#ifdef TRACE
			printf("New position is the same as start position. Break\n");
#endif
			return true;
		}

		// if distance from last added position to the next is larger than ...
		//	add the next position to the skeleton.
		// if not, don't add the next position to the skeleton
		if(Distance(LastAddedSkelPoint.x, LastAddedSkelPoint.y,	LastAddedSkelPoint.z,
					Nextpos.x, Nextpos.y, Nextpos.z) > SP_SMALL_DISTANCE)
		{
			(*Skel)[(*numSkel)].x = Nextpos.x;
			(*Skel)[(*numSkel)].y = Nextpos.y;
			(*Skel)[(*numSkel)].z = Nextpos.z;
			(*numSkel) = (*numSkel) + 1;

			if((*numSkel) >= MAX_NUM_SKEL ) {
				printf("UPS! - Too many skeleton points to output ! - Abort.\n");
				exit(1);
			}

			LastAddedSkelPoint.x = Nextpos.x;
			LastAddedSkelPoint.y = Nextpos.y;
			LastAddedSkelPoint.z = Nextpos.z;
			
		}

		// next position becomes current position
		Startpos.x = Nextpos.x;
		Startpos.y = Nextpos.y;
		Startpos.z = Nextpos.z;
	}

	return true;
}


inline double Distance(	double x1, double y1, double z1,
						double x2, double y2, double z2) {
	return (fabs(x1-x2) + fabs(y1-y2) + fabs(z1-z2));
}


bool CloseToSkel(	double x, double y, double z,
			VoxelPositionDouble **Skel, int nrAlreadyInSkel, 
			double maxDistance)
{
	int i;
	double a, b, c;

	// see if it's close to a skeleton point
	for(i=0; i < nrAlreadyInSkel; i++) {
		// incremental testing of the "close to" condition for higher speed
		if(fabs((*Skel)[i].x - x) < maxDistance) {
			a = fabs((*Skel)[i].x - x);
			if((a + fabs((*Skel)[i].y - y)) < maxDistance) {
				b = fabs((*Skel)[i].y - y);
				if((a + b + fabs((*Skel)[i].z - z)) < maxDistance) {
					// the point is "close" to a skeleton point
					return true;
				}
			}
		}
	}

	return false;
}


bool CloseToCP (	double x, double y, double z,
			CriticalPoint *CritPts, int numCritPts, int	originCP, 
			double maxDistance)
{
	int i;
	double a, b, c;
	
	// see if it's close to a critical point except the one specified in originCP
	for(i=0; i < numCritPts; i++) {
		// ignore the critical point that started this streamline
		if(i != originCP) {
			if(fabs(CritPts[i].position.x - x) < maxDistance) {
				a = fabs(CritPts[i].position.x - x);
				if((a + fabs(CritPts[i].position.y - y)) < maxDistance) {
					b = fabs(CritPts[i].position.y - y);
					if((a + b + fabs(CritPts[i].position.z - z)) < maxDistance) {
						// the point is "close" to a critical point
						return true;
					}
				}
			}
		}
	}
	
	return false;
}
