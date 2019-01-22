// Extract the boundary points from a input dataset: volume or polygonal mesh
// --- Input: 3D data set or polygonal mesh
// --- Output: boundary points with neighbourhood information and normals structure
//
// Nicu D. Cornea - Mon Jun  2 12:26:56 EDT 2003
#define TRACE


#include "BoundPoints.h"

#define MAX_BOUND_SIZE		20000
#define MAX_NEIGH_SIZE		40

#define	ALL_24_NEIGHBORS
// #define ONLY_FACE_NEIGHBORS

//#define EXTENDED_NEIGHBORHOOD		// for this you need to set MAX_NEIGH_SIZE to a larger number than 27.
//#define	EN_MINDIST			2
//#define EN_MAXDIST			3

bool GetBoundaryPoints(											\
	unsigned char* volume, int sX, int sY, int sZ, 			\
	bdElement **Boundary, int *numBound)
{
	unsigned int Xm1, Ym1, Zm1;
	unsigned int i, j, k;
	bool flagSurf;
	unsigned int idx, iidx;
	unsigned int slsz;
	unsigned char maxn;
	unsigned int* tmp;
	unsigned char pxm1, pxp1, pym1, pyp1, pzm1, pzp1;
	double length;

#ifdef _DEBUG
	PrintElapsedTime("GetBoundaryPoints (GBP)- start.");
#endif

	if(((*Boundary) = new bdElement[MAX_BOUND_SIZE]) == NULL) {
		printf("ERROR: allocating memory for boundary structure. Abort.\n");
		exit(1);
	}
	*numBound = 0;

#ifdef _DEBUG
	PrintElapsedTime("GBP-1: allocating memory.");
#endif


	Xm1 = sX - 1;
  	Ym1 = sY - 1;
  	Zm1 = sZ - 1;
	slsz = sX*sY;

	for (k = 1; k < Zm1; k++) {
		for (j = 1; j < Ym1; j++) {
			for (i = 1; i < Xm1; i++) {
				flagSurf = false;
				idx = k*slsz + j*sX + i;
				pxm1 = 1;	pym1 = 1;	pzm1 = 1;
				pxp1 = 1;	pyp1 = 1;	pzp1 = 1;

				if (volume[idx] == 0) continue;

				//consider six face neighbors, if anyone is zero, it is a boundary voxel

				if (volume[idx - 1] == 0) {
					flagSurf = true;
					pxm1 = 0;
				}

				if (volume[idx + 1] == 0) {
					flagSurf = true;
					pxp1 = 0;
				}

				if (volume[idx -sX] == 0) {
					flagSurf = true;
					pym1 = 0;
				}

				if (volume[idx + sX] == 0) {
					flagSurf = true;
					pyp1 = 0;
				}

				if (volume[idx -slsz] == 0) {
					flagSurf = true;
					pzm1 = 0;
				}

				if (volume[idx + slsz] == 0) {
					flagSurf = true;
					pzp1 = 0;
				}

#ifdef TRACE
				if((i == 76)	&&	(j == 17)	&& (k == 8)) {
					printf("Reached (76, 17, 8): flagsurf = %d, \n\tpxm1 = %d, pxp1 = %d, pym1 = %d, pyp1 = %d, pzm1 = %d, pzp1 = %d\n",
						flagSurf, pxm1, pxp1, pym1, pyp1, pzm1, pzp1);
				}
#endif

				if (flagSurf) {
					(*Boundary)[*numBound].point.x = i;
					(*Boundary)[*numBound].point.y = j;
					(*Boundary)[*numBound].point.z = k;
					(*Boundary)[*numBound].normal.xd = pxm1 - pxp1;
					(*Boundary)[*numBound].normal.yd = pym1 - pyp1;
					(*Boundary)[*numBound].normal.zd = pzm1 - pzp1;
					// normalize the normal vector
					length = 	sqrt(((*Boundary)[*numBound].normal.xd * (*Boundary)[*numBound].normal.xd) +
								((*Boundary)[*numBound].normal.yd * (*Boundary)[*numBound].normal.yd) +
								((*Boundary)[*numBound].normal.zd * (*Boundary)[*numBound].normal.zd));
					if(length != 0) {
						(*Boundary)[*numBound].normal.xd = (*Boundary)[*numBound].normal.xd / length;
						(*Boundary)[*numBound].normal.yd = (*Boundary)[*numBound].normal.yd / length;
						(*Boundary)[*numBound].normal.zd = (*Boundary)[*numBound].normal.zd / length;
					}
#ifdef _DEBUG
					else {

						printf("WARNING: Normal of point (%d, %d, %d) is (0, 0, 0) !!\n");
					}
#endif
					(*Boundary)[*numBound].curvature = 0;
					(*Boundary)[*numBound].nrOfNeighbors = 0;
					if(((*Boundary)[*numBound].neighbors = new unsigned int[MAX_NEIGH_SIZE]) == NULL) {
						printf("ERROR: allcoating memory for the boundary points !! - Abort.\n");
						exit(1);
					}


					(*numBound) = (*numBound) + 1;
					if(*numBound >= MAX_BOUND_SIZE) {
						printf("ERROR: too many boundary points detected !! - Abort.\n");
						exit(1);
					}
				}
			}
		}
	}

#ifdef _DEBUG
	PrintElapsedTime("GBP-2: finding boundary points.");
	printf("-- Found %d points on the boundary\n", *numBound);
#endif

	// completing the neighborhood information
	maxn = 0;
	for(i=0; i < (*numBound); i++) {
		// printf("** Voxel (%d, %d, %d)\n", (*Boundary)[i].point.x, (*Boundary)[i].point.y, (*Boundary)[i].point.z);
		for(j=i+1; j < (*numBound); j++) {

#ifdef ONLY_FACE_NEIGHBORS
			if(	((abs((*Boundary)[i].point.x	- 	(*Boundary)[j].point.x) == 1) 	&&
				 ((*Boundary)[i].point.y 	== 	(*Boundary)[j].point.y) 			&&
				 ((*Boundary)[i].point.z 	== 	(*Boundary)[j].point.z))
				||
			   	((abs((*Boundary)[i].point.y	- 	(*Boundary)[j].point.y) == 1) 	&&
				 ((*Boundary)[i].point.x 	== 	(*Boundary)[j].point.x) 			&&
				 ((*Boundary)[i].point.z 	== 	(*Boundary)[j].point.z))
				||
			   	((abs((*Boundary)[i].point.z	- 	(*Boundary)[j].point.z) == 1) 	&&
				 ((*Boundary)[i].point.x 	== 	(*Boundary)[j].point.x) 			&&
				 ((*Boundary)[i].point.y 	== 	(*Boundary)[j].point.y))
			  )
			{
#else
	#ifdef	EXTENDED_NEIGHBORHOOD
			length = sqrt(((*Boundary)[i].point.x - (*Boundary)[j].point.x)*((*Boundary)[i].point.x - (*Boundary)[j].point.x) +
						((*Boundary)[i].point.y - (*Boundary)[j].point.y)*((*Boundary)[i].point.y - (*Boundary)[j].point.y) +
						((*Boundary)[i].point.z - (*Boundary)[j].point.z)*((*Boundary)[i].point.z - (*Boundary)[j].point.z));
			if((length >= EN_MINDIST) && (length <= EN_MAXDIST))
			{
	#else
			// delfault - use  ALL_24_NEIGHBORS
			if(	(abs((*Boundary)[i].point.x - (*Boundary)[j].point.x) <= 1) &&
				(abs((*Boundary)[i].point.y - (*Boundary)[j].point.y) <= 1) &&
				(abs((*Boundary)[i].point.z - (*Boundary)[j].point.z) <= 1))
			{
	#endif
#endif
				// printf("  - neighbor at (%d, %d, %d)\n", (*Boundary)[j].point.x, (*Boundary)[j].point.y, (*Boundary)[j].point.z);
				(*Boundary)[i].neighbors[(*Boundary)[i].nrOfNeighbors] = j;
				(*Boundary)[i].nrOfNeighbors = (*Boundary)[i].nrOfNeighbors + 1;
				if((*Boundary)[i].nrOfNeighbors >= MAX_NEIGH_SIZE) {
					printf("ERROR: neighborhood size is too small at voxel %d (%d, %d, %d) - %d neighbors.!! - Abort.\n", i, (*Boundary)[i].point.x, (*Boundary)[i].point.y, (*Boundary)[i].point.z, (*Boundary)[i].nrOfNeighbors);
					exit(1);
				}

				(*Boundary)[j].neighbors[(*Boundary)[j].nrOfNeighbors] = i;
				(*Boundary)[j].nrOfNeighbors = (*Boundary)[j].nrOfNeighbors + 1;
				if((*Boundary)[j].nrOfNeighbors >= MAX_NEIGH_SIZE) {
					printf("ERROR: neighborhood size is too small at voxel %d (%d, %d, %d) - %d neighbors !! - Abort.\n", j, (*Boundary)[j].point.x, (*Boundary)[j].point.y, (*Boundary)[j].point.z, (*Boundary)[j].nrOfNeighbors);
					exit(1);
				}
			}
		}

		if((*Boundary)[i].nrOfNeighbors > maxn) {
			maxn = (*Boundary)[i].nrOfNeighbors;
		}

		// delete unnecessary allocated memory
		tmp = new unsigned int[(*Boundary)[i].nrOfNeighbors];
		for(k=0; k < (*Boundary)[i].nrOfNeighbors; k++) {
			tmp[k] = (*Boundary)[i].neighbors[k];
		}
		delete [] (*Boundary)[i].neighbors;
		(*Boundary)[i].neighbors = tmp;
		tmp = NULL;
	}

#ifdef TRACE
	printf("Looking for (76, 17, 8) ...");
	for(i=0; i < (*numBound); i++) {
		if(	((*Boundary)[i].point.x == 76)	&&
			((*Boundary)[i].point.y == 17)	&&
			((*Boundary)[i].point.z == 8))
		{
			printf("FOUND\n");
			printf("Nr of neighbors: %d\n", (*Boundary)[i].nrOfNeighbors);
			break;
		}
	}
#endif

#ifdef _DEBUG
	PrintElapsedTime("GBP-3: filling in neighborhood information.");
	printf("-- maximum nr of neighbors: %u\n", maxn);
#endif

	return true;
}
