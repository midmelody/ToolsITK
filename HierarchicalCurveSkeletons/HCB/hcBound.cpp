#include "hcBound.h"

#include "BoundPoints.h"
#include "hcBoundHaimeiri3.h"

bool GetHighCurvatureBoundaryPoints(
	unsigned char *vol, int L, int M, int N,
	VoxelPosition ** HCBPoints, int *numHCBPoints,
	double perc)
{
	bdElement* Bound;
	int numBound, i;


	Bound = NULL;
	GetBoundaryPoints(vol, L, M, N,	&Bound, &numBound);


#ifdef _DEBUG
	PrintElapsedTime("Finding boundary voxels.");
#endif

	(*HCBPoints) = NULL;
	(*numHCBPoints) = 0;
	GetHighCurvaturePoints(Bound, numBound, HCBPoints, numHCBPoints, perc);

	// free the memory occupied by the Bound array
	for(i=0; i < numBound; i++) {
		if(Bound[i].neighbors != NULL) {
			delete [] Bound[i].neighbors;
			Bound[i].nrOfNeighbors = 0;
		}
	}
	delete [] Bound;
	numBound = 0;

#ifdef _DEBUG
	PrintElapsedTime("Detecting high curvature points.");
#endif

	return true;
}
