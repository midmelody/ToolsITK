#include "def.h"


bool GetBoundaryPoints(										\
	unsigned char* volume, int sX, int sY, int sZ, 		\
	bdElement **Boundary, int *numBound);

// The volume should have 0 for the background voxels and != 0 for objet voxels.
