#include "../common.h"

bool GetHighCurvatureBoundaryPoints(
	unsigned char *vol, int L, int M, int N,
	VoxelPosition ** HCBPoints, int *numHCBPoints,
	double perc);
