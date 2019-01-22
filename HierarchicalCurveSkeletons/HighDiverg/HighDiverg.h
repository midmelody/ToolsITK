#include "../common.h"

bool GetHighDivergencePoints(
	Vector* ForceField, 	      // [in] vector field
	int L, int M, int N,	      // [in] size of vector field (X, Y and Z)
	unsigned char *flags,	      // [in] flags array
	float perc,		      // [in] percentage of high div. points 
	                              //         to be returned (top <perc> %)
	VoxelPositionDouble **HDPts,  // [out] high divergence point list
	int *numHDPts,		      // [out] number of points in the list
	bool inOut = false            // [in] flag specifying if we should look
	                              //    outside the object too (if true).
	                              // DEFAULT: false
);

