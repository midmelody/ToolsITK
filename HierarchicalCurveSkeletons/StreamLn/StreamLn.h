#include "../common.h"

///////////////////////////////////////////////////////////////////////////////
// Allocates a Skeleton data structure with a maximum points and segments
///////////////////////////////////////////////////////////////////////////////
bool AllocateSkeleton(Skeleton **Skel, int numPoints, int numSegments);

///////////////////////////////////////////////////////////////////////////////
// Frees a Skeleton data structure previously allocated with AllocateSkeleton
///////////////////////////////////////////////////////////////////////////////
bool FreeSkeleton(Skeleton **Skel);

///////////////////////////////////////////////////////////////////////////////
// Get level 1, 2 and 3 skeletons combined in one step
///////////////////////////////////////////////////////////////////////////////
bool GetStreamLines(
	Vector* ForceField, 		// [in] vector field
	int L, int M, int N,		// [in] vector field size (X, Y and Z)
	CriticalPoint *CritPts,         // [in] critical points array
	int numCritPts, 		// [in] number of critical points
	VoxelPosition *BoundarySeeds,   // [in] boundary seeds array
	int numBoundarySeeds,		// [in] number of boundary seeds
	VoxelPositionDouble *HDPoints,  // [in] high divergence points array
	int numHDPoints,		// [in] number of high divergence 
	                                //        points
	Skeleton *Skel                 // [out] skeleton points and segments
);


//////////////////////////////////////////////////////////////////////////////
// Get basic skeleton (level 1): connects only critical points
//////////////////////////////////////////////////////////////////////////////

bool GetLevel1Skeleton(
        Vector* ForceField, 		// [in] vector field
	int L, int M, int N,		// [in] vector field size (X, Y and Z)
	CriticalPoint *CritPts,         // [in] critical points array
	int numCritPts, 		// [in] number of critical points
	Skeleton *Skel                 // [out] skeleton points and segments
);

/////////////////////////////////////////////////////////////////////////////
// GetLevel2Skeleton: connects the high divergence points to the existing
//    skeleton
/////////////////////////////////////////////////////////////////////////////

bool GetLevel2Skeleton(
	Vector* ForceField, 		// [in] vector field
	int L, int M, int N,		// [in] vector field size (X, Y and Z)
	VoxelPositionDouble *HDPoints,  // [in] high divergence points array
	int numHDPoints,		// [in] number of high divergence 
	                                //        points
	Skeleton *Skel                  // [in,out] skeleton points and 
	                                //             segments
);


/////////////////////////////////////////////////////////////////////////////
// GetLevel3Skeleton: connects the boundary seeds points to the existing
//    skeleton
/////////////////////////////////////////////////////////////////////////////
bool GetLevel3Skeleton(
	Vector* ForceField, 		// [in] vector field
	int L, int M, int N,		// [in] vector field size (X, Y and Z)
	VoxelPosition *BoundarySeeds,   // [in] boundary seeds array
	int numBoundarySeeds,		// [in] number of boundary seeds
	Skeleton *Skel                  // [in, out] skeleton points and
	                                //              segments
);
