//
// Expands the volume by adding one more layer of voxels arround the object.
// !!! ATENTION !!! This function works on the flags array, not the original volume
//
#include "../common.h"


// !! ATENTION
// The flags passed here as parameters should be set to INTERIOR for any object voxel
//	and to EXTERIOR (0) for every background voxel
//	NO SURF or BOUNDARY flags !!!

bool ExpandNLayers(
	int L, int M, int N, 			// [in] size of volume
	unsigned char* f,				// [in, out] volume
	int nrLayers);					// [in] nr of layers to add

