///////////////////////////////////////////////////////////////////////////////////////////////
// ----  makes the 3D dataset a solid volume.
//		fills all the holes in the volume if there are any, and makes the volume binary
//
// Last change: Fri Aug 15 12:49:48 EDT 2003 by Nicu D. Cornea
//
///////////////////////////////////////////////////////////////////////////////////////////////

//
// Transforms a 3D object into a binary object: all the object voxels will have the same value,
//	and the background voxels will have another value.
//
bool MakeSolidVolume(
	int L, int M, int N, 				// [in] volume size
	unsigned char* vol, 				// [in, out] volume 
	unsigned char backgroundVal, 		// [in] desired background value
	unsigned char objectVal				// [in] desired object value
);
