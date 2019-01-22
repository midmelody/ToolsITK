// ----  Adds one layer arround the object
// ----  Input: flags, dimensions: X, Y, Z
// ----	 Output: new flags:
//
//

#include "expandVol.h"
////////////////////
// PeelVolume - removes the first layer of voxels from the volume
//   !! It may disconnect the volume !!
///////////////////
bool PeelVolume(
        unsigned char* f,
	int L, int M, int N
) {
  // mark all surface voxels as SURF
  // and then remove them
  
  long idx;
  
  //
  // mark all surface voxels as SURF
  //
  FlagVolume(f, L, M, N);
  
  //
  // remove all surface voxels
  //
  for (idx = 0; idx < L*M*N; idx++) {
    if(f[idx] == SURF) {
      f[idx] = EXTERIOR;
    }
  }

#ifdef _DEBUG
  PrintElapsedTime("\tPV-1: peeling one layer of object voxels.");
#endif

  return true;
}

bool ExpandVolume(
        int L, int M, int N, 
	unsigned char* f, 
	unsigned char padTarget, 
	unsigned char padCharacter
) {
  
  int i, j, k, ii;
  long idx, iidx;
  long slsz, sz;
  
  slsz = L*M;		// slice size
  sz = L*M*N;
  
  // neighbors:
  int ng[6];
  // face neighbors
  ng[0]	= + slsz + 0 + 0;
  ng[1]	= - slsz + 0 + 0;
  ng[2]	= +    0 + L + 0;
  ng[3]	= +    0 - L + 0;
  ng[4]	= +    0 + 0 + 1;
  ng[5]	= +    0 + 0 - 1;
	
  // look at face neighbors only because that gives the best result, 
  //    compared with looking at
  //	face and edge or face, edge and vertex.
	
  for (k = 0; k < N; k++) {
    for (j = 0; j < M; j++) {
      for (i = 0; i < L; i++) {
	idx = k*slsz + j*L + i;
	
	if(f[idx] == padTarget) {
	  // look at all the face neighbors of this voxel
	  //	and replace all the "blank" neighbors with the new value
	  for(ii = 0; ii < 6; ii++) {
	    iidx = idx + ng[ii];
	    if((iidx < 0) || (iidx >= sz)) {
	      printf("** Error padding the object. Volume bounds are too tight. Abort !\n");
	      exit(1);
	    }
	    if(f[iidx] == EXTERIOR) {
	      f[iidx] = padCharacter;
	    }
	  }
	}
      }
    }
  }

#ifdef _DEBUG
  PrintElapsedTime("\tEV-1: padding the object with one layer of extra voxels.");
#endif

  return true;
}


bool ExpandNLayers(
	int L, int M, int N, 			// [in] size of volume
	unsigned char* f,				// [in, out] volume
	int nrLayers)
{
  int i;
  
  if(nrLayers > NR_PAD_VALUES) {
    printf("Error padding object. Maximum %d layers allowed. Abort\n", 
	   NR_PAD_VALUES);
    exit(1);
  }
	
  if(nrLayers == 0) {
    return true;
  }

  if(nrLayers > 0) {
    //
    // exapand
    //
    // first layer attaches itself to the SURF voxels
    ExpandVolume(L, M, N, f, SURF, PADDING_MIN);
	
    // all the other layers attach themselves to the previous one.
    for(i=1; i < nrLayers; i++) {
      ExpandVolume(L, M, N, f, PADDING_MIN + (i - 1), PADDING_MIN + i);
    }
  }
  else {
    // nrLayers < 0
    //
    // peel
    //
    for(i=0; i > nrLayers; i--) {
      PeelVolume(f, L, M, N);
    }
  }
  
  return true;
}

