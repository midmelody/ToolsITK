///////////////////////////////////////////////////////////////////////////////////////////////
// ----  makes the 3D dataset a solid volume.
//		fills all the holes in the volume if there are any, and makes the volume binary
//
// Last change:  by Nicu D. Cornea
//
///////////////////////////////////////////////////////////////////////////////////////////////

#include "makeSolidVol.h"
#include <stdio.h>

#define OBJECT_VAL 1

// 3 values defining outisde voxels
#define OUTSIDE_1 2
#define OUTSIDE_2 3
#define OUTSIDE_3 4


bool FloodFillZPlane(int zPlane, int L, int M, int N, unsigned char* vol);

bool FloodFillYPlane(int yPlane, int L, int M, int N, unsigned char* vol);
bool FloodFillXPlane(int xPlane, int L, int M, int N, unsigned char* vol);

bool MakeSolidVolume(int L, int M, int N, unsigned char* vol, 
	unsigned char backgroundVal, unsigned char objectVal) 
{
	long idx;
	int i;
			
	// replace every object voxel with OBJECT_VAL
	// everything else remains 0
	for(idx=0; idx < L*M*N; idx++) {
		if(vol[idx] != 0) {
			vol[idx] = OBJECT_VAL;
		}
	}
	
	//printf("Floodfill Z planes ...");
	// floodfill the outside of the object with OUTSIDE_1 in the Z direction
	for(i=0; i < N; i++) {
	  printf("\r\tProcessing z plane %d out of %d.", i, N);
	  fflush(stdout);
	  FloodFillZPlane(i, L, M, N, vol);
	}
	printf("done.\n");

	//printf("Floodfill Y planes ...");
	// floodfill the outside of the object with OUTSIDE_2 in the Y direction
	for(i=0; i < M; i++) {
	  printf("\r\tProcessing y plane %d out of %d.", i, M);
	  fflush(stdout);
		FloodFillYPlane(i, L, M, N, vol);
	}
	printf("done.\n");

	//printf("Floodfill X planes ...");
	// floodfill the outside of the object with OUTSIDE_3 in the X direction
	for(i=0; i < L; i++) {
	  printf("\r\tProcessing x plane %d out of %d.", i, L);
	  fflush(stdout);
	  FloodFillXPlane(i, L, M, N, vol);
	}
	printf("done.\n");

	// replace any voxel that's not OUTSIDE_1, _2 or _3  with objectVal
	// and every OUTSIDE_1, _2 or _3  with backgroundVal
	for(idx=0; idx < L*M*N; idx++) {
		if((vol[idx] == OUTSIDE_1) || (vol[idx] == OUTSIDE_2) || (vol[idx] == OUTSIDE_3)) {
			vol[idx] = backgroundVal;
		}
		else {
			vol[idx] = objectVal;
		}
	}
	
	return true;	
}


// floodfill each Z plane with OUTSIDE_1 values
bool FloodFillZPlane(int zPlane, int L, int M, int N, unsigned char* vol) {
  long idx, idxS, idxN, ts;
  bool anyChange = false;
  int x, y;

  ts = L*M*N;
  // set point (0,0) to OUTSIZE_1
  idx = zPlane*L*M /* + 0*L + 0 */;
  vol[idx] = OUTSIDE_1;
  
  anyChange = true;
  while(anyChange) {
    
    anyChange = false;
    // loop from left to right and top to bottom
    for(x=0; x < L; x++) {
      for(y=0; y < M; y++) {
	idxS = idx + y*L + x;
	// if the point is set to OUTSIDE_1, the set all empty neightbors 
	// to OUTSIDE_1
	if(vol[idxS] == OUTSIDE_1) {
	  
	  idxN = idxS + L;
	  if((idxN >= 0) && (idxN < ts) && (vol[idxN] == 0)) { 
	    vol[idxN] = OUTSIDE_1;
	    anyChange = true;
	  }
	  
	  idxN = idxS - L;
	  if((idxN >= 0) && (idxN < ts) && (vol[idxN] == 0)) {
	    vol[idxN] = OUTSIDE_1;
	    anyChange = true;
	  }
	  
	  idxN = idxS + 1;
	  if((idxN >= 0) && (idxN < ts) && (vol[idxN] == 0)) {
	    vol[idxN] = OUTSIDE_1;
	    anyChange = true;
	  }
	  
	  idxN = idxS - 1;
	  if((idxN >= 0) && (idxN < ts) && (vol[idxN] == 0)) {
	    vol[idxN] = OUTSIDE_1;
	    anyChange = true;
	  }
	}
      }
    }
    
    if(anyChange) {
      // same loop but bottom to top and right to left
      anyChange = false;
      // loop from left to right and top to bottom
      for(x=L-1; x >=0; x--) {
	for(y=M-1; y >=0; y--) {
	  idxS = idx + y*L + x;
	  // if the point is set to OUTSIDE_1, the set all empty neightbors 
	  // to OUTSIDE_1
	  if(vol[idxS] == OUTSIDE_1) {
	    
	    idxN = idxS + L;
	    if((idxN >= 0) && (idxN < ts) && (vol[idxN] == 0)) { 
	      vol[idxN] = OUTSIDE_1;
	      anyChange = true;
	    }
	    
	    idxN = idxS - L;
	    if((idxN >= 0) && (idxN < ts) && (vol[idxN] == 0)) {
	      vol[idxN] = OUTSIDE_1;
	      anyChange = true;
	    }
	    
	    idxN = idxS + 1;
	    if((idxN >= 0) && (idxN < ts) && (vol[idxN] == 0)) {
	      vol[idxN] = OUTSIDE_1;
	      anyChange = true;
	    }
	    
	    idxN = idxS - 1;
	    if((idxN >= 0) && (idxN < ts) && (vol[idxN] == 0)) {
	      vol[idxN] = OUTSIDE_1;
	      anyChange = true;
	    }
	  }
	}
      } 
    }
  }
  
  return true;
}


// floodfill each Y plane with OUTSIDE_2 values
bool FloodFillYPlane(int yPlane, int L, int M, int N, unsigned char* vol) {
  long idx, idxS, idxN, ts;
  bool anyChange = false;
  int x, z;

  ts = L*M*N;
  // set point (0,0) to OUTSIZE_2
  idx = /*0*L*M  + */ yPlane*L /*+ 0 */;
  vol[idx] = OUTSIDE_2;
  
  anyChange = true;
  while(anyChange) {
    
    anyChange = false;
    // loop from left to right and top to bottom
    for(x=0; x < L; x++) {
      for(z=0; z < N; z++) {
	      idxS = z*L*M + idx + x;
	      // if the point is set to OUTSIDE_2, the set all empty neightbors 
	      // to OUTSIDE_2
	      if(vol[idxS] == OUTSIDE_2) {
      	  
	  idxN = idxS + L*M;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1))) { 
	    vol[idxN] = OUTSIDE_2;
	    anyChange = true;
	  }
	  
	  idxN = idxS - L*M;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1))) { 
	    vol[idxN] = OUTSIDE_2;
	    anyChange = true;
	  }
	  
	  idxN = idxS + 1;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1))) { 
	    vol[idxN] = OUTSIDE_2;
	    anyChange = true;
	  }
	  
	  idxN = idxS - 1;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1))) { 
	    vol[idxN] = OUTSIDE_2;
	    anyChange = true;
	  }
	}
      }
    }
    
    if(anyChange) {
      // same loop but bottom to top and right to left
      
      anyChange = false;
    // loop from left to right and top to bottom
    for(x=L-1; x >= 0; x--) {
      for(z=N-1; z >= 0; z--) {
	      idxS = z*L*M + idx + x;
	      // if the point is set to OUTSIDE_2, the set all empty neightbors 
	      // to OUTSIDE_2
	      if(vol[idxS] == OUTSIDE_2) {
      	  
	  idxN = idxS + L*M;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1))) { 
	    vol[idxN] = OUTSIDE_2;
	    anyChange = true;
	  }
	  
	  idxN = idxS - L*M;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1))) { 
	    vol[idxN] = OUTSIDE_2;
	    anyChange = true;
	  }
	  
	  idxN = idxS + 1;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1))) { 
	    vol[idxN] = OUTSIDE_2;
	    anyChange = true;
	  }
	  
	  idxN = idxS - 1;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1))) { 
	    vol[idxN] = OUTSIDE_2;
	    anyChange = true;
	  }
	}
      }
    }
    }
  }
  
  return true;
}


// floodfill each X plane with OUTSIDE_3 values
bool FloodFillXPlane(int xPlane, int L, int M, int N, unsigned char* vol) {
  long idx, idxS, idxN, ts;
  bool anyChange = false;
  int y, z;

  ts = L*M*N;
  // set point (0,0) to OUTSIZE_3
  idx = /*0*L*M  +  yPlane*L */+ xPlane ;
  vol[idx] = OUTSIDE_3;
  
  anyChange = true;
  while(anyChange) {
    
    anyChange = false;
    // loop from left to right and top to bottom
    for(y=0; y < M; y++) {
      for(z=0; z < N; z++) {
	      idxS = z*L*M + L*y + idx;
	      // if the point is set to OUTSIDE_3, the set all empty neightbors 
	      // to OUTSIDE_3
	      if(vol[idxS] == OUTSIDE_3) {
      	  
	  idxN = idxS + L*M;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1) || (vol[idxN] == OUTSIDE_2))) { 
	    vol[idxN] = OUTSIDE_3;
	    anyChange = true;
	  }
	  
	  idxN = idxS - L*M;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1) || (vol[idxN] == OUTSIDE_2))) { 
	    vol[idxN] = OUTSIDE_3;
	    anyChange = true;
	  }
	  
	  idxN = idxS + L;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1) || (vol[idxN] == OUTSIDE_2))) { 
	    vol[idxN] = OUTSIDE_3;
	    anyChange = true;
	  }
	  
	  idxN = idxS - L;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1) || (vol[idxN] == OUTSIDE_2))) { 
	    vol[idxN] = OUTSIDE_3;
	    anyChange = true;
	  }
	}
      }
    }
    
    if(anyChange) {
      // same loop but bottom to top and right to left
      
      anyChange = false;
    // loop from left to right and top to bottom
    for(y=M-1; y >= 0; y--) {
      for(z=N-1; z >= 0; z--) {
	      idxS = z*L*M + + L*y + idx;
	      // if the point is set to OUTSIDE_3, the set all empty neightbors 
	      // to OUTSIDE_3
	      if(vol[idxS] == OUTSIDE_3) {
      	  
	  idxN = idxS + L*M;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1)  || (vol[idxN] == OUTSIDE_2))) { 
	    vol[idxN] = OUTSIDE_3;
	    anyChange = true;
	  }
	  
	  idxN = idxS - L*M;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1) || (vol[idxN] == OUTSIDE_2))) { 
	    vol[idxN] = OUTSIDE_3;
	    anyChange = true;
	  }
	  
	  idxN = idxS + L;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1) || (vol[idxN] == OUTSIDE_2))) { 
	    vol[idxN] = OUTSIDE_3;
	    anyChange = true;
	  }
	  
	  idxN = idxS - L;
	  if((idxN >= 0) && (idxN < ts) && 
       ((vol[idxN] == 0) || (vol[idxN] == OUTSIDE_1) || (vol[idxN] == OUTSIDE_2))) { 
	    vol[idxN] = OUTSIDE_3;
	    anyChange = true;
	  }
	}
      }
    }
    }
  }
  
  return true;
}


