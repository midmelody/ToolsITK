// ----
// ----  Computes the potential field for a volume
// ----  Input: volume file, dimensions: X, Y, Z, output file name
// ----	 Output: normalized potential field:
//		 1 vector for each point in the volume
//
// Last change: Thu May 15 15:20:38 EDT 2003 by Nicu D. Cornea
//
//

// #define TRACE


#include "potVect.h"

#define BOUND_SIZE	1200000



bool GetIndexOfBPInXYZRange(
	short sx, short sy, short sz,
	short ex, short ey, short ez,
	VoxelPosition* Bound, int numBound, 
	int* startIndex, int* endIndex);
bool GetIndexOfBPInZRange(
	short z1, short z2, 
	VoxelPosition* Bound, int numBound, 
	int* startIndex, int* endIndex);
bool GetIndexOfBPInYRange(
	short y1, short y2, 
	VoxelPosition* Bound, int numBound, 
	int startAt, int endAt,
	int* startIndex, int* endIndex);
bool GetIndexOfBPInXRange(
	short x1, short x2, 
	VoxelPosition* Bound, int numBound, 
	int startAt, int endAt,
	int* startIndex, int* endIndex);

bool SortBoundaryArray(int numBound, VoxelPosition Bound[]);
bool SortByX(int startAt, int endAt, VoxelPosition Bound[]);
bool SortByY(int startAt, int endAt, VoxelPosition Bound[]);
bool SortByZ(int startAt, int endAt, VoxelPosition Bound[]);


bool CalculatePotentialField(
	int L, int M, int N, 	      // [in] size of volume
	unsigned char* f, 	      // [in] volume flags
	int fieldStrenght,	      // [in] potential field strenght
	Vector* force,		      // [out] force field	
	bool inOut                    // [in] flag indicating that we don't 
	                              //    know what the inside/outside of 
	                              //    the object is. We have only point 
	                              //    samples of the boundary.
	                              //  DEFAULT: false (only interior)
) {
  int Lm1, Mm1, Nm1;
  int i,j,k, s, p;
  long idx, iidx, slsz, sz;
  VoxelPosition* Bound;
  int numBound = 0;
  bool flagSurf, flagBound;
  double r, t;
  int v1, v2, v3;
  int startIndex, tmpStartIndex, endIndex, tmpEndIndex, zStartIndex, zEndIndex, yStartIndex, yEndIndex;

  //
  // check volume padding - fast version
  //
  if(!CheckVolumePadding(f, L, M, N)) {
    printf("** Error - Object touches bounding box. Abort.\n");
    exit(1);
  }


#ifdef _DEBUG
	printf("\t************ Potential Field calculation parameters: ******************\n");
#ifdef HALF_BOUNDARY_POINTS
	printf("\t** Using only HALF of the boundary points.\n");
#else
	printf("\t** Using ALL boundary points.\n");
#endif

#ifdef EUCLIDEAN_METRIC
	printf("\t** Using EUCLIDEAN metric.\n");
#else
	printf("\t** Using NON EUCLIDEAN metric.\n");
#endif	
	if(inOut) {
	  printf("\t** Inside and Outside.\n");
	}
	else {
	  printf("\t** Inside ONLY.\n");
	}
	printf("\t********* Potential Field calculation parameters - end ****************\n");

#endif

  if((Bound = new VoxelPosition[BOUND_SIZE]) == NULL) {
	printf("\nERROR allocating memory for boundary array! - Abort\n");
	exit(1);
  }

  Lm1 = L - 1;
  Mm1 = M - 1;
  Nm1 = N - 1;
  slsz = L*M;		// slice size
  sz = slsz*N;

  // save all the boundary voxels in array Bound[]
	for (k = 1; k < Nm1; k++) {
		for (j = 1; j < Mm1; j++) {
			for (i = 1; i < Lm1; i++) {
				flagSurf = false;
				flagBound = true;
				idx = k*slsz + j*L + i;

				// CASE 1: treat the inner layer
				if (f[idx] == 0) continue;

				//consider six face neighbors, if anyone is zero, it is a boundary voxel
				iidx = k*slsz + j*L + i-1;
				if (f[iidx] == 0) {
					flagSurf = true;
				}
#ifdef HALF_BOUNDARY_POINTS
				// consider only half of the boundary points
				else {
					if (f[iidx] == BOUNDARY) {
						// a neighbour of the point was already selected so we will not select this one as part of the boundary.
						flagBound = false;
					}

				}
#endif

				if(!flagSurf || flagBound) {

				iidx = k*slsz + j*L + i+1;
				if (f[iidx] == 0) {
					flagSurf = true;
				}
#ifdef HALF_BOUNDARY_POINTS
				// consider only half of the boundary points
				else {
					if (f[iidx] == BOUNDARY) {
						// a neighbour of the point was already selected so we will not select it as part of the boundary.
						flagBound = false;
					}
				}
#endif

				if(!flagSurf || flagBound) {

				iidx = k*slsz + (j-1)*L + i;
				if (f[iidx] == 0) {
					flagSurf = true;
				}
#ifdef HALF_BOUNDARY_POINTS
				// consider only half of the boundary points
				else {
					if (f[iidx] == BOUNDARY) {
						// a neighbour of the point was already selected so we will not select it as part of the boundary.
						flagBound = false;
					}
				}
#endif

				if(!flagSurf || flagBound) {

				iidx = k*slsz + (j+1)*L + i;
				if (f[iidx] == 0) {
					flagSurf = true;
				}
#ifdef HALF_BOUNDARY_POINTS
				// consider only half of the boundary points
				else {
					if (f[iidx] == BOUNDARY) {
						// a neighbour of the point was already selected so we will not select it as part of the boundary.
						flagBound = false;
					}
				}
#endif

				if(!flagSurf || flagBound) {

				iidx = (k-1)*slsz + j*L + i;
				if (f[iidx] == 0) {
					flagSurf = true;
				}
#ifdef HALF_BOUNDARY_POINTS
				// consider only half of the boundary points
				else {
					if (f[iidx] == BOUNDARY) {
						// a neighbour of the point was already selected so we will not select it as part of the boundary.
						flagBound = false;
					}
				}
#endif

				if(!flagSurf || flagBound) {

				iidx = (k+1)*slsz + j*L + i;
				if (f[iidx] == 0) {
					flagSurf = true;
				}
#ifdef HALF_BOUNDARY_POINTS
				// consider only half of the boundary points
				else {
					if (f[iidx] == BOUNDARY) {
						// a neighbour of the point was already selected so we will not select it as part of the boundary.
						flagBound = false;
					}
				}
#endif

				}
				}
				}
				}
				}

				// restore idx to the right value
				idx = k*slsz + j*L + i;
				if (flagSurf) {
					f[idx] = SURF;

					if(flagBound) {
							// if no neighbour of this voxel is already marked as boundary, then mark this one.
							// or if we are taking all the boundary voxels 
							// 	(in this case flagBound stays true)
						f[idx] = BOUNDARY;
						Bound[numBound].x = i;
						Bound[numBound].y = j;
						Bound[numBound].z = k;
						numBound++;
						if(numBound >= BOUND_SIZE) {
							printf("ERROR: too many boundary points detected !! - Abort.\n");
							exit(1);
						}
					}
				}
			}
		}
	}

	//printf("numBound = %d \n", numBound);

#ifdef _DEBUG
	PrintElapsedTime("\tPF-1: finding the boundary voxels.");
	printf("\t--Found %d boundary voxels.\n", numBound);
#endif

/*
// output boundary voxels
FILE *ff;
unsigned char a;
long int b;

ff=fopen("bound.vol", "w");
for(idx=0; idx < L*M*N; idx++) {
if(f[idx] == BOUNDARY) {
a = 255;
}
else {
a = 0;
}


  b = random();
  if(b > RAND_MAX/2) {
  a = 0;
  }


  fwrite(&a, sizeof(unsigned char), 1, ff);
  b = 0;
  }
  fclose(ff);
  exit(1);
*/
	
// sort the boundary array.
SortBoundaryArray(numBound, Bound);

#ifdef _DEBUG
	PrintElapsedTime("\tPF-2: sorting the boundary voxels.");
#ifdef TRACE
	// print the boundary voxels
	for(i=0; i < numBound; i++) {
		printf("%d %d %d 0.5\n", Bound[i].x, Bound[i].y, Bound[i].z);
	}
	exit(1);
#endif	

#endif


// Compute the potential field
	printf("Computing potential field.\n");
idx = -1;
for (k = 0; k < N; k++) {

	printf("\tProcessing plane %d out of %d\r", k, N-1);
	fflush(stdout);

	// find the boundary voxels that will influence this point
	// look at the Z coordinate
	zStartIndex = 0;
	zEndIndex = numBound- 1;
	for (s = 0; s < numBound; s++) {
		if((k - Bound[s].z) <= PF_THRESHOLD) {
			zStartIndex = s;
			break;
		}
	}
	for (s = numBound-1; s >= zStartIndex; s--) {
		if((Bound[s].z - k) <= PF_THRESHOLD) {
			zEndIndex = s;
			break;
		}
	}
	// printf("ZStart: %d\t ZEnd: %d\n", zStartIndex, zEndIndex);

	for (j = 0; j < M; j++) {
		// find the boundary voxels that will influence this point
		// look at the Y coordinate
		yStartIndex = zStartIndex;
		yEndIndex = zEndIndex;
		for (s = zStartIndex; s <= zEndIndex; s++) {
			if((j - Bound[s].y) <= PF_THRESHOLD) {
				yStartIndex = s;
				break;
			}
		}
		for (s = zEndIndex; s >= yStartIndex; s--) {
			if((Bound[s].y - j) <= PF_THRESHOLD) {
				yEndIndex = s;
				break;
			}
		}
		// printf("YStart: %d\t YEnd: %d\n", yStartIndex, yEndIndex);

		for (i = 0; i < L; i++) {
			// printf("Point: %d\t%d\t%d:\n", i, j, k);
			// idx = k*slsz + j*L + i;
			idx = idx + 1;

			force[idx].xd = 0.00;
			force[idx].yd = 0.00;
			force[idx].zd = 0.00;

			if(!inOut) {
			  if(f[idx] == 0) {
			    // outside voxels have null force
			    continue;
			  }
			}
			else {
			  // we don't know where the inside of the object is
			  // so we compute the vector field everywhere.
			  // NOTHING
			}
			    
			// surface voxels (including those selected for the 
			//    field calculation)
			//	are ignored for now. The force there will be 
			//    the average of their neighbors
			//	if we are to compute the force at boundary 
			//    voxels too, the force will point
			//	towards the exterior of the object 
			//    (example: a 30x30x100 box)

			if(f[idx] == SURF) continue;
			if(f[idx] == BOUNDARY) continue;
			
			// find the boundary voxels that will influence this point
			// look at the X coordinate
			startIndex = yStartIndex;
			endIndex = yEndIndex;
			for (s = yStartIndex; s <= yEndIndex; s++) {
				if((i - Bound[s].x) <= PF_THRESHOLD) {
					startIndex = s;
					break;
				}
			}
			for (s = yEndIndex; s >= startIndex; s--) {
				if((Bound[s].x - i) <= PF_THRESHOLD) {
					endIndex = s;
					break;
				}
			}

			// printf("Start at: %d, end at: %d\n", startIndex, endIndex);
			// exit(-1);

			if(endIndex < startIndex) {
				// no boundary point is close enough to this point - take all the boundary points
				startIndex = 0;
				endIndex = numBound - 1;
			}

			for (s = startIndex; s <= endIndex; s++) {
				// printf("%d %d %d\n", i - Bound[s].x, j - Bound[s].y, k - Bound[s].z);

				/*
				// visibility test - too slow
				if(GetIndexOfBPInXYZRange(i, j, k,	Bound[s].x, Bound[s].y, Bound[s].z, 
						Bound, numBound, &v1, &v2))
				{
					// check if this boundary pont is visible from the current position
					if (IsLineCrossingBoundary(i, j, k, Bound[s].x, Bound[s].y, Bound[s].z, L, M, N, f)) {
						// not visible
						continue;
					}
				}
				*/
				
				v1 = i - Bound[s].x;
				v2 = j - Bound[s].y;
				v3 = k - Bound[s].z;
#ifdef EUCLIDEAN_METRIC
				// euclidean metric
				r = sqrt(v1*v1 + v2*v2 + v3*v3);
#else
				// simpler metric
				r = abs(v1) + abs(v2) + abs(v3);
#endif

				// r CAN BE 0 if we are computing the force 
				//    at boundary voxels too
				// if the current point is a BOUNDARY point, 
				//    some r will be 0, and that should be 
				//    ignored
				if(r != 0.00) {
				
				  // raise r to the fieldStrenght+1 power
				  //   so that the force is 
				  //   1/(dist^fieldStrength)
				  t = 1.00;
				  for(p = 0; p <= fieldStrenght; p++) {
				    t = t * r;
				  }
				  r = t;
				  
				  force[idx].xd = force[idx].xd + (v1 / r);
				  force[idx].yd = force[idx].yd + (v2 / r);
				  force[idx].zd = force[idx].zd + (v3 / r);
				}

			}
/*
			printf("First point with force vector != 0\n");
			printf("%f\t%f\t%f: %d, %d, %d\n", force[idx].xd, force[idx].yd, force[idx].zd, i, j, k);
			exit(1);
*/

		}
	}
}

// delete the Bound array - don't need it anymore
delete [] Bound;


#ifdef _DEBUG
	PrintElapsedTime("\tPF-3: computing potential field for inside voxels.");
#endif

// normalize force vectors:
for(idx=0; idx < L*M*N; idx++) {
  if(!inOut) {
    // only for interior voxels we had calculated forces
    if(f[idx] == EXTERIOR) continue;
  }
  
  r = force[idx].xd*force[idx].xd + 
      force[idx].yd*force[idx].yd + 
      force[idx].zd*force[idx].zd;
    
  if(r > 0.00) {
    r = sqrt(r);
    
    force[idx].xd = force[idx].xd / r;
    force[idx].yd = force[idx].yd / r;
    force[idx].zd = force[idx].zd / r;
  }
}

#ifdef _DEBUG
  PrintElapsedTime("\tPF-4: normalizing force vectors for inside voxels.");
#endif

  // if we know the inside from the outside
  // calculate the force at the surface voxels as the average of the 
  // interior neighbors
  if (!inOut) {
    //neighbors:
    int ng[26];
    
    // face neighbors
    ng[0]	= + slsz + 0 + 0;
    ng[1]	= - slsz + 0 + 0;
    ng[2]	= +    0 + L + 0;
    ng[3]	= +    0 - L + 0;
    ng[4]	= +    0 + 0 + 1;
    ng[5]	= +    0 + 0 - 1;
    // v-neighbors
    ng[6]	= - slsz - L - 1;
    ng[7]	= - slsz - L + 1;
    ng[8]	= - slsz + L - 1;
    ng[9]	= - slsz + L + 1;
    ng[10]	= + slsz - L - 1;
    ng[11]	= + slsz - L + 1;
    ng[12]	= + slsz + L - 1;
    ng[13]	= + slsz + L + 1;
    // e-neighbors
    ng[14]	= + slsz + L + 0;
    ng[15]	= + slsz - L + 0;
    ng[16]	= - slsz + L + 0;
    ng[17]	= - slsz - L + 0;
    ng[18]	= + slsz + 0 + 1;
    ng[19]	= + slsz + 0 - 1;
    ng[20]	= - slsz + 0 + 1;
    ng[21]	= - slsz + 0 - 1;
    ng[22]	= +    0 + L + 1;
    ng[23]	= +    0 + L - 1;
    ng[24]	= +    0 - L + 1;
    ng[25]	= +    0 - L - 1;
    
    for (k = 1; k < Nm1; k++) {
      for (j = 1; j < Mm1; j++) {
	for (i = 1; i < Lm1; i++) {
	  
	  idx = k*slsz + j*L + i;
	  
	  if((f[idx] == SURF) ||
	     (f[idx] == BOUNDARY))
	    {
	      force[idx].xd = 0.00;
	      force[idx].yd = 0.00;
	      force[idx].zd = 0.00;
	      
	      // look at the neighbors and average the forces if not 0
	      //
	      v1 = 0;
	      for(s=0; s < 26; s++) {
		
		iidx = idx + ng[s];		// index of neighbor
		// take only neighbors that are not SURF or BOUNDARY 
		//	because those neighbors have force = 0 
		if(f[iidx] == SURF)		continue;
		if(f[iidx] == BOUNDARY)	continue;
		
		// if we know the interior of the object, take only interior
		// neighbors
		if(!inOut) {
		  if(f[iidx] == EXTERIOR)	continue;
		}
		
		force[idx].xd = force[idx].xd + force[iidx].xd;
		force[idx].yd = force[idx].yd + force[iidx].yd;
		force[idx].zd = force[idx].zd + force[iidx].zd;
		v1 = v1 + 1;
		
	      }
	      
	      // average
	      if(v1 != 0) {
		force[idx].xd = force[idx].xd / (double) v1;
		force[idx].yd = force[idx].yd / (double) v1;
		force[idx].zd = force[idx].zd / (double) v1;
	      }
	      else {
		printf("Boundary voxel has no interior neighbor !!! - Force = 0\n");
	      }
	      
	      // normalize
	      r = force[idx].xd*force[idx].xd + 
		force[idx].yd*force[idx].yd + 
		force[idx].zd*force[idx].zd;
	      
	      if(r > 0.00) {
		r = sqrt(r);
		
		force[idx].xd = force[idx].xd / r;
		force[idx].yd = force[idx].yd / r;
		force[idx].zd = force[idx].zd / r;	
	      }
	    }
	}
      }
    }
  }
  else {
    // we don't know the inside from the outside.
    // boundary points remain 0
    // nothing to do
  }

#ifdef _DEBUG
	PrintElapsedTime("\tPF-5: computing potential field for boundary voxels.");
#endif


  return true;
}


// Sort the boundary array so that we can speed up the potential field calculation: ZYX in that order
// selection sort
bool SortBoundaryArray(int numBound, VoxelPosition Bound[]) {
	int st, i;
	short zvst, yvst;

	// sort by Z
	SortByZ(0, numBound-1, Bound);

	// then by Y
	st = 0;
	zvst = Bound[st].z;
	for(i=0; i < numBound; i++) {
		if(Bound[i].z != zvst) {
			SortByY(st, i-1, Bound);

			st = i;
			zvst = Bound[st].z;
		}
	}
	SortByY(st, numBound-1, Bound);

	// then by X
	st = 0;
	zvst = Bound[st].z;
	yvst = Bound[st].y;
	for(i=0; i < numBound; i++) {
		if((Bound[i].y != yvst) || (Bound[i].z != zvst)) {
			SortByX(st, i-1, Bound);

			st = i;
			zvst = Bound[st].z;
			yvst = Bound[st].y;
		}
	}
	SortByX(st, numBound-1, Bound);

	return true;
}

bool SortByX(int startAt, int endAt, VoxelPosition Bound[]) {
	int i, j, minIndex, crtMin;
	short tmp;

	for(i=startAt; i <= endAt; i++) {
		minIndex = -1;
		crtMin = Bound[i].x;

		for(j=i+1; j <= endAt; j++) {
			if(Bound[j].x < crtMin) {
				minIndex = j;
				crtMin = Bound[j].x;
			}
		}
		if(minIndex != -1) {
			// swap values.
			tmp = Bound[i].x;
			Bound[i].x = Bound[minIndex].x;
			Bound[minIndex].x = tmp;

			tmp = Bound[i].y;
			Bound[i].y = Bound[minIndex].y;
			Bound[minIndex].y = tmp;

			tmp = Bound[i].z;
			Bound[i].z = Bound[minIndex].z;
			Bound[minIndex].z = tmp;
		}
	}

	return true;
}

bool SortByY(int startAt, int endAt, VoxelPosition Bound[]) {
	int i, j, minIndex, crtMin;
	short tmp;

	for(i=startAt; i <= endAt; i++) {
		minIndex = -1;
		crtMin = Bound[i].y;

		for(j=i+1; j <= endAt; j++) {
			if(Bound[j].y < crtMin) {
				minIndex = j;
				crtMin = Bound[j].y;
			}
		}
		if(minIndex != -1) {
			// swap values.
			tmp = Bound[i].x;
			Bound[i].x = Bound[minIndex].x;
			Bound[minIndex].x = tmp;

			tmp = Bound[i].y;
			Bound[i].y = Bound[minIndex].y;
			Bound[minIndex].y = tmp;

			tmp = Bound[i].z;
			Bound[i].z = Bound[minIndex].z;
			Bound[minIndex].z = tmp;
		}
	}

	return true;
}


bool SortByZ(int startAt, int endAt, VoxelPosition Bound[]) {
	int i, j, minIndex, crtMin;
	short tmp;

	for(i=startAt; i <= endAt; i++) {
		minIndex = -1;
		crtMin = Bound[i].z;

		for(j=i+1; j <= endAt; j++) {
			if(Bound[j].z < crtMin) {
				minIndex = j;
				crtMin = Bound[j].z;
			}
		}
		if(minIndex != -1) {
			// swap values.
			tmp = Bound[i].x;
			Bound[i].x = Bound[minIndex].x;
			Bound[minIndex].x = tmp;

			tmp = Bound[i].y;
			Bound[i].y = Bound[minIndex].y;
			Bound[minIndex].y = tmp;

			tmp = Bound[i].z;
			Bound[i].z = Bound[minIndex].z;
			Bound[minIndex].z = tmp;
		}
	}

	return true;
}



// returns the start and endindex of boundary point found in a region
//	in space bound by a box defined by the 2 points.
// it doesn't change startIndex or endIndex if it returns false;
// returns true if it finds any boundary point in that region, or false otherwise.
bool GetIndexOfBPInXYZRange(
	short sx, short sy, short sz,
	short ex, short ey, short ez,
	VoxelPosition* Bound, int numBound, 
	int* startIndex, int* endIndex)
{
	int si1, ei1, si2, ei2;	// temporary start and end indexes
	// 
	if(GetIndexOfBPInZRange(sz, ez, Bound, numBound, &si1, &ei1)) {
		if(GetIndexOfBPInYRange(sy, ey, Bound, numBound, si1, ei1, &si2, &ei2)) {
			if(GetIndexOfBPInXRange(sx, ex, Bound, numBound, si2, ei2, &si1, &ei1)) {
				(*startIndex) = si1;
				(*endIndex) = ei1;
				return true;
			}
		}
	}
	return false;
}

bool GetIndexOfBPInZRange(
	short z1, short z2, 
	VoxelPosition* Bound, int numBound, 
	int* startIndex, int* endIndex)
{
	short minz, maxz;
	int s;
	int si;
	
	// sort the 2 z values;
	if(z1 < z2) {
		minz = z1; maxz = z2;
	}
	else {
		minz = z2; maxz = z1;
	}
	
	si 	= -1;
	for (s = 0; s < numBound; s++) {
		if((minz - Bound[s].z) < 0) {
			si = s;
			break;
		}
	}
	
	if(si == -1) {
		// couldn't find any boundary voxel
		return false;
	}
	
	(*startIndex) = si;
	
	for (s = numBound-1; s >= (*startIndex); s--) {
		if((Bound[s].z - maxz) < 0) {
			(*endIndex) = s;
			break;
		}
	}
	
	return true;
}
	
bool GetIndexOfBPInYRange(
	short y1, short y2, 
	VoxelPosition* Bound, int numBound, 
	int startAt, int endAt,
	int* startIndex, int* endIndex)
{
	short miny, maxy;
	int s;
	int si;
	
	// sort the 2 y values;
	if(y1 < y2) {
		miny = y1; maxy = y2;
	}
	else {
		miny = y2; maxy = y1;
	}
	
	// start the search at startAt and end it endAt
	si 	= -1;
	for (s = startAt; s <= endAt; s++) {
		if((miny - Bound[s].y) < 0) {
			si = s;
			break;
		}
	}
	
	if(si == -1) {
		// couldn't find any boundary voxel
		return false;
	}
	
	(*startIndex) = si;
	for (s = endAt; s >= (*startIndex); s--) {
		if((Bound[s].y - maxy) < 0) {
			(*endIndex) = s;
			break;
		}
	}
	
	return true;
}

bool GetIndexOfBPInXRange(
	short x1, short x2, 
	VoxelPosition* Bound, int numBound, 
	int startAt, int endAt,
	int* startIndex, int* endIndex)
{
	short minx, maxx;
	int s;
	int si;
	
	// sort the 2 x values;
	if(x1 < x2) {
		minx = x1; maxx = x2;
	}
	else {
		minx = x2; maxx = x1;
	}
	
	// start the search at startAt and end it endAt
	si 	= -1;
	for (s = startAt; s <= endAt; s++) {
		if((minx - Bound[s].x) < 0) {
			si = s;
			break;
		}
	}
	
	if(si == -1) {
		// couldn't find any boundary voxel
		return false;
	}
	
	(*startIndex) = si;
	for (s = endAt; s >= (*startIndex); s--) {
		if((Bound[s].x - maxx) < 0) {
			(*endIndex) = s;
			break;
		}
	}
	
	return true;
}


