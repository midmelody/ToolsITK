#include "common.h"


#ifndef WIN32
	struct timeval tv;
	unsigned long startTime, crtTime, prevTime;
	unsigned long phaseCompletionTime, elapsedTime;
#else
	DWORD	startTime, crtTime, prevTime;
	DWORD	phaseCompletionTime, elapsedTime;
#endif

void SetStartTime() {
	#ifdef WIN32
		startTime = GetTickCount();
	#else
		gettimeofday(&tv, NULL);
		startTime = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
	#endif
	prevTime = startTime;
	return;
}

void PrintElapsedTime(const char* message) {
	#ifdef WIN32
		crtTime = GetTickCount();
	#else
		gettimeofday(&tv, NULL);
		crtTime = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
	#endif
	phaseCompletionTime = crtTime - prevTime;
	elapsedTime = crtTime - startTime;
	prevTime = crtTime;

	if(strlen(message) == 0) {
		printf("Total time: %d ms.\n", elapsedTime);
	}
	else {
		if(strcmp(message, " ") == 0) {
			printf("completed in: %d ms. Total time: %d ms.\n", phaseCompletionTime, elapsedTime);
		}
		else {
			printf("%s\n"
				"\tcompleted in: %d ms. Total time: %d ms.\n", message, phaseCompletionTime, elapsedTime);
		}
	}
	return;
}



// for the following function to give accurate results, the following conditions must apply:
//	1. the volume has no holes in it
//	2. the object is padded by at least one plane of empty voxels in every direction
//
bool SetFlags(unsigned char *vol, int L, int M, int N, unsigned char **flags) {
	long idx, idx2, slsz, upperlimit;
	int i, j, k, s;

	if(((*flags) = new unsigned char[L*M*N]) == NULL) {
		printf("Error allocating memory for the flags array !! Abort.\n");
		exit(1);
	}

	slsz = L*M;		// slice size
	long neighborhood[6] = {-1, +1, -L, +L, -slsz, +slsz};	// only face neighbors

	// 1
	// set flag to INTERIOR for all the voxels that have a non 0 value
	//	and to EXTERIOR for all the 0 voxels.
	upperlimit = L*M*N;
	for(i=0; i < upperlimit; i++) {
		(*flags)[i] = INTERIOR;
		if(vol[i] == 0) {
			(*flags)[i] = EXTERIOR;
		}
	}


	// 2
	// look at the INTERIOR voxels. If an INTERIOR voxel has an EXTERIOR neighbor,
	//		then it is a SURFace voxel
	// look only at the neighbors defined by neighborhood.

	for (k = 1; k < N-1; k++) {
		for (j = 1; j < M-1; j++) {
			for (i = 1; i < L-1; i++) {
			
				idx = k*slsz + j*L + i;
				if((*flags)[idx] == INTERIOR) {
					for(s=0; s < 6; s++) {
						idx2 = idx + neighborhood[s];
						if((*flags)[idx2] == EXTERIOR) {
							(*flags)[idx] = SURF;
							break;
						}
					}
				}
			}
		}
	}

	return true;
}



bool ReadVolume(char *filename, int L, int M, int N, unsigned char **vol) {
	FILE* fvol;

	if ((fvol = fopen(filename,"r")) == NULL) {
		printf("Cannot open input file %s\n", filename);
		exit(1);
	}

	if(((*vol) = new unsigned char[L*M*N]) == NULL) {
		printf("UPS! - Error allocating memory for the volume. Abort.\n");
		exit(1);
	}

	if ( fread((*vol), sizeof(unsigned char), L*M*N, fvol) < L*M*N) {
	    printf("UPS! - File size is not the same as volume size\n");
	    delete [] vol;
	    exit(1);
  	}

	fclose(fvol);
	return true;
}

bool IsLineCrossingBoundary(
	short x1, short y1, short z1,
	short x2, short y2, short z2,
	int sX, int sY, int sZ,
	unsigned char* flags)
{
	double t, step, tmp;
	double x, y, z;
	long slsz, idx;
	
	
	// if the 2 voxels are the same, return false;
	if(	(x1 == x2) && 
		(y1 == y2) && 
		(z1 == z2))
	{
		return false;
	}
	
	// calculate optimum step that will take us to another voxel each time we move
	step = 2.00;
	
	// step for the X direction
	t = abs(x2 - x1);
	if(t > 0.00) {
		tmp = 1.00 / t;
		if(tmp < step) {
			step = tmp;
		}
	}
	
	// step for the Y direction
	t = abs(y2 - y1);
	if(t > 0.00) {
		tmp = 1.00 / t;
		if(tmp < step) {
			step = tmp;
		}
	}
	// step for the Z direction
	t = abs(z2 - z1);
	if(t > 0.00) {
		tmp = 1.00 / t;
		if(tmp < step) {
			step = tmp;
		}
	}
	
#ifdef _DEBUG
	// the 2 voxels are not identical so there will be a step value of at most 1.00
	//	but just to make sure, I will check
	if(step > 1.00) {
		printf("Vizibility test: UPS - step value is > 1.00. Abort\n");
		exit(1);
	}
#endif
	
	slsz = sX*sY;		// slice size
	// sample the line between the 2 points at <step> intervals, and check each
	//	of the sample points whether they are in one of the "blank" voxels
	t = 0.00;
	while (t <= 1.00) {
		x = x1 + (x2 - x1)*t;
		y = y1 + (y2 - y1)*t;
		z = z1 + (z2 - z1)*t;
		
		// index of voxel in the flags array
		idx = ((int)z)*slsz + ((int)y)*sX + ((int)x);
#ifdef _DEBUG		
		// this should not happen, but just in case...
		if((idx < 0) || (idx > sX*sY*sZ)) {
			printf("Vizibility test: UPS - line gets out of volume bounds ! Abort\n");
			exit(1);
		}
#endif		
		// 
		if(flags[idx] == EXTERIOR) {
			return true;
		}
		
		t = t + step;
	}
	
	// line did not intersect any "blank" voxels
	return false;
}


//
// function SaveSkeleton - saves the skeleton to a file
//   mode - 0 (default) saves skeleton as points, with the segment specified 
//            for each point
//            format: X Y Z segment 0.5\n
//          1 saves the skeleton as line segments (2 points per line)
//            format: X1 Y1 Z1 X2 Y2 Z2 segment\n 
//  
bool SaveSkeleton(Skeleton *Skel, char *file, char mode /*=0*/) {
  int i, j;
  FILE *fskelout;

#ifdef TRACE
  printf("Starting Save Skeleton ...n");
  printf("Skeleton:\n");
  printf("\tnumPoints = %d; sizePoints = %d\n", Skel->numPoints, 
	 Skel->sizePoints);
  printf("\tnumSegments = %d; sizeSegments = %d\n", Skel->numSegments, 
	 Skel->sizeSegments);
  printf("Segments:\n");
  printf("\tLEFT\tFIRST\tLAST\tRIGHT\n");
  for(i=0; i < Skel->numSegments; i++) {
    printf("\t%d\t%d\t%d\t%d\n", Skel->Segments[i][SKEL_SEG_LEFT],
	   Skel->Segments[i][SKEL_SEG_FIRST],
	   Skel->Segments[i][SKEL_SEG_LAST],
	   Skel->Segments[i][SKEL_SEG_RIGHT]);
  }
  printf("-----\n");
  fflush(stdout);
#endif

// open the file
  if ((fskelout = fopen(file,"w")) == NULL) {
    printf("Cannot open output file %s for writing\n", file);
    exit(1);
  }

  switch(mode) {
  case 0:
    {
      //
      // write out the skeleton points
      //
      // printed is an array that specifies for each skeleton point whether it
      //   was printed to the output file or not
      // used because intersection points belong to more than one segment and 
      //   might be printed more than one time.

      bool *printed = NULL;
      printed = new bool[Skel->numPoints];
      if(printed == NULL) {
	printf("UPS! - Error allocating memory for the working data structures \
- Abort.\n");
	exit(1);
      }
      
      // initialize to false
      for(i=0; i < Skel->numPoints; i++) {
	printed[i] = false;
      }
      
      // output each segment
      for(i=0; i < Skel->numSegments; i++) {
	
	// output the left end point of the segment
	if(!printed[Skel->Segments[i][SKEL_SEG_LEFT]]) {
	  fprintf(fskelout,"%f %f %f\n", 
		  Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].x, 
		  Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].y, 
		  Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].z);
	  printed[Skel->Segments[i][SKEL_SEG_LEFT]] = true;
	}
	
	// output the interior points
	for(j = Skel->Segments[i][SKEL_SEG_FIRST]; 
	    j <=  Skel->Segments[i][SKEL_SEG_LAST]; j++) 
	  {
	    if( (j != Skel->Segments[i][SKEL_SEG_LEFT]) && 
		(j != Skel->Segments[i][SKEL_SEG_RIGHT]))
	      {
		if(!printed[j]) {
		  fprintf(fskelout,"%f %f %f\n", 
			  Skel->Points[j].x, 
			  Skel->Points[j].y, 
			  Skel->Points[j].z);
		  printed[j] = true;
		}
	      }
	  }
	
	// output the right endpoint
	if(!printed[Skel->Segments[i][SKEL_SEG_RIGHT]]) {
	  fprintf(fskelout,"%f %f %f\n", 
		  Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].x, 
		  Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].y, 
		  Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].z);
	  printed[Skel->Segments[i][SKEL_SEG_RIGHT]] = true;
	}
	
      }

    }
    break;
  case 1:
    {
      // output line segments
      for(i=0; i < Skel->numSegments; i++) {
	
	// output the left and right end points of the segment
	// and the segment number
	fprintf(fskelout,"%f %f %f %f %f %f %d\n", 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].x, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].y, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].z, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].x, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].y, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].z, 
		i);
      }
      break;
    }
  default:
    printf("Wrong parameter to SaveSkeleton: %d ! Skeleton was NOT saved.\n", 
	   mode);
    break;
  }

  // close the file
  fclose(fskelout);

  return true;
}


///////////////////////////////////////////////////////////////////////////////
bool ReadVectorField(Vector *field, int L, int M, int N, char *fileName) {
  FILE *finVF;
  int i, j, k;
  long idx;

  // open the file
  finVF = fopen(fileName, "r");
  if (finVF == NULL)  {
    printf("Couldn't open vector field file for reading: %s. Abort\n", 
	   fileName);
    return false;
  }

  // read in force vectors
  for (k = 0; k < N; k++) {
    for (j = 0; j < M; j++) {
      for (i = 0; i < L; i++) {
	idx = k*L*M + j*L + i;
	if(fscanf(finVF, "%lf %lf %lf", 
		  &(field[idx].xd), &(field[idx].yd), &(field[idx].zd)) != 3) 
	{
	  printf("Error reading vector field input file. Abort\n");
	  return false;
	}
      }
    }
  }

  // close the file
  fclose(finVF);

  return true;
}


///////////////////////////////////////////////////////////////////////////////
bool SaveVectorField(Vector *field, int L, int M, int N, char *fileName) {
  FILE *foutVF;
  int i, j, k;
  long idx;
  
  // open the file
  if ((foutVF = fopen(fileName,"w")) == NULL) {
    printf("Cannot open output file %s for writing\n", fileName);
    return false;
  }
  
  // write force vectors
  for (k = 0; k < N; k++) {
    for (j = 0; j < M; j++) {
      for (i = 0; i < L; i++) {
	idx = k*L*M + j*L + i;
	fprintf(foutVF, "%lf %lf %lf\n", 
		field[idx].xd, field[idx].yd, field[idx].zd);
      }
    }
  }

  // close the file
  fclose(foutVF);

  return true;
}


///////////////////////////////////////////////////////////////////////////////
// for volume vol, sets the voxel values to SURF, INTERIOR or EXTERIOR
///////////////////////////////////////////////////////////////////////////////
bool FlagVolume(unsigned char *vol, int L, int M, int N) {
  long idx, idx2, slsz, upperlimit;
  int i, j, k, s;
  
  if(vol == NULL) {
    return false;
  }

  slsz = L*M;		// slice size
  long neighborhood[6] = {-1, +1, -L, +L, -slsz, +slsz};	
  // only face neighbors

  // 1
  // set flag to INTERIOR for all the voxels that have a non 0 value
  //	and to EXTERIOR for all the 0 voxels.
  upperlimit = L*M*N;
  for(i=0; i < upperlimit; i++) {
    if(vol[i] == 0) {
      vol[i] = EXTERIOR;
    }
    else {
      vol[i] = INTERIOR;
    }
  }


  // 2
  // look at the INTERIOR voxels. 
  // If an INTERIOR voxel has an EXTERIOR neighbor,
  //    then it is a SURFace voxel
  // Look only at the neighbors defined by neighborhood.

  for (k = 1; k < N-1; k++) {
    for (j = 1; j < M-1; j++) {
      for (i = 1; i < L-1; i++) {
	idx = k*slsz + j*L + i;
	if(vol[idx] == INTERIOR) {
	  for(s=0; s < 6; s++) {
	    idx2 = idx + neighborhood[s];
	    if(vol[idx2] == EXTERIOR) {
	      vol[idx] = SURF;
	      break;
	    }
	  }
	}
      }
    }
  }

  return true;
}


// 
// Function CheckVolumePadding 
// Make sure the volume is padded with at lease 1 empty plane 
//    in all 3 directions
//
bool CheckVolumePadding(
        unsigned char *vol, 	 // [in] volume to be checked
	int L, int M, int N   // [in] volume size (X, Y and Z). 
) {
  
  long slsz, idx;
  int i, j, k;

  slsz = L*M;
  
  k = 0;
  for (j = 0; j < M; j++) {
    for (i = 0; i < L; i++) {
      idx = k*slsz + j*L + i;
      if(vol[idx] != EXTERIOR) {
	return false;
      }
    }
  }
  k = N-1;
  for (j = 0; j < M; j++) {
    for (i = 0; i < L; i++) {
      idx = k*slsz + j*L + i;
      if(vol[idx] != EXTERIOR) {
	return false;
      }
    }
  }
  
  j = 0;
  for (k = 0; k < N; k++) {
    for (i = 0; i < L; i++) {
      idx = k*slsz + j*L + i;
      if(vol[idx] != EXTERIOR) {
	return false;
      }
    }
  }
  j = M-1;
  for (k = 0; k < N; k++) {
    for (i = 0; i < L; i++) {
      idx = k*slsz + j*L + i;
      if(vol[idx] != EXTERIOR) {
	return false;
      }
    }
  }

  i = 0;
  for (k = 0; k < N; k++) {
    for (j = 0; j < M; j++) {
      idx = k*slsz + j*L + i;
      if(vol[idx] != EXTERIOR) {
	return false;
      }
    }
  }
  i = L-1;
  for (k = 0; k < N; k++) {
    for (j = 0; j < M; j++) {
      idx = k*slsz + j*L + i;
      if(vol[idx] != EXTERIOR) {
	return false;
      }
    }
  }

  return true;
}
