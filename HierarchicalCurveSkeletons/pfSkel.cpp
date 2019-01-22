//
// Computes the skeleton of a 3D object
//

#include "common.h"
#include "pfSkel.h"

// #define TRACE

bool CheckVolumePadding(unsigned char **vol, 
			int *L, int *M, int*N, 
			int distCharges);


bool pfSkel(
	char* volFileName, 		     // [in] volume file name
	int L, int M, int N, 		     // [in] volume size (x, y and z)
	int distCharges,		     // [in] distance from object 
	                                     //   boundary where to place the 
	                                     //   charges (>=0)
	int fieldStrenght,		     // [in] potential field strenght 
	                                     //   (4 .. 9)
	float pHDPts,		             // [in] percentage of high 
	                                     //   divergence points to use
	Skeleton *Skel, 	             // [out] array containing 
	                                     //   skeleton points
	cbChangeParameters pfnChangeParams,  // [in]  callback function to 
	                                     //   change parameters on the run.
	void *other,                         // [in] value to be passed on to
	                                     //    the callback function
	char *vfInFile,                      // [in] vector field input file
	char *vfOutFile                      // [in] vector field output file
){
  int i,j,k;
  unsigned char *f;
  Vector *force;
  long idx, sz, slsz;
  int xs, ys, zs;

  VoxelPosition *HCBPoints = NULL;
  int numHCBPoints = 0;

  CriticalPoint *CritPts = NULL;
  int numCritPts = 0;

  VoxelPositionDouble *HDPts = NULL;
  int numHDPts = 0;

  // this variable should be removed and added as a parameters of this function
  float percHCPts = 0;
  float percHDPts = pHDPts;

  int **segmentsLevel1;
  int numSegmentsLevel1;

  // distance of charges should be >= 0 but it is ignored if vfInFile 
  //   is not NULL
  if((distCharges < 0) && (vfInFile == NULL))  {
    printf("pfSkel.cpp: <distCharges> parameter cannot be < 0.\n");
    return false;
  }

  // field strength should be between 1 and 10 but it is ignored if vfInFile 
  //    is not NULL
  if(((fieldStrenght < 1) || (fieldStrenght > 10)) && (vfInFile == NULL)) {
    printf("pfSkel.cpp: <fieldStrenght> parameter must be in the [1..10] \
range. Recommended values [4..9].\n");
    return false;
  }
  
  // read in the volume
  // volume array f is allocated inside this function
#ifdef _DEBUG
  printf("PFS-1: Reading volume data...\n");
#endif

  if(!ReadVolume(volFileName, L, M, N, &f)) {
    return false;
  }

#ifdef _DEBUG
  PrintElapsedTime(" ");
#endif

  if(vfInFile == NULL) {
    //
    // Make sure the volume is padded by enough empty planes in all 3 
    // directions
    //    so that placing charges at specified distance from boundary will 
    // still leave a plane of empty voxels in each direction
    //
    if(!CheckVolumePadding(&f, &L, &M, &N, distCharges)) {
      printf("Volume needs to be padded by at least 1 layer of empty voxels in every direction.\n (The object should not touch the volume's bounding box)\n Program will now exit.\n");
      return false;
    }
  }
  
  sz = L*M*N;
  slsz = L*M;

#ifdef _DEBUG
  printf(" PFS-2: Make solid volume...\n ");
#endif

  // 
  // Make sure the volume does not have holes in it
  //
  MakeSolidVolume(L, M, N, f, EXTERIOR, INTERIOR); 

#ifdef _DEBUG 
  PrintElapsedTime(" ");
#endif

#ifdef _DEBUG
  printf("PFS-3: Calculating high curvature boundary points - NOT DONE.\n");
#endif

  //
  // detecting the high curvature boundary points, and retrieve the first 1%
  //
  // Module is not working yet, ...
  HCBPoints = NULL;
  numHCBPoints = 0;
  // GetHighCurvatureBoundaryPoints(cf, L, M, N, &HCBPoints, &numHCBPoints,	1.00);


  // flag the volume
  FlagVolume(f, L, M, N);


/*
// output the flags array
printf("output flags array\n");
if ((fskelout = fopen(argv[5],"w")) == NULL) {
	printf("Cannot open output file %s for writing\n",argv[5]);
    exit(1);
}

fwrite(f, sizeof(unsigned char), L*M*N, fskelout);
fclose(fskelout);
exit(1);
*/

  force = new Vector[L*M*N];			// potential field
  if(force == NULL) {
    printf("Error allocating memory for the working data structures. Abort\n");
    exit(1);
  }

  if(vfInFile == NULL) {
    //
    // thicken the object with <distCharges> layers of extra voxels and 
    //    compute the potential field 
    //
#ifdef _DEBUG
    printf("PFS-4: Placing charges outside the object (at %d layers)...\n", 
	   distCharges);
#endif
    
    ExpandNLayers(L, M, N, f, distCharges);
    
#ifdef _DEBUG
    PrintElapsedTime(" ");
#endif
    
    //
    // Calculating potential field...
    //
#ifdef _DEBUG
    printf("Phase 5: Calculating potential field...\n");
#endif
    
    // Compute force vectors
    CalculatePotentialField(L, M, N, f,	fieldStrenght, force);
    
#ifdef _DEBUG
    PrintElapsedTime(" ");
#endif
  }
  else {
    // do not compute the potential field, just read it from the specified file
    printf("Reading vector field from file %s...", vfInFile);
    ReadVectorField(force, L, M, N, vfInFile);
    printf("done.\n");
  }


  //
  // If vfOutFile id not null, dump the vector field to that file
  //
  if(vfOutFile != NULL) {
    SaveVectorField(force, L, M, N, vfOutFile);
  }

//
// Detecting critical points
//

#ifdef _DEBUG
  printf("Phase 6: Detecting critical points...\n");
#endif

  CritPts = NULL;
  numCritPts = 0;
  GetCriticalPoints(force, L, M, N, f, &CritPts, &numCritPts);

#ifdef _DEBUG
  PrintElapsedTime(" ");
#endif


  //
  // generating the skeleton
  //
#ifdef _DEBUG
  printf("Phase 7: Detecting skeleton points...\n");
#endif

  bool stop = false;
  int numPointsInLevel1Skeleton = 0;
  pfSkelCommand cmd;

  // in interactive mode, we allow the user to change parameters and then run
  // the streamlines part again, until he/she is satisfied with the result.

  // first we compute the level 1 skeleton that connects only critical points.
  // The user has no control over the computation of this part, so it is done
  //   only once and then reused in case the parameters change.
  //
  GetLevel1Skeleton(force, L, M, N, CritPts, numCritPts, Skel);

#ifdef _DEBUG
  printf("Phase 7.1: Level 1 skeleton generation completed.\n");
#endif
  //
  // save the level 1 skeleton - only need to save the segments array and the
  // save the number of skeleton points computed at this point 
  // so we can restart at this point if the parameters change
  //
  numPointsInLevel1Skeleton = Skel->numPoints;
  numSegmentsLevel1 = Skel->numSegments;
  // allocate space for saving the segments
  if(numSegmentsLevel1 > 0) {
    if((segmentsLevel1 = new int*[numSegmentsLevel1]) == NULL) {
      printf("Error allocating memory for the working data structures. \
Abort\n");
      exit(1);
    }
    
    for(i = 0; i < numSegmentsLevel1; i++) {
      if((segmentsLevel1[i] = new int[4]) == NULL) {
	printf("Error allocating memory for the working data structures. \
Abort\n");
	exit(1);
      }
    }
  }

  // save the segments array
  for(i=0; i < numSegmentsLevel1; i++) {
    for(j=0; j < 4; j++) {
      segmentsLevel1[i][j] = Skel->Segments[i][j];
    }
  }

  // compute level 2 skeleton using as basis the level 1 skeleton
  stop = false;
  while(!stop) {
    //
    // Get top ... % of high divergence points
    //
    HDPts = NULL;
    numHDPts = 0;
    GetHighDivergencePoints(force, L, M, N, f, percHDPts, &HDPts, &numHDPts);
    
    //
    // get the skeleton
    //
    GetLevel2Skeleton(force, L, M, N, HDPts, numHDPts, Skel);
    
    // at this point we give the user a chance to save the skeleton and
    // change the parameters if he/she wishes.
    //
    cmd.HDP = percHDPts;
    cmd.HCP = percHCPts;
    
    if((*pfnChangeParams)(Skel, &cmd, other)) {
      switch(cmd.cmdCode[0]) {
      case 'q':
	// quit
	stop = true;
	break;
      case 'p':
	// change parameters request
	percHDPts = cmd.HDP;
	percHCPts = cmd.HCP;
	break;
      default:
	// what could this be ???
	stop = true;
	break;
      }
    }
    else {
      stop = true;
    }

    printf("Deleting HD points array...\n");
    delete [] HDPts;
    HDPts = NULL;
    printf("done\n");

    // return to level1 skeleton
    Skel->numPoints = numPointsInLevel1Skeleton;
    for(i=0; i < numSegmentsLevel1; i++) {
      for(j=0; j < 4; j++) {
	Skel->Segments[i][j] = segmentsLevel1[i][j];
      }
    }
    Skel->numSegments = numSegmentsLevel1;
  }

  // free the memory wasted with saving the level 1 segments 
  if(numSegmentsLevel1 > 0) {
    for(i=0; i < numSegmentsLevel1; i++) {
      delete [] segmentsLevel1[i];
    }
    delete [] segmentsLevel1;
    segmentsLevel1 = NULL;
  }
  numSegmentsLevel1 = 0;

#ifdef _DEBUG
  PrintElapsedTime(" ");
#endif



/*
  //
  // THIS STEP IS NO LONGER PERFORMED FOR NOW ...
  //   Removing points from the skeleton might break some segments
  //      and we have to figure out what we should do in that case
  // last step:
  // removing skeleton points that are outside the original object
  // read in the volume again

#ifdef _DEBUG
  printf("PFS-8: Checking skeleton points...\n");
#endif

  ReadVolume(volFileName, L, M, N, &f);
  MakeSolidVolume(L, M, N, f, EXTERIOR, INTERIOR);
  j = 0;
  for(i=0; i < (*numSkel); i++) {
    xs = (int)(*Skel)[i].x;
    ys = (int)(*Skel)[i].y;
    zs = (int)(*Skel)[i].z;
    
    idx = zs*slsz + ys*L + xs;
    if(f[idx] == EXTERIOR) {
      // delete skeleton point
      printf("** Removing skeleton point (%d): (%lf, %lf, %lf)\n", i, (*Skel)[i].x, (*Skel)[i].y, (*Skel)[i].z);
      (*Skel)[i].x = -1.00;
      (*Skel)[i].y = -1.00;
      (*Skel)[i].z = -1.00;
      j++;
    }	
  }

  // clean up the gaps in the Skel array if any.
  if(j > 0) {
    k = 0;
    for(i=0; i < (*numSkel); i++) {
      if(	((*Skel)[i].x == -1.00) &&
		((*Skel)[i].y == -1.00) &&
		((*Skel)[i].z == -1.00)) 
	{
	  k++;
	  continue;
	}
      if(k > 0) {
	(*Skel)[i-k].x = (*Skel)[i].x;
	(*Skel)[i-k].y = (*Skel)[i].y;
	(*Skel)[i-k].z = (*Skel)[i].z;
      }
    }
    (*numSkel) = (*numSkel) - k;
  }
  
#ifdef _DEBUG
	PrintElapsedTime(" ");
#endif

#ifdef _DEBUG
	// printf("*** Removed %d skeleton points that are outside the object ???\n", j);
#endif
*/

/*
// write out the skeleton points
for(i=0; i < numSkel; i++) {
	printf("%f %f %f 0.5 1\n", Skel[i].x, Skel[i].y, Skel[i].z);
}
*/

/*
// write out the force vectors
for(idx=0; idx < L*M*N; idx++) {
	printf("%f %f %f\n",force[idx].xd,force[idx].yd, force[idx].zd);
}
*/

/*
// write out the critical points
for(i=0; i < numCritPts; i++) {
   printf("%f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f\n", 
   CritPts[i].position.x, CritPts[i].position.y, CritPts[i].position.z, 
   CritPts[i].type, 
   CritPts[i].eval[0], CritPts[i].eval[1], CritPts[i].eval[2], 
   CritPts[i].evect[0].xd, CritPts[i].evect[0].yd, CritPts[i].evect[0].zd,
   CritPts[i].evect[1].xd, CritPts[i].evect[1].yd, CritPts[i].evect[1].zd,
   CritPts[i].evect[2].xd, CritPts[i].evect[2].yd, CritPts[i].evect[2].zd);
}
*/


  // free the allocated memory
  //
  delete [] f;
  delete [] force;
  // delete [] HCBPoints;
  delete [] CritPts;

  // do not free the skeleton structure here. The calling program should do 
  //   that

  return true;
}
		


// 
// Function CheckVolumePadding 
// Make sure the volume is padded by enough empty planes in all 3 directions
//	so that placing charges at specified distance from boundary will still 
//    leave
//	a plane of empty voxels in each direction
//
bool CheckVolumePadding(
        unsigned char **vol, 	 // [in,out] volume to be checked
	int *L, int *M, int*N,   // [in,out] volume size (X, Y and Z). 
	                         //    In case the volume needs 
				 //    to be padded with more empty planes, 
	                         //    the size will change
	int distCharges		 // distance from boundary where charges 
	                         //    will be placed
) {

  int i, j, k;
  long idx, slsz;
  int minx, maxx, miny, maxy, minz, maxz;
  
  maxx = 0;		maxy = 0;		maxz = 0;
  minx = (*L);	miny = (*M);	minz = (*N);
  
  slsz = (*L)*(*M);		// slice size
  
  for (k = 0; k < (*N); k++) {
    for (j = 0; j < (*M); j++) {
      for (i = 0; i < (*L); i++) {
	
	idx = k*slsz + j*(*L) + i;
	if((*vol)[idx] != EXTERIOR) {
	  if(i < minx) minx = i;
	  if(j < miny) miny = j;
	  if(k < minz) minz = k;
	  
	  if(i > maxx) maxx = i;
	  if(j > maxy) maxy = j;
	  if(k > maxz) maxz = k;
	}
      }
    }
  }
#ifdef TRACE	
  printf("CheckVolumePadding: Volume extents: (%d, %d, %d) - (%d, %d, %d)\n", minx, miny, minz, maxx, maxy, maxz);
#endif
  if((minx - distCharges) <= 0) return false;
  if((miny - distCharges) <= 0) return false;
  if((minz - distCharges) <= 0) return false;
  if((maxx + distCharges) >= (*L))	return false;
  if((maxy + distCharges) >= (*M))	return false;
  if((maxz + distCharges) >= (*N))	return false;
  
  return true;
}
