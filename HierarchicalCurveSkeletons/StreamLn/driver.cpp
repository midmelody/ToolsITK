//
// driver program for the StreamLn module
//

#include "../common.h"
#include "StreamLn.h"

// #define TRACE

#define MAX_NUM_HDPTS		1000

main(int argc, char *argv[])
{
  FILE *fout, *finCP, *fseeds, *finVF;
  char *infilename, *seedsfilename;
  int L,M,N;         // Sizes in x,y,z dimensions

  int i,j,k;
  Vector *fv;

  VoxelPositionDouble *HDPoints;
  VoxelPosition *BdSeeds;
  Skeleton *Skel;
  int numBoundSeeds, numHDPoints;

  long idx, slsz;

  unsigned char *cf;	// the volume
  unsigned char *f;		// flags

  char line[512];
  CriticalPoint CritPts[1000];
  int numCritPts;

  double v1, v2, v3;
  CriticalPointType tp;
  int read;

#ifdef _DEBUG
  SetStartTime();
#endif

  if (argc < 8)
  {
    printf("\n\
Usage: %s <vector file> <xs> <ys> <zs> <critical points file> \n\
          <seeds file|-noseeds> <out skel>.\n",argv[0]);
    printf("\tif -noseeds is specified intead of the seeds file name, \n\tno seeds are read.\n");
    exit(1);
  }


  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);
  slsz = L*M;		// slice size


  if((fv = new Vector[L*M*N]) == NULL) {
    printf("UPS! - Error allocating memory for the input vector field. \
Abort.\n");
    exit(1);
  }

  // read in force vectors
  ReadVectorField(fv, L, M, N, argv[1]);

#ifdef _DEBUG
  PrintElapsedTime("Phase2: reading in force vectors");
#endif

  float f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15;

// reading in critical points
  numCritPts = 0;
  
  if ((finCP = fopen(argv[5],"r")) == NULL) {
    printf("Cannot open %s for reading\n",argv[5]);
    exit(1);
  }
  
  while(!feof(finCP)) {
    if(fgets(line, 500, finCP) != NULL) {
      //printf("%s\n", line);
      
      if(sscanf(line, "%f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f", 
		&f1, &f2, &f3, &CritPts[numCritPts].type, 
		&f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13, &f14, &f15) == 16) 
	{	
	  CritPts[numCritPts].position.x = f1;
	  CritPts[numCritPts].position.y = f2;
	  CritPts[numCritPts].position.z = f3; 
	  CritPts[numCritPts].eval[0] = f4;
	  CritPts[numCritPts].eval[1] = f5;
	  CritPts[numCritPts].eval[2] = f6;
	  CritPts[numCritPts].evect[0].xd = f7;
	  CritPts[numCritPts].evect[0].yd = f8;
	  CritPts[numCritPts].evect[0].zd = f9;
	  CritPts[numCritPts].evect[1].xd = f10;
	  CritPts[numCritPts].evect[1].yd = f11;
	  CritPts[numCritPts].evect[1].zd = f12;
	  CritPts[numCritPts].evect[2].xd = f13;
	  CritPts[numCritPts].evect[2].yd = f14;
	  CritPts[numCritPts].evect[2].zd = f15;
	  
	  numCritPts++;
	  if(numCritPts >= 1000) {
	    printf("Too many critical points read. Maximum 1000. Abort\n");
	    exit(1);
	  }
	}
      else {
	printf("Error reading critical points file. - Abort\n");
	exit(1);
      }
    }
  }	

  fclose(finCP);

#ifdef _DEBUG
  PrintElapsedTime("Phase3: reading in critical points");
  printf("** Read %d critical points\n", numCritPts);
  int nrs, nra, nrr;
  nrs = 0;	nra = 0;	nrr = 0;
  for(i=0; i < numCritPts; i++) {
    
    switch(CritPts[i].type) {
    case CPT_SADDLE:
      nrs++;
      break;
    case CPT_ATTRACTING_NODE:
      nra++;
      break;
    case CPT_REPELLING_NODE:
      nrr++;
      break;
    }
  }
  printf("%d saddles, %d attracting nodes, %d repelling nodes.\n", nrs, nra, nrr);
#endif

// reading seed points
  numHDPoints = 0;
  HDPoints = NULL;
  
  if(strncmp(argv[6], "-noseeds", 8) == 0) {
    // do not open the seeds file
  }
  else {
    if ((fseeds = fopen(argv[6], "r")) == NULL)  {
      printf("Error: couldn't open %s for input\n", argv[6]);
      exit(1);
    }
    
    if((HDPoints = new VoxelPositionDouble[MAX_NUM_HDPTS]) == NULL) {
      printf("UPS! - Error allocating memory for the seeds array! - Abort.\n");
      exit(1);
    }
    
    numHDPoints = 0;
    while(!feof(fseeds)) {
      read =  fscanf(fseeds, "%lf %lf %lf %d\n", 
		     &(HDPoints[numHDPoints].x),  
		     &(HDPoints[numHDPoints].y), 
		     &(HDPoints[numHDPoints].z), &j);
      if(read!= 4)
	{
	  if(read != EOF) {
	    printf("Error reading seeds file. Abort\n");
	    exit(1);
	  }
	  else {
	    printf("Reached end of file.\n");
	    break;
	  }
	}
      numHDPoints++;
    }
    
    fclose(fseeds);
  }

#ifdef _DEBUG
  PrintElapsedTime("Phase3: reading in seed points");
#endif

#ifdef TRACE	  
  printf("Read %d seed points:\n", numHDPoints);
  for(i=0; i < numHDPoints; i++) {
    printf("%f %f %f\n", HDPoints[i].x, HDPoints[i].y, HDPoints[i].z);
  }
#endif

  Skel = NULL;
  AllocateSkeleton(&Skel, 100000, 10000);
  
#ifdef TRACE
  printf("Skeleton structure allocated.\n");
#endif

  BdSeeds = NULL;
  numBoundSeeds = 0;

  GetStreamLines(fv, L, M, N, 
		 CritPts, numCritPts, 
		 BdSeeds, numBoundSeeds, 
		 HDPoints, numHDPoints, 
		 Skel);


#ifdef _DEBUG
  PrintElapsedTime("Phase4: Getting the streamlines\n");
#endif
  if(numBoundSeeds > 0) {
    delete [] BdSeeds;
  }
  if(numHDPoints > 0) {
    delete [] HDPoints;
  }

  SaveSkeleton(Skel, argv[7], 0);

  FreeSkeleton(&Skel);

  
#ifdef _DEBUG
  PrintElapsedTime("Phase5: Writing output file");
#endif
  printf("End \n");


}
