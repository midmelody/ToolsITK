//
// driver program for the CritPts module
//

#include "../common.h"
#include "CritPts.h"

#include <fstream>
#include <iostream>

int main(int argc, char *argv[])
{
  FILE *fin, *fseeds, *fvol;
  FILE *fout;
  char *infilename;
  int L,M,N;         // Sizes in x,y,z dimensions

  int i,j,k;
  Vector *fv;

  CriticalPoint *CritPts;
  int numCritPts;

  long idx, slsz;

  unsigned char *f;		// flags
  bool inOut = false;
  
  int sd, at, re;

#ifdef _DEBUG
  SetStartTime();
#endif

  if (argc < 7)
  {
    printf("\n\
Usage: %s <volume file> <xs> <ys> <zs> <vector file> \n\
          <critical points file> [<inOut>].\n\
  <inOut> flag specifies whether you want to look for critical points \n\
          outside the object too (set to 1), \n\
          not just in the interior (set to 0).\n\
          Defalut value is 0 (inside only).\n\n",
	   argv[0]);
    exit(1);
  }


  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);
  if(argc > 7) {
    if(atoi(argv[7]) == 0) {
      inOut = false;
    }
    else {
      inOut = true;
    }
  }

  slsz = L*M;		// slice size

  // read volume
  ReadVolume(argv[1], L, M, N, &f);

  // set flags
  FlagVolume(f, L, M, N);

#ifdef _DEBUG
  PrintElapsedTime("Phase1: reading volume and setting flags");
#endif

  if((fv = new Vector[L*M*N]) == NULL) {
    printf("UPS! - Error allocating memory for the input vector field. Abort.\n");
    exit(1);
  }

  // read in force vectors
  ReadVectorField(fv, L, M, N, argv[5]);


#ifdef _DEBUG
  PrintElapsedTime("Phase2: reading in force vectors");
#endif

  CritPts = NULL;
  numCritPts = 0;

  GetCriticalPoints(fv, L, M, N, f, &CritPts, &numCritPts, inOut);

#ifdef _DEBUG
  PrintElapsedTime("Phase4: Getting the critical points.");
  printf("**  %d critical points returned.\n", numCritPts);
#endif
  
  delete [] f;

  if ((fout = fopen(argv[6],"w")) == NULL) {
    printf("Cannot open %s for writing\n",argv[6]);
    exit(1);
  }

  sd = 0;
  at = 0;
  re = 0;
  
  for(i=0; i < numCritPts; i++) {
#ifdef TRACE
    printf("writing point (%f, %f, %f) to output file.\n", 
	   CritPts[i].position.x, CritPts[i].position.y, CritPts[i].position.z);
#endif

    fprintf(fout,"%f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f\n", 
	    CritPts[i].position.x, CritPts[i].position.y, CritPts[i].position.z, 
	    CritPts[i].type, 
	    CritPts[i].eval[0], CritPts[i].eval[1], CritPts[i].eval[2], 
	    CritPts[i].evect[0].xd, CritPts[i].evect[0].yd, CritPts[i].evect[0].zd, 
	    CritPts[i].evect[1].xd, CritPts[i].evect[1].yd, CritPts[i].evect[1].zd, 
	    CritPts[i].evect[2].xd, CritPts[i].evect[2].yd, CritPts[i].evect[2].zd);
    
    switch(CritPts[i].type) {
    case CPT_SADDLE:
      sd++;
      break;
    case CPT_ATTRACTING_NODE:
      at++;
      break;
    case CPT_REPELLING_NODE:
      re++;
      break;
    }
  }
  fclose(fout);

  if(numCritPts > 0) {
    delete [] CritPts;
  }
  
#ifdef _DEBUG
  printf("\
Critical points: %d saddles, %d attracting nodes, %d repelling nodes.\n", 
	 sd, at, re);
  PrintElapsedTime("Phase5: Writing output file");
#endif
  printf("End \n");
  
  return 0;
}
