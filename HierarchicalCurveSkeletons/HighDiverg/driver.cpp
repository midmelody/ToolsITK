//
// driver program for the HighDiverg module
//

#include "../common.h"
#include "HighDiverg.h"

#include <fstream>
#include <iostream>

main(int argc, char *argv[])
{
  FILE *fout, *fvec;
  int L,M,N;         // Sizes in x,y,z dimensions

  int i,j,k;
  Vector *fv;

  VoxelPositionDouble *HDPts;
  int numHDPts;

  long idx, slsz;

  unsigned char *cf, *f;
  int perc;
  bool inOut = false;

#ifdef _DEBUG
  SetStartTime();
#endif

  if (argc < 8)
  {
    printf("\n\
Usage: %s <volume file> <xs> <ys> <zs> <vector file> <perc> \n\
          <high divergence points file> [inOut].\n\
<inOut> flag specifies whether you want to look for high divergence points \n\
        outside the object too (set to 1), \n\
        not just in the interior (set to 0).\n\
        Defalut value is 0 (inside only).\n",
	   argv[0]);
    exit(1);
  }
  

  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);
  perc = atoi(argv[6]);

  if(argc > 8) {
    if(atoi(argv[8]) == 0) {
      inOut = false;
    }
    else {
      inOut = true;
    }
  }

  
  slsz = L*M;		// slice size

  // read volume
  ReadVolume(argv[1], L, M, N, &cf);

  // set flags
  SetFlags(cf, L, M, N, &f);
  delete [] cf;

#ifdef _DEBUG
  PrintElapsedTime("Phase1: reading volume and setting flags");
#endif


 if((fv = new Vector[L*M*N]) == NULL) {
	printf("UPS! - Error allocating memory for the input vector field. Abort.\n");
	exit(1);
  }

 ReadVectorField(fv, L, M, N, argv[5]);

  for(idx = 0; idx < L*M*N; idx++) {
    if(fv[idx].xd !=0) {
    	printf("NOT zero idx = %d!!!\n", idx);
    	break;
    }
  }
  
#ifdef _DEBUG
  PrintElapsedTime("Phase2: reading in force vectors");
#endif


	HDPts = NULL;
	numHDPts = 0;

	GetHighDivergencePoints(fv, L, M, N, f, perc, &HDPts, &numHDPts, inOut);

#ifdef _DEBUG
	  PrintElapsedTime("Phase4: Getting the high divergence points.");
	  printf("**  %d high divergence points returned.\n", numHDPts);
#endif

	delete [] f;
	delete [] fv;

	if ((fout = fopen(argv[7],"w")) == NULL) {
		printf("Cannot open %s for writing\n",argv[7]);
		exit(1);
	}

	for(i=0; i < numHDPts; i++) {
		// printf("writing point (%f, %f, %f) to output file.\n", 
		//	HDPts[i].x, HDPts[i].y, HDPts[i].z);
		fprintf(fout,"%f %f %f 1\n", 
			HDPts[i].x, HDPts[i].y, HDPts[i].z);
	}
	fclose(fout);

	delete [] HDPts;

#ifdef _DEBUG
   	  PrintElapsedTime("Phase5: Writing output file");
#endif
	printf("End \n");

}
