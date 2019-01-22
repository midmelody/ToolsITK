////////////////////////////////////////////////////////////////////////////
// driver program for the expandVol module
////////////////////////////////////////////////////////////////////////////

#include "expandVol.h"

int main(int argc, char *argv[])
{
  FILE *fout;
  unsigned char *cf;
  int L,M,N;         // Sizes in x,y,z dimensions
  int i,j,k, s;
  unsigned char *f;
  long idx, slsz, sz;
  int layers;

#ifdef _DEBUG
	SetStartTime();
#endif

  if (argc < 7) {
    printf("\nUsage: \n\t%s <volfile> <xs> <ys> <zs> <nroflayers> <paddedvolfile>\n\n",argv[0]);
    exit(1);
  }

  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);
  layers = atoi(argv[5]);

#ifdef _DEBUG
	PrintElapsedTime("Phase 1: opening input file and other initializations.");
#endif

  f = new unsigned char[L*M*N];			// flags - interior, boundary, surface, ...

  slsz = L*M;		// slice size
  sz = slsz*N;

#ifdef _DEBUG
  PrintElapsedTime("Phase 2: allocating memory.");
#endif

  ReadVolume(argv[1], L, M, N, &f);


#ifdef _DEBUG
  PrintElapsedTime("Phase 3: reading volume from file.");
#endif

  // flag the volume: INTERIOR, EXTERIOR, SURF
  FlagVolume(f, L, M, N);


#ifdef _DEBUG
  PrintElapsedTime("Phase 4: flag the volume.");
#endif

  printf("Padding the object with %d layers of extra voxels\n", layers);
  ExpandNLayers(L, M, N, f, layers);


#ifdef _DEBUG
  PrintElapsedTime("Phase 5: expand the object.");
#endif

  // make the volume binary
  for(idx=0; idx < L*M*N; idx++) {
    if(f[idx] != EXTERIOR) {
      f[idx] = 255;
    }
  }

// save the new volume
  if ((fout = fopen(argv[6],"w")) == NULL) {
    printf("Cannot open output file %s for writing\n",argv[6]);
    exit(1);
  }

  if(fwrite(f, sizeof(unsigned char), sz, fout) != sz) {
    printf("Error writing to the output file\n");
    exit(1);
}

#ifdef _DEBUG
  PrintElapsedTime("Phase 6: saving new volume to file.");
#endif
  fclose(fout);


#ifdef _DEBUG
  PrintElapsedTime("");
#endif

// free the allocated memory
  delete [] f;

  printf("Done.\n");
  return 0;
}



