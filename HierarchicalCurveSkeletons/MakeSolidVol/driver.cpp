#include <stdio.h>
#include <stdlib.h>

#include "makeSolidVol.h"

int main(int argc, char *argv[])
{
  FILE *fin, *fout;
  unsigned char *vol;
  int L,M,N;     		// Sizes in x,y,z dimensions
  long idx, slsz;
  int i, j, k;
  
  if (argc < 6) {
    printf("Makes a 3D object solid, filling all the holes in the object.\n");
    printf("\nUsage: \n\t%s <VolFile> <xs> <ys> <zs> <newVolFile>.\n\n",argv[0]);
    printf("	<VolFile> - volume data file name\n");
    printf("	<xs>, <ys>, <zs> - volume data size\n");
    printf("	<newVolFile> - output file name\n");
    exit(1);
  }

  if ((fin = fopen(argv[1],"r")) == NULL) {
    printf("Cannot open input file %s\n",argv[1]);
    exit(1);
  }

  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);
  
  printf("Phase 1: Allocating memory for the volume (%dx%dx%d)...\t\t", L, M, N);
  vol = new unsigned char[L*M*N];		// volume

  if((vol == NULL)) {
	  printf("Error allocating memory for the volume. Abort.\n");
	  exit(1);
  }

  printf("done.\n");

  printf("Phase 2: Reading the volume ...\t\t");  
  if (fread(vol, sizeof(unsigned char), L*M*N, fin) < L*M*N) {
  	printf("Error reading input file. Abort.\n");
   	exit(1);
  }	
  fclose(fin);

printf("done.\n");

  printf("Phase 3: Processing volume ...\t\t");  
  MakeSolidVolume(L, M, N, vol, 0, 255);
printf("done.\n");

  // output new volume to file
  printf("Phase 4: Writing the volume to output file...\t\t");
  if ((fout = fopen(argv[5],"w")) == NULL) {
    printf("Cannot open output file %s for writing\n",argv[5]);
    exit(1);
  }

  if(fwrite(vol, sizeof(unsigned char), L*M*N, fout) < L*M*N) {
  	printf("Error writing to output file. Abort.\n");
  	exit(1);
  }
  fclose(fout);

  // free the allocated memory
  delete [] vol;
  printf("done.\n");
  printf("ALL OK.\n");
  return 0;
}

