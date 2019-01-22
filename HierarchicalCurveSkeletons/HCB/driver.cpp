//
//
// driver program for the high curvature detection module
//
//

// #include "../common.h"
#include "hcBound.h"


void main(int argc, char *argv[]) {
	FILE* fin, *finv, *fout;
	int L, M, N;
	unsigned char *vol;
	VoxelPosition	*HighCurvaturePoints;
	int numHC;
	int i, j, k, kk;

#ifdef _DEBUG
	SetStartTime();
#endif

	if (argc < 6) {
    		printf("Usage: %s <volume> <xs> <ys> <zs> <Output scalars>.\n", argv[0]);
   	 	exit(1);
  	}

  	if ((finv = fopen(argv[1],"r")) == NULL) {
		printf("Cannot open %s for reading\n",argv[5]);
		exit(1);
	}

  	L = atoi(argv[2]);
  	M = atoi(argv[3]);
  	N = atoi(argv[4]);

	vol = new unsigned char[L*M*N];
	if (fread(vol, sizeof(unsigned char), L*M*N, finv) < L*M*N) {
    		printf("File size is not the same as volume size\n");
    		exit(1);
  	}

#ifdef _DEBUG
	PrintElapsedTime("Phase 1: reading volume from file.");
#endif

	HighCurvaturePoints = NULL;
	numHC = 0;
	GetHighCurvatureBoundaryPoints(vol, L, M, N, &HighCurvaturePoints, &numHC,1.00);


	// OUTPUT
	if ((fout = fopen(argv[5],"w")) == NULL) {
		printf("Cannot open %s for writing\n",argv[6]);
		exit(1);
	}

/*
	// write all boundary points to output file
	for(i=0; i < numBound; i++) {
		fprintf(fout, "%d %d %d\n", Bound[i].point.x, Bound[i].point.y, Bound[i].point.z);
	}
*/
/*
	// write normals to output file

	for (k = 1; k < N-1; k++) {
		for (j = 1; j < M-1; j++) {
			for (i = 1; i < L-1; i++) {
				found = false;
				for(kk=0; kk < numBound; kk++) {
					if(	(Bound[kk].point.x == i) &&
						(Bound[kk].point.y == j) &&
						(Bound[kk].point.z == k))
					{
						fprintf(fout, "%f %f %f\n", Bound[kk].normal.xd, Bound[kk].normal.yd, Bound[kk].normal.zd);
						found = true;
					}
				}
				if(!found) {
					fprintf(fout, "0.00 0.00 0.00\n");
				}
			}
		}
	}
*/

	// write the high curvature points to the output file.
	for(i=0; i < numHC; i++) {
		fprintf(fout, "%d  %d  %d\n", HighCurvaturePoints[i].x, HighCurvaturePoints[i].y, HighCurvaturePoints[i].z);
	}


	fclose(fout);

#ifdef _DEBUG
	PrintElapsedTime("Phase 3: writing output file.");
#endif

	printf("Done\n");
}

