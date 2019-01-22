///////////////////////////////////////////////////////////////////////////////////////////////
// driver program for the potVect module
///////////////////////////////////////////////////////////////////////////////////////////////

// #include "common.h"
#include "potVect.h"

int main(int argc, char *argv[])
{
  FILE *fin;
  int L,M,N;         // Sizes in x,y,z dimensions
  int i,j,k, s;
  unsigned char *f;
  Vector *force;
  long idx, slsz, sz;
  //int measureTime = 0;
  int distCharges, fieldStrenght;
  bool inOut = false;

#ifdef _DEBUG
	SetStartTime();
#endif

  if (argc < 7) {
    printf("\n\
Usage: \n\t%s <volfile> <xs> <ys> <zs> <fieldStrenght> <pfFile> [<inOut>]\n\
  <inOut> flag specifies whether you want the field computed in the exterior\n\
          of the object too (set to 1), not just in the interior (set to 0).\n\
          Defalut value is 0 (inside only).\n\n",
	   argv[0]);
    exit(1);
  }

  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);
  fieldStrenght = atoi(argv[5]);

  if(argc > 7) {
    if(atoi(argv[7]) == 0) {
      inOut = false;
    }
    else {
      inOut = true;
    }
  }

  if ((fin = fopen(argv[1],"r")) == NULL) {
    printf("Cannot open input file %s\n",argv[1]);
    exit(1);
  }

 

#ifdef _DEBUG
	PrintElapsedTime("Phase 1: opening input file and other initializations.");
#endif

  f = new unsigned char[L*M*N];			// flags - interior, boundary, surface, ...

  slsz = L*M;		// slice size
  sz = slsz*N;

#ifdef _DEBUG
	PrintElapsedTime("Phase 2: allocating memory.");
#endif


  if ( fread(f, sizeof(unsigned char), sz, fin) < sz) {
    printf("File size is not the same as volume size\n");
    exit(1);
  }

  fclose(fin);

#ifdef _DEBUG
	PrintElapsedTime("Phase 3: reading volume from file.");
#endif

// initialize the f array
  for (idx = 0; idx < sz; idx++) {
	if (f[idx] > 0) {
	   f[idx] = INTERIOR;
	}
	else {
	   f[idx] = 0;
	}
  }



#ifdef _DEBUG
	PrintElapsedTime("Phase 5: initialize the flags.");
#endif

force = new Vector[L*M*N];			// potential field

// Compute force vectors

CalculatePotentialField(L, M, N, f, fieldStrenght, force, inOut);



// print the force vectors.
 SaveVectorField(force, L, M, N, argv[6]);
 
#ifdef _DEBUG
	PrintElapsedTime("Phase 6: saving potential field to file.");
#endif



  #ifdef _DEBUG
	PrintElapsedTime("");
  #endif

// analyse the vector field
// find the max and min vector lenghts:
double maxvl, minvl, vl;
int nrmax, nrmin;

maxvl = -0.00;
minvl = 9999.99;
nrmax = 0;
nrmin = 0;
for(i=0; i < L*M*N; i++) {
	if(f[i] == INTERIOR) {
		vl = 	sqrt((force[i].xd * force[i].xd) + (force[i].yd * force[i].yd) + (force[i].zd * force[i].zd));

		if(vl > maxvl) {
			maxvl = vl;
			nrmax = 1;
		}
		else {
			if(vl == maxvl) {
				nrmax = nrmax + 1;
			}
		}

		if(vl < minvl) {
			minvl = vl;
			nrmin = 1;
		}
		else {
			if(vl == minvl) {
				nrmin = nrmin + 1;
			}
		}
	}
}

/*
printf("Vector field analysis:\n");
printf("\tMaximum vector length: %f. Number of points with maximum value: %d\n", maxvl, nrmax);
printf("\tMinimum vector length: %f. Number of points with minimum value: %d\n", minvl, nrmin);
printf("--\n");

double perc[7] = {0.00001, 0.00005, 0.00010, 0.0005, 0.0010, 0.0020, 0.0030};
int percnr[7] = {0, 0, 0, 0, 0, 0, 0};
VoxelPosition pts[7][20];

for (k = 0; k < N; k++) {
	for (j = 0; j < M; j++) {
		for (i = 0; i < L; i++) {
			idx = k*slsz+j*L+i;

			if(f[idx] == INTERIOR) {
				vl = 	sqrt((force[idx].xd * force[idx].xd) + (force[idx].yd * force[idx].yd) + (force[idx].zd * force[idx].zd));

				for(s=0; s < 7; s++) {
					if(vl <= (minvl + ((maxvl - minvl) * perc[s] / 100.00))) {
						if(percnr[s] < 20) {
							pts[s][percnr[s]].x = i;
							pts[s][percnr[s]].y = j;
							pts[s][percnr[s]].z = k;
						}
						percnr[s] = percnr[s] + 1;
					}
				}
			}
		}
	}
}

for(j=0; j < 7; j++) {
	printf("\tNumber of points with vector value < %f %% (%f): %d\n", perc[j],
		minvl + ((maxvl - minvl) * perc[j] / 100.00), percnr[j]);
	for(k=0; k < MIN(percnr[j], 20); k++) {
		printf("\t\t(%d, %d, %d)\n", pts[j][k].x, pts[j][k].y, pts[j][k].z);
	}
}
printf("-------------------------------------\n");
*/


// free the allocated memory
delete [] f;
delete [] force;

printf("Done.\n");
 return 0; 
}



