// Extract the boundary points from a input dataset: volume or polygonal mesh
// --- Input: 3D data set or polygonal mesh
// --- Output: boundary points with neighbourhood information structure
// 
// Nicu D. Cornea - Mon Jun  2 12:26:56 EDT 2003


#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <math.h>

#include "../common.h"


struct  Vector3D
{
	float x;
	float y;
	float z;
};

void PartialDerivative(float *Is, float *Isd, int direc, int L, int M, int N);
void PartialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N);
void ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v);
void RotMatrixFromAngle(float RMatrix[3][3], float phi, float theta, float psi);
void Transpose(float MatTransp[3][3], float Mat[3][3]);
void Matrix3Multiply(float Mat[3][3], float Mat1[3][3], float Mat2[3][3]);

main (int argc, char *argv[])
{
  ifstream fin;
  FILE *fout, *finv;
  char *infilename;
  Vector3D vecin, gradient0;
  float *Iu, *Iv, *Iw;
  float *Iuu, *Iuv, *Iuw, *Ivv, *Ivw, *Iww;
  float *curv;
  int *f;
  long idx, iidx, slsz, sz;
  int i, j, k, c;
  int ii, jj, kk;
  int cc;
  int L, M, N;
  float k1, k2, gLength;
  float RotMatrix[3][3];
  float RotMatrixTransp[3][3];
  float Hessian[3][3];
  float HessianPrime[3][3];
  float HessianPrime1[3][3];
	
  int Lm1, Mm1, Nm1;
  int numBound = 0;
  bool flagSurf;
  unsigned char *cf;

  SetStartTime();

  if (argc < 7)
  {
    printf("Usage: %s <volume> <Input vectors> <xs> <ys> <zs> <Output scalars>.\n", argv[0]);
    exit(1);
  }

  infilename = new char[80];
  infilename = argv[2];
  fin.open(infilename);
  if (!fin)  {
     cerr << "couldn't open " << infilename << " for input" << endl;
     return -1;
  }

  L = atoi(argv[3]);
  M = atoi(argv[4]);
  N = atoi(argv[5]);

  slsz = L*M;        // slice size
  sz = slsz*N;

  if ((fout = fopen(argv[6],"w")) == NULL)
  {
    printf("Cannot open %s for writing\n",argv[5]);
    exit(1);
  }

///////////////////
// get the boundary voxels
  if ((finv = fopen(argv[1],"r")) == NULL)
  {
    printf("Cannot open %s for reading\n",argv[1]);
    exit(1);
  }

  cf = new unsigned char[L*M*N];		// volume
  f = new int[L*M*N];

  if ( fread(cf, sizeof(unsigned char), sz, finv) < sz) {
    printf("File size is not the same as volume size\n");
    exit(1);
  }

 
#ifdef _DEBUG
	PrintElapsedTime("Phase 3: reading volume from file.");
#endif


  // initialize the f array
  for (idx = 0; idx < sz; idx++) {
	if (cf[idx] > 0) {
	   f[idx] = INTERIOR;
	}
	else {
		f[idx] = 0;
	}
  }

  // no longer need the volume 
  delete [] cf;
  fclose(finv);

#ifdef _DEBUG
	PrintElapsedTime("Phase 4: initialize the f array.");
#endif



  Lm1 = L - 1;
  Mm1 = M - 1;
  Nm1 = N - 1;

  // save all the boundary voxels in array Bound[]
	for (k = 1; k < Nm1; k++) {
		for (j = 1; j < Mm1; j++) {
			for (i = 1; i < Lm1; i++) {
				flagSurf = false;
				idx = k*slsz + j*L + i;

				// CASE 1: treat the inner layer
				if (f[idx] == 0) continue;

				//consider six face neighbors, if anyone is zero, it is a boundary voxel
				iidx = k*slsz + j*L + i-1;
				if (f[iidx] == 0) {
					flagSurf = true;
				}

				if(!flagSurf) {

				iidx = k*slsz + j*L + i+1;
				if (f[iidx] == 0) {
					flagSurf = true;
				}

				if(!flagSurf) {

				iidx = k*slsz + (j-1)*L + i;
				if (f[iidx] == 0) {
					flagSurf = true;
				}
		
				if(!flagSurf) {

				iidx = k*slsz + (j+1)*L + i;
				if (f[iidx] == 0) {
					flagSurf = true;
				}

				if(!flagSurf) {

				iidx = (k-1)*slsz + j*L + i;
				if (f[iidx] == 0) {
					flagSurf = true;
				}

				if(!flagSurf) {

				iidx = (k+1)*slsz + j*L + i;
				if (f[iidx] == 0) {
					flagSurf = true;
				}
				
				}
				}
				}
				}
				}

				// restore idx to the right value
				idx = k*slsz + j*L + i;
				if (flagSurf) {
					f[idx] = SURF;
					numBound++;
				}
			}
		}
	}

	//printf("numBound = %d \n", numBound);

#ifdef _DEBUG
	PrintElapsedTime("Phase PF-1: finding the boundary voxels.");
	printf("--Found %d boundary voxels.\n", numBound);
#endif

///////////////////
PrintElapsedTime("phase 1: opening input and output files and reading in the potential field.");

  Iu = new float[L*M*N];
  Iv = new float[L*M*N];
  Iw = new float[L*M*N];
  Iuu = new float[L*M*N];
  Ivv = new float[L*M*N];
  Iww = new float[L*M*N];
  Iuv = new float[L*M*N];
  Iuw = new float[L*M*N];
  Ivw = new float[L*M*N];
  curv = new float[L*M*N];
  

PrintElapsedTime("phase 2: allocating memory.");

  for (k = 0; k < N; k++)
     for (j = 0; j < M; j++)
        for (i = 0; i < L; i++) {
	    fin >> vecin.x >> vecin.y >> vecin.z;
	    idx = k*slsz + j*L +i;
	    Iu[idx] = vecin.x;
	    Iv[idx] = vecin.y;
	    Iw[idx] = vecin.z;
	    //f[L*M*N] indicates whether a voxel is in the foreground
	    //f[] -- 0: background; 1: surface 2: under surface
	    // if(vecin.x ==0 && vecin.y==0 && vecin.z==0) f[idx]=0; else f[idx]=2;
	}


PrintElapsedTime("phase 3: reading the petential field vectors.");

  PartialDerivative1(Iu, Iuu, 1, L, M, N);
  PartialDerivative1(Iv, Ivv, 2, L, M, N);
  PartialDerivative1(Iw, Iww, 3, L, M, N);
  PartialDerivative1(Iu, Iuv, 2, L, M, N);
  PartialDerivative1(Iu, Iuw, 3, L, M, N);
  PartialDerivative1(Iv, Ivw, 3, L, M, N);

PrintElapsedTime("phase 4: calculating the partial derivatives.");

  double maxCurvature = 0;
  int DisAway = 2;
  for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) {
	     idx = k*slsz + j*L + i;
	     if (f[idx] !=SURF) continue;
	     gradient0.x = Iu[idx];
	     gradient0.y = Iv[idx];
	     gradient0.z = Iw[idx];
	     Hessian[0][0] = Iuu[idx];
	     Hessian[0][1] = Iuv[idx];
	     Hessian[0][2] = Iuw[idx];
	     Hessian[1][0] = Iuv[idx];
	     Hessian[1][1] = Ivv[idx];
	     Hessian[1][2] = Ivw[idx];
	     Hessian[2][0] = Iuw[idx];
	     Hessian[2][1] = Ivw[idx];
	     Hessian[2][2] = Iww[idx];
	     ComputeRotMatrix(RotMatrix, gradient0);
	     Transpose(RotMatrixTransp, RotMatrix);
	     Matrix3Multiply(HessianPrime1, RotMatrixTransp, Hessian);
	     Matrix3Multiply(HessianPrime, HessianPrime1, RotMatrix);
	     k1 = 0.5*(HessianPrime[1][1]+HessianPrime[2][2])
	         +0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2)+4*Hessian[1][2]*Hessian[2][1]);
	     k2 = 0.5*(HessianPrime[1][1]+HessianPrime[2][2])
	         -0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2)+4*Hessian[1][2]*Hessian[2][1]);
	     gLength = sqrt(gradient0.x*gradient0.x +gradient0.y*gradient0.y +gradient0.z*gradient0.z);
	     //if (k2>k1) k1=k2;
	     //if(gLength !=0) curv[idx] = -(k1)/gLength;   else curv[idx]=500;//set large to be included in skeleton
	     if(gLength !=0) curv[idx] = -(k1)/pow(gLength,1.5);   else curv[idx]=10000;//set large to be included in skeleton

	     if (curv[idx] > maxCurvature)  maxCurvature = curv[idx];
	     //case 1: get neg maximum
	     if (curv[idx]>10000) {
		curv[idx]=10000;
	     }

	     if (curv[idx]<=1) curv[idx]=0;// first phase: threshold the curvature value (larger->fewer)
	                                    // thres = 0.5 for cube19x29x41
					    // thres = 11  for knight 52x100x52
					    // thres = 5 for monster
					    // thres = 10 for tetrahedral
					    // thres = 8 for mushroom
					    // thres = 3 for 767
					    // thres = 10 for human
					    // thres = 18 for TightSegColon
					    // thres = 11 for knightNS
					    // thres = 16 /11/ 6 for herachical cow
					    // thres = 16  for shark
					    // thres = 10  for hammerhead
					    // thres = 14  for dogDilate
					    // thres = 14  for woodman
					    // thres =  8  for teapot
					    // thres =     for temple
					    // thres =  9  for toilet
					    // thres = 50  for cannon
						// thres = 11.19 for heart
						// thres = 10 for anchor
						// thres = 20 for fighterplane
						// thres = 11 for trout
						// thres = 14 for bird
						// thres = 14.45 for aircar
						// thres= 15.44 for trout (100x52x52)
						// thres= 643.9 for thron

//	     fprintf(fout,"%d %d %d %f %f\n", i, j, k, curv[idx], curv[idx]);
	}
PrintElapsedTime("phase 4: computing curvatures.");

cout << "The highest Curvature is : " << maxCurvature << endl;

  // only select the local maximum voxel of curv[idx] as skeleton

  DisAway = 5;
  for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) {
	     idx = k*slsz + j*L + i;
	     //consider the six face neighbors
	     iidx = k*slsz + j*L + i-1;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = k*slsz + j*L + i+1;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = k*slsz + (j-1)*L + i;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = k*slsz + (j+1)*L + i;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = (k-1)*slsz + j*L + i;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = (k+1)*slsz + j*L + i;
	     if (curv[iidx]> curv[idx]) continue;

	     if (curv[idx] != 0)
	       //  fprintf(fout,"%d %d %d %f %f\n", i, j, k, curv[idx], curv[idx]);
		  fprintf(fout,"%d %d %d\n", i, j, k);
  }

PrintElapsedTime("phase5: getting the highest curvatures.");

  fclose(fout);
  PrintElapsedTime("");
}


void PartialDerivative(float *Is, float *Isd, int direc, int L, int M, int N)  {
  // direc = 1: x-direction   2: y-direction   3: z-direction
  long idx;
  int i,k,j;
  long slsz = L*M;
  float lessXsum, moreXsum, lessYsum, moreYsum, lessZsum, moreZsum;
  int disaway = 1;
  for (k = disaway; k < N-disaway; k++)
     for (j = disaway; j < M-disaway; j++)
        for (i = disaway; i < L-disaway; i++) {
	    idx = k*slsz + j*L +i;
	    if (direc == 1) {
	         lessXsum = Is[(k-1)*slsz +(j-1)*L +(i-1)] + Is[(k-1)*slsz +j*L +(i-1)] +Is[(k-1)*slsz +(j+1)*L +(i-1)]
		           +Is[  k*slsz + (j-1)*L +(i-1)] + Is[   k*slsz + j*L +(i-1)] + Is[  k*slsz + (j+1)*L +(i-1)]
			   +Is[(k+1)*slsz +(j-1)*L +(i-1)] + Is[(k+1)*slsz +j*L +(i-1)] +Is[(k+1)*slsz +(j+1)*L +(i-1)];
		 moreXsum = Is[(k-1)*slsz +(j-1)*L +(i+1)] + Is[(k-1)*slsz +j*L +(i+1)] +Is[(k-1)*slsz +(j+1)*L +(i+1)]
		           +Is[  k*slsz + (j-1)*L +(i+1)] + Is[   k*slsz + j*L +(i+1)] + Is[  k*slsz + (j+1)*L +(i+1)]
			   +Is[(k+1)*slsz +(j-1)*L +(i+1)] + Is[(k+1)*slsz +j*L +(i+1)] +Is[(k+1)*slsz +(j+1)*L +(i+1)];
		 Isd[idx] = moreXsum - lessXsum;
	    }
	    else if(direc == 2) {
	         lessYsum = Is[(k-1)*slsz +(j-1)*L +(i-1)] + Is[(k-1)*slsz +(j-1)*L + i] +Is[(k-1)*slsz +(j-1)*L +(i+1)]
		           +Is[  k*slsz + (j-1)*L +(i-1)] + Is[   k*slsz + (j-1)*L +i] + Is[  k*slsz + (j-1)*L +(i+1)]
			   +Is[(k+1)*slsz +(j-1)*L +(i-1)] + Is[(k+1)*slsz +(j-1)*L +i] +Is[(k+1)*slsz +(j-1)*L +(i+1)];
	         moreYsum = Is[(k-1)*slsz +(j+1)*L +(i-1)] + Is[(k-1)*slsz +(j+1)*L + i] +Is[(k-1)*slsz +(j+1)*L +(i+1)]
		           +Is[  k*slsz + (j+1)*L +(i-1)] + Is[   k*slsz + (j+1)*L +i] + Is[  k*slsz + (j+1)*L +(i+1)]
			   +Is[(k+1)*slsz +(j+1)*L +(i-1)] + Is[(k+1)*slsz +(j+1)*L +i] +Is[(k+1)*slsz +(j+1)*L +(i+1)];
		 Isd[idx] = moreYsum - lessYsum;
	    }
	    else {
	         lessZsum = Is[(k-1)*slsz +(j-1)*L +(i-1)] + Is[(k-1)*slsz +(j-1)*L + i] +Is[(k-1)*slsz +(j-1)*L +(i+1)]
		           +Is[ (k-1)*slsz + j *L + (i-1)] + Is[(k-1)*slsz + j*L + i ] + Is[(k-1)*slsz + j*L + (i+1)]
			   +Is[(k-1)*slsz +(j+1)*L +(i-1)] + Is[(k-1)*slsz +(j+1)*L +i] +Is[(k-1)*slsz +(j+1)*L +(i+1)];
	         moreZsum = Is[(k+1)*slsz +(j-1)*L +(i-1)] + Is[(k+1)*slsz +(j-1)*L + i] +Is[(k+1)*slsz +(j-1)*L +(i+1)]
		           +Is[ (k+1)*slsz + j *L + (i-1)] + Is[(k+1)*slsz + j*L + i ] + Is[(k+1)*slsz + j*L + (i+1)]
			   +Is[(k+1)*slsz +(j+1)*L +(i-1)] + Is[(k+1)*slsz +(j+1)*L +i] +Is[(k+1)*slsz +(j+1)*L +(i+1)];
		 Isd[idx] = moreZsum - lessZsum;
	    }
        }
}


void PartialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N)  {
  // direc = 1: x-direction   2: y-direction   3: z-direction
  long idx;
  int i,k,j;
  long slsz = L*M;
  float lessXsum, moreXsum, lessYsum, moreYsum, lessZsum, moreZsum;
  int disaway = 1;
  for (k = disaway; k < N-disaway; k++)
     for (j = disaway; j < M-disaway; j++)
        for (i = disaway; i < L-disaway; i++) {
	    idx = k*slsz + j*L +i;
	    if (direc == 1) {
	         lessXsum = Is[k*slsz + j*L +(i-1)];
		 moreXsum = Is[k*slsz + j*L +(i+1)];
		 Isd[idx] = moreXsum - lessXsum;
	    }
	    else if(direc == 2) {
	         lessYsum = Is[k*slsz + (j-1)*L +i];
	         moreYsum = Is[k*slsz + (j+1)*L +i];
		 Isd[idx] = moreYsum - lessYsum;
	    }
	    else {
	         lessZsum = Is[(k-1)*slsz + j*L + i ];
	         moreZsum = Is[(k+1)*slsz + j*L + i ];
		 Isd[idx] = moreZsum - lessZsum;
	    }
        }
}




void ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v) {
   // Find the matrix that can rotate any vector v to x-axis direction
   float phi, theta;
   Vector3D vector1;
   float RotateMatrix1[3][3];

   theta = atan2(v.z, v.y);
   RotMatrixFromAngle(RotateMatrix1, 0, -theta, 0);
   vector1.x = RotateMatrix1[0][0]*v.x +RotateMatrix1[0][1]*v.y +RotateMatrix1[0][2]*v.z;
   vector1.y = RotateMatrix1[1][0]*v.x +RotateMatrix1[1][1]*v.y +RotateMatrix1[1][2]*v.z;
   vector1.z = RotateMatrix1[2][0]*v.x +RotateMatrix1[2][1]*v.y +RotateMatrix1[2][2]*v.z;
   phi = atan2(vector1.y, vector1.x);
   RotMatrixFromAngle(RotateMatrix, -phi, -theta, 0);
}


void RotMatrixFromAngle(float RMatrix[3][3], float phi, float theta, float psi) {
  // rotation matrix is      R(0,0) R(0,1) R(0,2)
  //	        	     R(1,0) R(1,1) R(1,2)
  //   			     R(2,0) R(2,1) R(2,2)
   RMatrix[0][0] = cos(phi)*cos(psi) - cos(theta)*sin(phi)*sin(psi);
   RMatrix[1][0] = cos(psi)*sin(phi) + cos(phi)*cos(theta)*sin(psi);
   RMatrix[2][0] = sin(psi)*sin(theta);
   RMatrix[0][1] = -(cos(psi)*cos(theta)*sin(phi)) - cos(phi)*sin(psi);
   RMatrix[1][1] = cos(phi)*cos(psi)*cos(theta) - sin(phi)*sin(psi);
   RMatrix[2][1] = cos(psi)*sin(theta);
   RMatrix[0][2] = sin(phi)*sin(theta);
   RMatrix[1][2] = -(cos(phi)*sin(theta));
   RMatrix[2][2] = cos(theta);
}


void Transpose(float MatTransp[3][3], float Mat[3][3])  {
   for(int j=0; j<=2; j++)
        for(int i=0; i<=2; i++)  {
             MatTransp[j][i] = Mat[i][j];
	}
}


void Matrix3Multiply(float Mat[3][3], float Mat1[3][3], float Mat2[3][3]) {
   for(int j=0; j<=2; j++)
        for(int i=0; i<=2; i++)  {
             Mat[j][i] = Mat1[j][0]*Mat2[0][i] +Mat1[j][1]*Mat2[1][i] +Mat1[j][2]*Mat2[2][i];
	}
}

