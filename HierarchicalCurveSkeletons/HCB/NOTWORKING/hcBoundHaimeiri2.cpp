//
// This version uses the Householder transformation and the Givens rotation to find the eigenvectors
//

#include "BoundPoints.h"

#define TRACE

void main(int argc, char *argv[]) {
	FILE* fin, *finv, *fout;
	int L, M, N;
	unsigned char *vol;
	bdElement* Bound;
	int numBound;
	int i, j, k, kk;
	bool found;

	Vector Ti, pmq, pd1, pd2, Wv;
	double sl, ki, length, wi, c, ca, cb, cc, t, delta, e1, e2, A, B, C, sinasq, cosa, k1, k2, maxk1, mink1, l1, l2;
	double Ma[3][3], Qv[3][3];

#ifdef _DEBUG
	SetStartTime();
#endif


	if (argc < 7) {
    		printf("Usage: %s <volume> <Input vectors> <xs> <ys> <zs> <Output scalars>.\n", argv[0]);
   	 	exit(1);
  	}

  	if ((finv = fopen(argv[1],"r")) == NULL) {
		printf("Cannot open %s for reading\n",argv[5]);
		exit(1);
	}

  	L = atoi(argv[3]);
  	M = atoi(argv[4]);
  	N = atoi(argv[5]);

	vol = new unsigned char[L*M*N];
	if (fread(vol, sizeof(unsigned char), L*M*N, finv) < L*M*N) {
    		printf("File size is not the same as volume size\n");
    		exit(1);
  	}

#ifdef _DEBUG
	PrintElapsedTime("Phase 1: reading volume from file.");
#endif

	Bound = NULL;
	GetBoundaryPoints(vol, L, M, N,	&Bound, &numBound);

#ifdef _DEBUG
	PrintElapsedTime("");
#endif

// get the maximum principal direction at each point on the boundary, using the modified Taubin method,
//	as described in
//		Eyal Haimeri, Ilan Shimshoni
//		"Estimating teh Principal Curvatures and the Darboux Frame from Real 3D Range Data", 2002

	maxk1 = -999999999;
	mink1 = 999999999;
	for(i=0; i < numBound; i++) {
		Bound[i].curvature = 0.0;
#ifdef TRACE
		printf("i = %d (%d, %d, %d). neighbors: %d.\n", i, Bound[i].point.x, Bound[i].point.y, Bound[i].point.z, Bound[i].nrOfNeighbors);
		printf("\tNormal: (%f, %f, %f)\n", Bound[i].normal.xd, Bound[i].normal.yd, Bound[i].normal.zd);
#endif
		// calculate the sum of all distances from this point to it's neighbors
		//	- used for the weights.
		sl = 0;
		for(j=0; j < Bound[i].nrOfNeighbors; j++) {
			sl = sl +  (double) 1.0 / sqrt(
			(Bound[Bound[i].neighbors[j]].point.x - Bound[i].point.x) * (Bound[Bound[i].neighbors[j]].point.x - Bound[i].point.x) +
			(Bound[Bound[i].neighbors[j]].point.y - Bound[i].point.y) * (Bound[Bound[i].neighbors[j]].point.y - Bound[i].point.y) +
			(Bound[Bound[i].neighbors[j]].point.z - Bound[i].point.z) * (Bound[Bound[i].neighbors[j]].point.z - Bound[i].point.z)
			);
		}
#ifdef TRACE
		printf("\ttotal distances: %f\n", sl);
#endif

		Ma[0][0] = 0;	Ma[0][1] = 0;	Ma[0][2] = 0;
		Ma[1][0] = 0;	Ma[1][1] = 0;	Ma[1][2] = 0;
		Ma[2][0] = 0;	Ma[2][1] = 0;	Ma[2][2] = 0;

		for(j=0; j < Bound[i].nrOfNeighbors; j++) {
			// estimating the ki term
#ifdef TRACE
			printf("\tNeighbor%d: (%d, %d, %d)\n", j, Bound[Bound[i].neighbors[j]].point.x, Bound[Bound[i].neighbors[j]].point.y, Bound[Bound[i].neighbors[j]].point.z);
#endif
			pmq.xd = Bound[Bound[i].neighbors[j]].point.x - Bound[i].point.x;
			pmq.yd = Bound[Bound[i].neighbors[j]].point.y - Bound[i].point.y;
			pmq.zd = Bound[Bound[i].neighbors[j]].point.z - Bound[i].point.z;
#ifdef TRACE
			printf("\t\tpmq: (%f, %f, %f)\n", pmq.xd, pmq.yd, pmq.zd);
#endif

			ki = (double) 2.0 * 
				 (double) ((Bound[i].normal.xd * (pmq.xd)) + (Bound[i].normal.yd * (pmq.yd)) + (Bound[i].normal.zd * (pmq.zd))) /
				 (double) (pmq.xd * pmq.xd + pmq.yd * pmq.yd + pmq.zd * pmq.zd);
#ifdef TRACE
			printf("\t\tki = %f\n", ki);
#endif

			// calculation the normalized projection of pmq onto the tangent plane = Ti
			// project the neighbor on to the tangent plane
			t = -(Bound[i].normal.xd * pmq.xd + Bound[i].normal.yd * pmq.yd + Bound[i].normal.zd * pmq.zd);
#ifdef TRACE
			printf("\t\tt = %f\n", t);
#endif
			Ti.xd = t*Bound[i].normal.xd + pmq.xd;
			Ti.yd = t*Bound[i].normal.yd + pmq.yd;
			Ti.zd = t*Bound[i].normal.zd + pmq.zd;

			// normalize Ti
			length = sqrt(Ti.xd * Ti.xd + Ti.yd * Ti.yd + Ti.zd * Ti.zd);
			Ti.xd = Ti.xd / length;
			Ti.yd = Ti.yd / length;
			Ti.zd = Ti.zd / length;
#ifdef TRACE
			printf("\t\tprojection of pmq on the tangent plane: (%f, %f, %f)\n", Ti.xd, Ti.yd, Ti.zd);
#endif

			// the weights wi
			wi = ((double) 1.0 / sqrt(pmq.xd * pmq.xd + pmq.yd * pmq.yd + pmq.zd * pmq.zd)) / sl;
#ifdef TRACE
			printf("\t\tw%d = %f\n", j, wi);
#endif

			// the matrix M
			c = wi * ki;

			Ma[0][0] = Ma[0][0] + (c * Ti.xd * Ti.xd);
			Ma[0][1] = Ma[0][1] + (c * Ti.xd * Ti.yd);
			Ma[0][2] = Ma[0][2] + (c * Ti.xd * Ti.zd);
			Ma[1][0] = Ma[1][0] + (c * Ti.yd * Ti.xd);
			Ma[1][1] = Ma[1][1] + (c * Ti.yd * Ti.yd);
			Ma[1][2] = Ma[1][2] + (c * Ti.yd * Ti.zd);
			Ma[2][0] = Ma[2][0] + (c * Ti.zd * Ti.xd);
			Ma[2][1] = Ma[2][1] + (c * Ti.zd * Ti.yd);
			Ma[2][2] = Ma[2][2] + (c * Ti.zd * Ti.zd);
		}

#ifdef TRACE
		printf("\tThe Ma matrix:\n");
		printf("\t\t%f  %f  %f\n", Ma[0][0], Ma[0][1], Ma[0][2]);
		printf("\t\t%f  %f  %f\n", Ma[1][0], Ma[1][1], Ma[1][2]);
		printf("\t\t%f  %f  %f\n", Ma[2][0], Ma[2][1], Ma[2][2]);
#endif

		l1 = 	((1.00 - Bound[i].normal.xd) * (1.00 - Bound[i].normal.xd)) + 
				((0.00 - Bound[i].normal.yd) * (0.00 - Bound[i].normal.yd)) + 
				((0.00 - Bound[i].normal.zd) * (0.00 - Bound[i].normal.zd));

		l2 = 	((1.00 + Bound[i].normal.xd) * (1.00 + Bound[i].normal.xd)) + 
				((0.00 + Bound[i].normal.yd) * (0.00 + Bound[i].normal.yd)) + 
				((0.00 + Bound[i].normal.zd) * (0.00 + Bound[i].normal.zd));
		
		length = 0;
		if(l1 > l2) {
			Wv.xd = 1.00 - Bound[i].normal.xd;
			Wv.yd = 0.00 - Bound[i].normal.yd;
			Wv.zd = 0.00 - Bound[i].normal.zd;
			length = sqrt(l1);
		}
		else {
			Wv.xd = 1.00 + Bound[i].normal.xd;
			Wv.yd = 0.00 + Bound[i].normal.yd;
			Wv.zd = 0.00 + Bound[i].normal.zd;
			length = sqrt(l2);
		}
		if(length == 0) {
			printf("UPS! - length of vector = 0! - Abort.\n");
			exit(1);
		}
		Wv.xd = Wv.xd / length;
		Wv.yd = Wv.yd / length;
		Wv.zd = Wv.zd / length;

#ifdef TRACE
		printf("\tWv = (%f  %f  %f)\n", Wv.xd, Wv.yd, Wv.zd);
#endif

		Qv[0][0] = 1.00 - 2*Wv.xd*Wv.xd;
		Qv[0][1] = 0.00 - 2*Wv.xd*Wv.yd;
		Qv[0][2] = 0.00 - 2*Wv.xd*Wv.zd;
		Qv[1][0] = 0.00 - 2*Wv.yd*Wv.xd;
		Qv[1][1] = 1.00 - 2*Wv.yd*Wv.yd;
		Qv[1][2] = 0.00 - 2*Wv.yd*Wv.zd;
		Qv[2][0] = 0.00 - 2*Wv.zd*Wv.xd;
		Qv[2][1] = 0.00 - 2*Wv.zd*Wv.yd;
		Qv[2][2] = 1.00 - 2*Wv.zd*Wv.zd;

#ifdef TRACE
		printf("\tThe Qv matrix:\n");
		printf("\t\t%f  %f  %f\n", Qv[0][0], Qv[0][1], Qv[0][2]);
		printf("\t\t%f  %f  %f\n", Qv[1][0], Qv[1][1], Qv[1][2]);
		printf("\t\t%f  %f  %f\n", Qv[2][0], Qv[2][1], Qv[2][2]);
#endif

		
		

		ca = 	Ma[0][0] + Ma[1][1] + Ma[2][2];
		cb = 	Ma[2][0]*Ma[0][2] + Ma[1][0]*Ma[0][1] + Ma[2][1]*Ma[1][2] -
				Ma[0][0]*Ma[1][1] - Ma[0][0]*Ma[2][2] - Ma[1][1]*Ma[2][2];
		cc = 	Ma[0][1]*Ma[1][2]*Ma[2][0] + Ma[0][2]*Ma[1][0]*Ma[2][1] + Ma[0][0]*Ma[1][1]*Ma[2][2] -
				Ma[2][0]*Ma[0][2]*Ma[1][1] - Ma[1][0]*Ma[0][1]*Ma[2][2] - Ma[0][0]*Ma[2][1]*Ma[1][2];

#ifdef TRACE
		printf("\tca = %f; cb = %f; cc = %f\n", ca, cb, cc);
#endif

		// cc has to be 0 or else something is wrong
		if(fabs(cc) > 0.0000001) {
			printf("UPS cc is not 0 at i = %d. cc = %f\n", i, cc);
			exit(1);
		}

		delta = ca*ca + 4*cb;
		if(fabs(delta) < 0.0000001) {
			delta = 0.0;
		}
		if(delta < 0) {
			printf("UPS delta is < 0 at i = %d. delta = %f, ca = %f, cb = %f\n", i, delta, ca, cb);
			exit(1);
		}
#ifdef TRACE
		printf("\tdelta = %f\n", delta);
#endif

		e1 = (ca + sqrt(delta)) / 2.0;
		e2 = (ca - sqrt(delta)) / 2.0;
		if(e1 < e2) {
			printf("UPS e1 < e2 at i = %d. delta = %f, e1 = %f, e2 = %f\n", i, delta, e1, e2);
			exit(1);
		}
#ifdef TRACE
		printf("\te1 = %f, e2 = %f\n", e1, e2);
#endif

		// find the principal directions (pd1, pd2) = the eigenvectors of matrix M associated with eigenvalues e1 and e2
		// the maximum principal direction should be the one coresponding to the eigenvalue with the maximum absolute value

		// use the eigenvalue with the highest absolute value 
		if(fabs(e1) < fabs(e2)) {
			t = e1;
			e1 = e2;
			e2 = t;
		}

		// pd1
		pd1.xd = 1.0;
		t = Ma[0][1]*Ma[1][2] - Ma[0][2]*(Ma[1][1] - e1);
		if(t == 0.0) {
			pd1.zd = 0.0;
		}
		else {
			pd1.zd = ((Ma[1][1] - e1)*(Ma[0][0] - e1) - (Ma[1][0]*Ma[0][1])) / t;
		}
		if(Ma[0][1] == 0) {
			pd1.yd = 0.0;
		}
		else {
			pd1.yd = (e1 - Ma[0][0] - Ma[0][2]*pd1.zd) / Ma[0][1];
		}
		// normalize the pd1 vector
		length = sqrt(pd1.xd * pd1.xd + pd1.yd * pd1.yd + pd1.zd * pd1.zd);
		if(length == 0) {
			printf("UPS! - length of vector = 0! - Abort.\n");
			exit(1);
		}
		pd1.xd = pd1.xd / length;
		pd1.yd = pd1.yd / length;
		pd1.zd = pd1.zd / length;
#ifdef TRACE
		printf("\tpd1 = (%f, %f, %f)\n", pd1.xd, pd1.yd, pd1.zd);
#endif

		// find the A, B, C coeficients
		A = 0;
		B = 0;
		C = 0;
		for(j=0; j < Bound[i].nrOfNeighbors; j++) {
			pmq.xd = Bound[Bound[i].neighbors[j]].point.x - Bound[i].point.x;
			pmq.yd = Bound[Bound[i].neighbors[j]].point.y - Bound[i].point.y;
			pmq.zd = Bound[Bound[i].neighbors[j]].point.z - Bound[i].point.z;

			// the weights wi
			wi = ((double) 1.0 / sqrt(pmq.xd * pmq.xd + pmq.yd * pmq.yd + pmq.zd * pmq.zd)) / sl;

			// normalize pmq
			length = sqrt(pmq.xd * pmq.xd + pmq.yd * pmq.yd + pmq.zd * pmq.zd);
			if(length == 0) {
				printf("UPS! - length of vector = 0! - Abort.\n");
				exit(1);
			}
			pmq.xd = pmq.xd / length;
			pmq.yd = pmq.yd / length;
			pmq.zd = pmq.zd / length;

			// cosa, sinasq = the angle between the pd1 and pmq vectors
			// cos = dot product; sinsq = 1 - cos*cos;
			cosa = pmq.xd*pd1.xd + pmq.yd*pd1.yd + pmq.zd*pd1.zd;
			sinasq = 1 - (cosa*cosa);

			A = A + wi*(cosa*cosa)*(cosa*cosa);
			B = B + wi*(sinasq)*(cosa*cosa);
			C = C + wi*(sinasq)*(sinasq);
		}
#ifdef TRACE
		printf("\tA = %f, B = %f, C = %f\n", A, B, C);
#endif

		// the principal curvatures
		if((B*B - C*A) == 0) {
			// printf("UPS! - (B*B - C*A) = 0 !?? - Abort.\n");
			k1 = 0;
		}
		else {
			k1 = (B*e2 - C*e1) / (B*B - C*A);
		}
#ifdef TRACE
		printf("\tk1 = %f", k1);
#endif
		if(B == 0) {
			// printf("UPS! - B = 0 !?? - Abort.\n");
			k2 = 0;
		}
		else {
			k2 = (e1 - A*k1) / B;
		}
#ifdef TRACE
		printf("k2 = %f;\n", k2);
#endif
		k1 = fabs(k1);
		k2 = fabs(k2);

		if(k1 < k2) {
			// THIS SHOULDN'T HAPPEN !!!
			printf("UPS! |k1| < |k2| !?? - Try again.\n");

			//printf("i = %d (%d, %d, %d)\n", i, Bound[i].point.x, Bound[i].point.y, Bound[i].point.z);
			//printf("\tNormal: (%f, %f, %f)\n", Bound[i].normal.xd, Bound[i].normal.yd, Bound[i].normal.zd);
			//printf("\tdelta = %f, e1 = %f, e2 = %f\n", delta, e1, e2);
			//exit(1);

			// try e2 and pd2.
			// pd2
			pd2.xd = 1.0;
			t = Ma[0][1]*Ma[1][2] - Ma[0][2]*(Ma[1][1] - e2);
			if(t == 0.0) {
				pd2.zd = 0.0;
			}
			else {
				pd2.zd = ((Ma[1][1] - e2)*(Ma[0][0] - e2) - (Ma[1][0]*Ma[0][1])) / t;
			}
			if(Ma[0][1] == 0) {
				pd2.yd = 0.0;
			}
			else {
				pd2.yd = (e2 - Ma[0][0] - Ma[0][2]*pd2.zd) / Ma[0][1];
			}
			// normalize the pd2 vector
			length = sqrt(pd2.xd * pd2.xd + pd2.yd * pd2.yd + pd2.zd * pd2.zd);
			if(length == 0) {
				printf("UPS! - length of vector = 0! - Abort.\n");
				exit(1);
			}
			pd2.xd = pd2.xd / length;
			pd2.yd = pd2.yd / length;
			pd2.zd = pd2.zd / length;
#ifdef TRACE	
			printf("\tpd2 = (%f, %f, %f)\n", pd2.xd, pd2.yd, pd2.zd);
#endif
			// find the A, B, C coeficients
			A = 0;
			B = 0;
			C = 0;
			for(j=0; j < Bound[i].nrOfNeighbors; j++) {
				pmq.xd = Bound[Bound[i].neighbors[j]].point.x - Bound[i].point.x;
				pmq.yd = Bound[Bound[i].neighbors[j]].point.y - Bound[i].point.y;
				pmq.zd = Bound[Bound[i].neighbors[j]].point.z - Bound[i].point.z;

				// the weights wi
				wi = ((double)1.0 / sqrt(pmq.xd * pmq.xd + pmq.yd * pmq.yd + pmq.zd * pmq.zd)) / sl;

				// normalize pmq
				length = sqrt(pmq.xd * pmq.xd + pmq.yd * pmq.yd + pmq.zd * pmq.zd);
				if(length == 0) {
					printf("UPS! - length of vector = 0! - Abort.\n");
					exit(1);
				}
				pmq.xd = pmq.xd / length;
				pmq.yd = pmq.yd / length;
				pmq.zd = pmq.zd / length;

				// cosa, sina = the angled between the pd1 and pmq vectors
				// cos = dot product; sin - 1 - cos*cos;
				cosa = pmq.xd*pd2.xd + pmq.yd*pd2.yd + pmq.zd*pd2.zd;
				sinasq = 1 - cosa*cosa;

				A = A + wi*(cosa*cosa)*(cosa*cosa);
				B = B + wi*(sinasq)*(cosa*cosa);
				C = C + wi*(sinasq)*(sinasq);
			}
#ifdef TRACE
			printf("\tA = %f, B = %f, C = %f\n", A, B, C);
#endif

			// the principal curvatures
			if((B*B - C*A) == 0) {
				// printf("UPS! - (B*B - C*A) = 0 !?? - Abort.\n");
				k1 = 0;
			}
			else {
				k1 = (B*e2 - C*e1) / (B*B - C*A);
			}
#ifdef TRACE
			printf("\tk1 = %f", k1);
#endif
			if(B == 0) {
				// printf("UPS! - B = 0 !?? - Abort.\n");
				k2 = 0;
			}
			else {
				k2 = (e1 - A*k1) / B;
			}
#ifdef TRACE
			printf("k2 = %f;\n", k2);
#endif

			k1 = fabs(k1);
			k2 = fabs(k2);
			
			if(k1 < k2) {
				printf("UPS! k1 is still < k2 !?? - Abort.\n");
				exit(1);
			}
		}
		
		Bound[i].curvature = k1;
		Bound[i].curvk2 = k2;

		if(k1 > maxk1) maxk1 = k1;
		if(k1 < mink1) mink1 = k1;
	}

#ifdef _DEBUG
	PrintElapsedTime("Finding principal curvature");
	printf(" -- Max curvature = %f, min curvature = %f\n", maxk1, mink1);
#endif

	// OUTPUT
	if ((fout = fopen(argv[6],"w")) == NULL) {
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
	// threshold = ...% of the highest curvature value
	// output all the points with curvature value > threshold
	t = 0.0;
	for(i=0; i < numBound; i++) {
		if(Bound[i].curvature > t) {
			fprintf(fout, "%d  %d  %d \t k1 = %f, k2 = %f\n", Bound[i].point.x, Bound[i].point.y, Bound[i].point.z, Bound[i].curvature, Bound[i].curvk2);
		}
	}


	fclose(fout);

#ifdef _DEBUG
	PrintElapsedTime("Phase 3: writing output file.");
#endif

	printf("Done\n");
}

