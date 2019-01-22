//
// hcBoundHaimeri3
//
// get the maximum principal direction at each point on the boundary, using the modified Taubin method,
//	as described in
//		Eyal Haimeri, Ilan Shimshoni
//		"Estimating teh Principal Curvatures and the Darboux Frame from Real 3D Range Data", 2002
//
//	uses the TNT and JAMA/C++ libraries for the eigenvalue/eigenvector calculations
//
// Created: Saturday, July 19, 2003, 10:15 AM - Nicu D. Cornea
//

#include "hcBoundHaimeiri3.h"

#include "tnt_array2d.h"
#include "jama_eig.h"


#define TRACE


#define EPSILON				0.0000001
#define IsZero(nr) 			(((nr) < EPSILON) && ((nr) > -EPSILON))
#define IsEqual(nr1, nr2)	(IsZero((nr1) - (nr2)))

#define MAX_HIGH_CURVATURE_POINTS	100


// returns k1 and k2 given e1, e2, A, B and C
bool FindK(double e1, double e2, double A, double B, double C, double *k1, double *k2) {
	if(!IsZero(B*B - A*C)) {
		*k1 = (B*e2 - C*e1) / (B*B - A*C);
	}
	else {
		if(!IsZero(A)) {
			*k1 = e1 / A;
		}
		else {
			if(	(IsZero(B*B - A*C) && !IsZero(B*e2 - C*e1)) ||
				(IsZero(A) && IsZero(B) && !IsZero(e1)) 	||
				(IsZero(B) && IsZero(C) && !IsZero(e2))
			  )
			{
				printf("UPS ! cannot solve for k1! - Abort!\n");
				exit(1);
			}
			else {
				*k1 = 999.99;
			}
		}
	}

	if(!IsZero(B)) {
		*k2 = (e1 - A*(*k1)) / B;
	}
	else {
		if(!IsZero(C)) {
			*k2 = e2 / C;
		}
		else {
			*k2 = 999.99;
		}
	}

	return true;
}


bool IsProportional(double v1x, double v1y, double v1z,
					double v2x, double v2y, double v2z)
{
	double k;

#ifdef TRACE
	// printf("Is proportional ?\n");
#endif

	if(!IsZero(v1x)) {
		k = v2x / v1x;

#ifdef TRACE
		//printf("v1x !=0\n\tk = %f\n", k);
#endif

		if(!IsZero(k)) {
			if(IsEqual(v2y, (v1y * k)) && IsEqual(v2z, (v1z * k)) ) {
				return true;
			}
		}
	}

	// v1x == 0
	// in that case v2x has to be = 0, for the 2 vectors to be proportional
	if(!IsZero(v2x)) {
		return false;
	}

	if(!IsZero(v1y)) {
			k = v2y / v1y;
			if(!IsZero(k)) {
				if(IsEqual(v2z, (v1z * k)) ) {
					return true;
				}
			}
	}

	// v1y == 0 too
	// in that case v2y has to be = 0, for the 2 vectors to be proportional
	if(!IsZero(v2y)) {
		return false;
	}

	//at this point the 2 vectors are proportional, regardless of the z value
	// they are of the following form: V1(0, 0, z1), V2(0, 0, z2)

	return true;
}


bool GetHighCurvaturePoints(
	bdElement* Bound, int numBound,
	VoxelPosition **HighCurvaturePoints, int *numHC,
	double perc)
{
	int i, j, k, kk;
	bool found, foundthenormal;

	Vector Ti, qmp, pd1, pd2, Wv;
	double sl, ki, length, wi, c, t, e1, e2, A, B, C, sinasq, cosa, k1, k2, maxk1, mink1;

	TNT::Array2D<double>	Ma(3, 3, 0.0);
	TNT::Array2D<double>	EigVals(3, 3, 0.0);
	TNT::Array2D<double>	EigVects(3, 3, 0.0);

#ifdef _DEBUG
	SetStartTime();
#endif

	(*numHC) = 0;
// get the maximum principal direction at each point on the boundary

	maxk1 = -999999999;
	mink1 = 999999999;
	for(i=0; i < numBound; i++) {
		Bound[i].curvature = 0.0;
#ifdef TRACE
		printf("\n\ni = %d (%d, %d, %d). neighbors: %d.\n", i, Bound[i].point.x, Bound[i].point.y, Bound[i].point.z, Bound[i].nrOfNeighbors);
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

		Ma[0][0] = 0.0;	Ma[0][1] = 0.0;	Ma[0][2] = 0.0;
		Ma[1][0] = 0.0;	Ma[1][1] = 0.0;	Ma[1][2] = 0.0;
		Ma[2][0] = 0.0;	Ma[2][1] = 0.0;	Ma[2][2] = 0.0;

		for(j=0; j < Bound[i].nrOfNeighbors; j++) {
			// estimating the ki term
#ifdef TRACE
			printf("\tNeighbor%d: (%d, %d, %d)\n", j, Bound[Bound[i].neighbors[j]].point.x, Bound[Bound[i].neighbors[j]].point.y, Bound[Bound[i].neighbors[j]].point.z);
#endif
			qmp.xd = Bound[Bound[i].neighbors[j]].point.x - Bound[i].point.x;
			qmp.yd = Bound[Bound[i].neighbors[j]].point.y - Bound[i].point.y;
			qmp.zd = Bound[Bound[i].neighbors[j]].point.z - Bound[i].point.z;
#ifdef TRACE
			printf("\t\tqmp: (%f, %f, %f)\n", qmp.xd, qmp.yd, qmp.zd);
#endif

			ki = (double) 2.0 *
				 (double) ((Bound[i].normal.xd * (qmp.xd)) + (Bound[i].normal.yd * (qmp.yd)) + (Bound[i].normal.zd * (qmp.zd))) /
				 (double) (qmp.xd * qmp.xd + qmp.yd * qmp.yd + qmp.zd * qmp.zd);
#ifdef TRACE
			printf("\t\tki = %f\n", ki);
#endif

			// calculation the normalized projection of qmp onto the tangent plane = Ti
			// project the neighbor on to the tangent plane
			t = -(Bound[i].normal.xd * qmp.xd + Bound[i].normal.yd * qmp.yd + Bound[i].normal.zd * qmp.zd);
#ifdef TRACE
			printf("\t\tt = %f\n", t);
#endif
			Ti.xd = t*Bound[i].normal.xd + qmp.xd;
			Ti.yd = t*Bound[i].normal.yd + qmp.yd;
			Ti.zd = t*Bound[i].normal.zd + qmp.zd;

			// normalize Ti
			length = sqrt(Ti.xd * Ti.xd + Ti.yd * Ti.yd + Ti.zd * Ti.zd);
			Ti.xd = Ti.xd / length;
			Ti.yd = Ti.yd / length;
			Ti.zd = Ti.zd / length;
#ifdef TRACE
			printf("\t\tprojection of qmp on the tangent plane: (%f, %f, %f)\n", Ti.xd, Ti.yd, Ti.zd);
#endif

			// the weights wi
			wi = ((double) 1.0 / sqrt(qmp.xd * qmp.xd + qmp.yd * qmp.yd + qmp.zd * qmp.zd)) / sl;
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

		// find the eigenvalues and eigenvectors of the Ma matrix
		//
		//	e1 = ...
		//	e2 = ...
		//
		// pd1 = ...
		// pd2 = ...
		// the third eigenvector is the normal.
		//
		//

		JAMA::Eigenvalue<double>		ev(Ma);
		ev.getD(EigVals);
		ev.getV(EigVects);

#ifdef TRACE
		printf("\tThe eigenvalues of Ma:\n");
		printf("\t\t%f  %f  %f\n", EigVals[0][0], EigVals[0][1], EigVals[0][2]);
		printf("\t\t%f  %f  %f\n", EigVals[1][0], EigVals[1][1], EigVals[1][2]);
		printf("\t\t%f  %f  %f\n", EigVals[2][0], EigVals[2][1], EigVals[2][2]);

		printf("\tThe eigenvectors of Ma:\n");
		printf("\t\t%f  %f  %f\n", EigVects[0][0], EigVects[0][1], EigVects[0][2]);
		printf("\t\t%f  %f  %f\n", EigVects[1][0], EigVects[1][1], EigVects[1][2]);
		printf("\t\t%f  %f  %f\n", EigVects[2][0], EigVects[2][1], EigVects[2][2]);
#endif

		// eigenvalues should be real ??
		// maybe do a check here...
#ifdef TRACE
		//printf("Looking for the normal among the eigenvectors...\n");
#endif
		// one of the eigenvalues is 0, and coresponds to the normal (which is an eigenvector)
		foundthenormal = false;
		for(k=0; k < 3; k++) {

#ifdef TRACE
			//printf("%d", k);
#endif

			if(IsZero(EigVals[k][k])) {
#ifdef TRACE
				//printf("\tnr %d  is = 0.0. Check the coresponding eigenvector\n", k);
#endif

				// check if the coresponding eigenvector is proportional with the normal
				if(IsProportional(Bound[i].normal.xd, Bound[i].normal.yd, Bound[i].normal.zd,
									EigVects[0][k], EigVects[1][k], EigVects[2][k]))
				{
					// save the other 2 eigenvalues in e1 and e2, and the other two eigenvectors in pd1 and pd2
					switch(k) {
						case 0:
							e1 = EigVals[1][1];
							e2 = EigVals[2][2];
							pd1.xd = EigVects[0][1];	pd1.yd = EigVects[1][1];	pd1.zd = EigVects[2][1];
							pd2.xd = EigVects[0][2];	pd2.yd = EigVects[1][2];	pd2.zd = EigVects[2][2];
							break;
						case 1:
							e1 = EigVals[0][0];
							e2 = EigVals[2][2];
							pd1.xd = EigVects[0][0];	pd1.yd = EigVects[1][0];	pd1.zd = EigVects[2][0];
							pd2.xd = EigVects[0][2];	pd2.yd = EigVects[1][2];	pd2.zd = EigVects[2][2];
							break;
						case 2:
							e1 = EigVals[0][0];
							e2 = EigVals[1][1];
							pd1.xd = EigVects[0][0];	pd1.yd = EigVects[1][0];	pd1.zd = EigVects[2][0];
							pd2.xd = EigVects[0][1];	pd2.yd = EigVects[1][1];	pd2.zd = EigVects[2][1];
							break;
						default:
							printf("UPS ! How did I get here ? - Abort.");
							exit(1);
							break;
					}
					foundthenormal = true;
					break;
				}
			}
		}


		if(!foundthenormal) {
			printf("UPS ! - Couldn't find the normal among the eigenvectors. Abort");
			printf("\tNormal: (%f, %f, %f)\n", Bound[i].normal.xd, Bound[i].normal.yd, Bound[i].normal.zd);

			printf("\tThe Ma matrix:\n");
			printf("\t\t%f  %f  %f\n", Ma[0][0], Ma[0][1], Ma[0][2]);
			printf("\t\t%f  %f  %f\n", Ma[1][0], Ma[1][1], Ma[1][2]);
			printf("\t\t%f  %f  %f\n", Ma[2][0], Ma[2][1], Ma[2][2]);

			printf("\tThe eigenvalues of Ma:\n");
			printf("\t\t%f  %f  %f\n", EigVals[0][0], EigVals[0][1], EigVals[0][2]);
			printf("\t\t%f  %f  %f\n", EigVals[1][0], EigVals[1][1], EigVals[1][2]);
			printf("\t\t%f  %f  %f\n", EigVals[2][0], EigVals[2][1], EigVals[2][2]);

			printf("\tThe eigenvectors of Ma:\n");
			printf("\t\t%f  %f  %f\n", EigVects[0][0], EigVects[0][1], EigVects[0][2]);
			printf("\t\t%f  %f  %f\n", EigVects[1][0], EigVects[1][1], EigVects[1][2]);
			printf("\t\t%f  %f  %f\n", EigVects[2][0], EigVects[2][1], EigVects[2][2]);

			exit(1);
		}

		if(IsZero(e1))	e1 = 0.0;
		if(IsZero(e2))	e2 = 0.0;

		if(IsZero(pd1.xd))	pd1.xd = 0.0;
		if(IsZero(pd1.yd))	pd1.yd = 0.0;
		if(IsZero(pd1.zd))	pd1.zd = 0.0;

		if(IsZero(pd2.xd))	pd2.xd = 0.0;
		if(IsZero(pd2.yd))	pd2.yd = 0.0;
		if(IsZero(pd2.zd))	pd2.zd = 0.0;

		// pd1 should be the eigenvector coresponding to the eigenvalue with the greatest absolute value
		if(fabs(e1) < fabs(e2)) {
			// swap the 2 values and vectors
			t = e1;			e1 = e2;			e2 = t;

			t = pd1.xd;		pd1.xd = pd2.xd;	pd2.xd = t;
			t = pd1.yd;		pd1.yd = pd2.yd;	pd2.yd = t;
			t = pd1.zd;		pd1.zd = pd2.zd;	pd2.zd = t;
		}

#ifdef TRACE
		printf("\tPrincipal directions: \n");
		printf("\t\te1 = %f, pd1 = (%f, %f, %f)\n", e1, pd1.xd, pd1.yd, pd1.zd);
		printf("\t\te2 = %f, pd2 = (%f, %f, %f)\n", e2, pd2.xd, pd2.yd, pd2.zd);
#endif

		// find the A, B, C coeficients
		A = 0;
		B = 0;
		C = 0;
		for(j=0; j < Bound[i].nrOfNeighbors; j++) {
			qmp.xd = Bound[Bound[i].neighbors[j]].point.x - Bound[i].point.x;
			qmp.yd = Bound[Bound[i].neighbors[j]].point.y - Bound[i].point.y;
			qmp.zd = Bound[Bound[i].neighbors[j]].point.z - Bound[i].point.z;

			// the weights wi
			wi = ((double) 1.0 / sqrt(qmp.xd * qmp.xd + qmp.yd * qmp.yd + qmp.zd * qmp.zd)) / sl;

			/*
			// normalize qmp
			length = sqrt(qmp.xd * qmp.xd + qmp.yd * qmp.yd + qmp.zd * qmp.zd);
			if(length == 0) {
				printf("UPS! - length of vector = 0! - Abort.\n");
				exit(1);
			}
			qmp.xd = qmp.xd / length;
			qmp.yd = qmp.yd / length;
			qmp.zd = qmp.zd / length;
			*/

			// calculation the normalized projection of qmp onto the tangent plane = Ti
			// project the neighbor on to the tangent plane
			t = -(Bound[i].normal.xd * qmp.xd + Bound[i].normal.yd * qmp.yd + Bound[i].normal.zd * qmp.zd);

			Ti.xd = t*Bound[i].normal.xd + qmp.xd;
			Ti.yd = t*Bound[i].normal.yd + qmp.yd;
			Ti.zd = t*Bound[i].normal.zd + qmp.zd;

			// normalize Ti
			length = sqrt(Ti.xd * Ti.xd + Ti.yd * Ti.yd + Ti.zd * Ti.zd);
			Ti.xd = Ti.xd / length;
			Ti.yd = Ti.yd / length;
			Ti.zd = Ti.zd / length;

			// cosa, sinasq = the angle between the pd1 and Ti vectors
			// cos = dot product; sinsq = 1 - cos*cos;
			cosa = Ti.xd*pd1.xd + Ti.yd*pd1.yd + Ti.zd*pd1.zd;
			sinasq = 1 - (cosa*cosa);

			A = A + wi*(cosa*cosa)*(cosa*cosa);
			B = B + wi*(sinasq)*(cosa*cosa);
			C = C + wi*(sinasq)*(sinasq);
		}

		if(IsZero(A))	A = 0.0;
		if(IsZero(B))	B = 0.0;
		if(IsZero(C))	C = 0.0;

#ifdef TRACE
		printf("\tA = %f, B = %f, C = %f\n", A, B, C);
#endif

		// the principal curvatures
		FindK(e1, e2, A, B, C, &k1, &k2);

		if(IsZero(k1))	k1 = 0.0;
		if(IsZero(k2))	k2 = 0.0;

#ifdef TRACE
		printf("\tk1 = %f", k1);
		printf("k2 = %f;\n", k2);
#endif
		//k1 = fabs(k1);
		//k2 = fabs(k2);

		if(fabs(k1) < fabs(k2)) {
			// THIS SHOULDN'T HAPPEN !!!
			printf("UPS! |k1| < |k2| !?? - Try again.\n");

			//printf("i = %d (%d, %d, %d)\n", i, Bound[i].point.x, Bound[i].point.y, Bound[i].point.z);
			//printf("\tNormal: (%f, %f, %f)\n", Bound[i].normal.xd, Bound[i].normal.yd, Bound[i].normal.zd);
			//printf("\tdelta = %f, e1 = %f, e2 = %f\n", delta, e1, e2);
			//exit(1);

			// try e2 and pd2.

			// swap the 2 values and vectors
			t = e1;			e1 = e2;			e2 = t;

			t = pd1.xd;		pd1.xd = pd2.xd;	pd2.xd = t;
			t = pd1.yd;		pd1.yd = pd2.yd;	pd2.yd = t;
			t = pd1.zd;		pd1.zd = pd2.zd;	pd2.zd = t;

			// and try again

#ifdef TRACE
			printf("\tpd1 = (%f, %f, %f)\n", pd1.xd, pd1.yd, pd1.zd);
#endif
			// find the A, B, C coeficients
			A = 0;
			B = 0;
			C = 0;
			for(j=0; j < Bound[i].nrOfNeighbors; j++) {
				qmp.xd = Bound[Bound[i].neighbors[j]].point.x - Bound[i].point.x;
				qmp.yd = Bound[Bound[i].neighbors[j]].point.y - Bound[i].point.y;
				qmp.zd = Bound[Bound[i].neighbors[j]].point.z - Bound[i].point.z;

				// the weights wi
				wi = ((double)1.0 / sqrt(qmp.xd * qmp.xd + qmp.yd * qmp.yd + qmp.zd * qmp.zd)) / sl;

				/*
				// normalize qmp
				length = sqrt(qmp.xd * qmp.xd + qmp.yd * qmp.yd + qmp.zd * qmp.zd);
				if(length == 0) {
					printf("UPS! - length of vector = 0! - Abort.\n");
					exit(1);
				}
				qmp.xd = qmp.xd / length;
				qmp.yd = qmp.yd / length;
				qmp.zd = qmp.zd / length;
				*/

				// calculation the normalized projection of qmp onto the tangent plane = Ti
				// project the neighbor on to the tangent plane
				t = -(Bound[i].normal.xd * qmp.xd + Bound[i].normal.yd * qmp.yd + Bound[i].normal.zd * qmp.zd);

				Ti.xd = t*Bound[i].normal.xd + qmp.xd;
				Ti.yd = t*Bound[i].normal.yd + qmp.yd;
				Ti.zd = t*Bound[i].normal.zd + qmp.zd;

				// normalize Ti
				length = sqrt(Ti.xd * Ti.xd + Ti.yd * Ti.yd + Ti.zd * Ti.zd);
				Ti.xd = Ti.xd / length;
				Ti.yd = Ti.yd / length;
				Ti.zd = Ti.zd / length;

				// cosa, sinasq = the angle between the pd1 and Ti vectors
				// cos = dot product; sinsq = 1 - cos*cos;
				cosa = Ti.xd*pd1.xd + Ti.yd*pd1.yd + Ti.zd*pd1.zd;
				sinasq = 1 - cosa*cosa;

				A = A + wi*(cosa*cosa)*(cosa*cosa);
				B = B + wi*(sinasq)*(cosa*cosa);
				C = C + wi*(sinasq)*(sinasq);
			}

			if(IsZero(A))	A = 0.0;
			if(IsZero(B))	B = 0.0;
			if(IsZero(C))	C = 0.0;

#ifdef TRACE
			printf("\tA = %f, B = %f, C = %f\n", A, B, C);
#endif

			// the principal curvatures
			FindK(e1, e2, A, B, C, &k1, &k2);

			if(IsZero(k1))	k1 = 0.0;
			if(IsZero(k2))	k2 = 0.0;

#ifdef TRACE
			printf("\tk1 = %f", k1);
			printf("k2 = %f;\n", k2);
#endif

			//k1 = fabs(k1);
			//k2 = fabs(k2);

			if(fabs(k1) < fabs(k2)) {
				printf("UPS! abs(k1) is still < abs(k2) !?? - Abort.\n");
				exit(1);
			}
		}

		Bound[i].curvature = k1;
		Bound[i].curvk2 = k2;

		if(!IsZero(k1)) {
			if(k1 > maxk1) maxk1 = k1;
			if(k1 < mink1) mink1 = k1;
		}
	}

#ifdef _DEBUG
	PrintElapsedTime("Finding principal curvature");
	printf(" -- Max curvature = %f, min curvature = %f\n", maxk1, mink1);
#endif

	// store the high curvature points into the HighCurvaturePoints array
	if(((*HighCurvaturePoints) = new VoxelPosition[MAX_HIGH_CURVATURE_POINTS]) == NULL) {
		printf("ERROR: allocating memory for high curvature points storage. Abort.\n");
		exit(1);
	}
	(*numHC) = 0;

	// write the high curvature points to the output file.
	// threshold = 1 % of the highest curvature value
	// output all the points with curvature value > threshold
	t = perc * (fabs(maxk1) / 100.00);
	t = maxk1 - t;

#ifdef TRACE
	printf("\tThreshold for outputed values: %f\n", t);
#endif

	for(i=0; i < numBound; i++) {
		if(!IsZero(Bound[i].curvature) && (Bound[i].curvature > t)) {
			(*HighCurvaturePoints)[(*numHC)].x = Bound[i].point.x;
			(*HighCurvaturePoints)[(*numHC)].y = Bound[i].point.y;
			(*HighCurvaturePoints)[(*numHC)].z = Bound[i].point.z;

			(*numHC) = (*numHC) + 1;
			if((*numHC) >= MAX_HIGH_CURVATURE_POINTS) {
				printf("UPS! Too many high curvature points selected (top %f percent). Stopping at %d\n", perc, MAX_HIGH_CURVATURE_POINTS);
				break;
			}
		}
	}

#ifdef _DEBUG
	PrintElapsedTime("Phase 3: generating output.");
#endif

	return true;
}

