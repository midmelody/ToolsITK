// Find critical points of a vector field
// --- Input: 1. normalized 3D vector field
//
// --- Output: critical point list
// --- Author: Nicu D. Cornea, Vizlab, Rutgers University
// --- Date: Fri Aug  1 13:14:18 EDT 2003
//

#include "CritPts.h"


//#define TRACE

#define CELL_SUBDIVISION_FACTOR	1048576.00

#define MAX_NUM_CRITPTS	2000

#define NR_OF_SIGN_CHANGES	3


// #define SIGN(nr)		(((nr) < 0.00) ? -1 : 1)
// #define SIGN(nr)		(((nr) < 0.00) ? -1 : (((nr) > 0.00) ? 1 : 0))
#define SIGN(nr)		(IS_ZERO(nr) ? 0 : ((nr) < 0.00 ? -1 : 1))


inline Vector interpolation(double x, double y, double z, 
			    int sizx, int sizy, int sizz, Vector *forcevec);

bool FindCriticalPointInFloatCell(double x, double y, double z, 
				  double cellsize,
				  int sX, int sY, int sZ, 
				  Vector* ForceField, 
				  VoxelPositionDouble *critPt);
bool FindCriticalPointInIntCell(int x, int y, int z, int* inds,
				int sX, int sY, int sZ, Vector* ForceField, 
				VoxelPositionDouble *critPt);

bool ChangeInSign(Vector* forces, int numForces);
bool ChangeInSign(Vector* forces, int *indices, int numIndices);

bool GetCriticalPoints(
	Vector* ForceField,         // [in] vector field
	int L, int M, int N,        // [in] volume size
	unsigned char *flags,       // [in] volume flags
	CriticalPoint **CritPts,    // [out] detected critical points
	int *numCritPts,            // [out] critical points count
	bool inOut                  // [in] flag specifying if we know the 
                                    //    inside of the object from the outside
                                    //    false - can distinguish the inside
	                            //    true - cannot distinguish the inside
	                            //    DEFAULF: false
) {


#ifdef TRACE
  printf("TRACE: Starting GetCriticalPoints function... flags = %p\n", flags);
#endif

  // Vector *Normforce;
  
  long idx, slsz;
  int i,j,k, ii;
  
  double vecLength;
  VoxelPositionDouble critPt;
  
  
  slsz = L*M;		// slice size
  
  
  (*CritPts) = NULL;
  (*numCritPts) = 0;
  
  if(((*CritPts) = new CriticalPoint[MAX_NUM_CRITPTS]) == NULL) {
    printf("UPS! - Error allocating memory for the output array. Abort.\n");
    exit(1);
  }
  
  
  // find critical points
  
  (*numCritPts) = 0;
  int inds[27];
  bool skipThisPoint = false;
  
#ifdef _DEBUG
  printf("\
** Critical points on the boundary \n\
      or a voxel away from the boundary are ignored **\n"); 
  if(inOut) {
    printf("** Inside and Outside. **\n"); 
  }
  else {
    printf("** Inside ONLY. **\n"); 
  }
#endif

/////////////////////////////////////

  for (k = 1; k < N-1; k++) {
    printf("\tProcessing plane %d out of %d\r", k, N-1);
    fflush(stdout);
    
    for (j = 1; j < M-1; j++) {
      for (i = 1; i < L-1; i++) {
	
	idx = k*slsz + j*L +i;

	// - ignore voxels that are BOUNDARY or SURF or have a neighbor that is
	//
	// - also, if we know the inside vs/ the outside of the object, ignore 
	//    voxels that are EXTERIOR or have a neighbor that is

	// - we have to check the neighbors that will be touched by the 
	//   interpolation

	// - the flags array indicates the state of all the face neighbors,
	//   but we are not interested in all the face neighbors, 
	//   just 3 of them, so we will not be using the flag SURF or BOUNDARY 
	//   at the current point to find out if it has any external neighbors.
	
	// these are the neighbors that will be touched by the interpolation
	inds[0] = idx;
	inds[1] = idx + 1;
	inds[2] = idx + L;
	inds[3] = idx + L + 1;
	inds[4] = idx + slsz;
	inds[5] = idx + slsz + 1;
	inds[6] = idx + slsz + L;
	inds[7] = idx + slsz + L + 1;

	skipThisPoint = false;

	if(!inOut) {
	  // if we know the inside from the outside, 
	  for(ii=0; ii < 8; ii++) { 
	    // ignore voxels close to the boundary
	    if(ii > 0) {
	      if((flags[inds[ii]] == SURF) ||
		 (flags[inds[ii]] == BOUNDARY)) 
	      {
		skipThisPoint = true;
		break;
	      }
	    }

	    // ignore voxels on the boundary (neighbor to an exterior point)
	    // or in the exterior 
	    if(flags[inds[ii]] == EXTERIOR) {
	      skipThisPoint = true;
	      break;
	    }
	  }
	}
	else {
	  // printf("here\n");
	  // we don't know the interior from the exterior
	  
	  // ignore points close to the boundary
	  for(ii=0; ii < 8; ii++) { 
	    // ignore voxels close to the boundary
	    if((flags[inds[ii]] == SURF) ||
	       (flags[inds[ii]] == BOUNDARY)) 
	    {
	      skipThisPoint = true;
	      break;
	    }
	  }
	  
	  /*
	  // ignore only boundary points
	  if((flags[idx] == SURF) ||
	     (flags[idx] == BOUNDARY)) 
	  {
	    skipThisPoint = true;
	  }
	  */
	  
	}
	
	
	if(skipThisPoint)  continue;

	if(
	   FindCriticalPointInIntCell(i, j, k, inds, L, M, N, 
				      ForceField, &critPt)) 
	{
	  (*CritPts)[(*numCritPts)].position.x = critPt.x;
	  (*CritPts)[(*numCritPts)].position.y = critPt.y;
	  (*CritPts)[(*numCritPts)].position.z = critPt.z;
	  (*CritPts)[(*numCritPts)].type = CPT_UNKNOWN;
	  (*numCritPts) = (*numCritPts) + 1;
	  if((*numCritPts) >= MAX_NUM_CRITPTS) {
	    printf("UPS! - Too many critical points found! - Abort.\n");
	    exit(1);
	  }
	}
      }
    }
  }

#ifdef _DEBUG
  PrintElapsedTime("\n\tCP-2: Locating critical points.");
  printf("\t** Number of critical points is: %d\n", (*numCritPts));
#endif

#ifdef TRACE
  // print the critical points
  printf("Critical points found:\n");
  for(i = 0; i < (*numCritPts); i++) {
    printf("%f %f %f\n", 
	   (*CritPts)[i].position.x, 
	   (*CritPts)[i].position.y, 
	   (*CritPts)[i].position.z);
  }
#endif

  // clasify the critical points as: atracting/repelling nodes or saddles

  Vector cv[6];

  // distance to a neighbor will be  1 / CELL_SUBDIVISION_FACTOR

  double vdist = 1.00 / CELL_SUBDIVISION_FACTOR;

  TNT::Array2D<double> 	Jac(3, 3, 0.00);
  TNT::Array2D<double>	EigVals(3, 3, 0.0);
  TNT::Array2D<double>	EigVects(3, 3, 0.0);
  
  // since critical points are at least 1 voxel inside the boundaing box, 
  //	we will not check that the neighbor coordinates are inside the 
  //    volume bounds, because they should be
  
  for(i=0; i < (*numCritPts); i++) {
    // interpolate the force vector at 6 neighbors of this point.

#ifdef TRACE
    printf("Point %d:\n", i);
#endif

    cv[0] = interpolation((*CritPts)[i].position.x + vdist, 
			  (*CritPts)[i].position.y, 
			  (*CritPts)[i].position.z,
			  L, M, N, ForceField);
    cv[1] = interpolation((*CritPts)[i].position.x - vdist, 
			  (*CritPts)[i].position.y, 
			  (*CritPts)[i].position.z,
			  L, M, N, ForceField);
    cv[2] = interpolation((*CritPts)[i].position.x, 
			  (*CritPts)[i].position.y + vdist, 
			  (*CritPts)[i].position.z,
			  L, M, N, ForceField);
    cv[3] = interpolation((*CritPts)[i].position.x, 
			  (*CritPts)[i].position.y - vdist, 
			  (*CritPts)[i].position.z,
			  L, M, N, ForceField);
    cv[4] = interpolation((*CritPts)[i].position.x, 
			  (*CritPts)[i].position.y, 
			  (*CritPts)[i].position.z + vdist,
			  L, M, N, ForceField);
    cv[5] = interpolation((*CritPts)[i].position.x, 
			  (*CritPts)[i].position.y, 
			  (*CritPts)[i].position.z - vdist,
			  L, M, N, ForceField);

#ifdef TRACE
    for(j=0; j < 6; j++) {
      printf("\t\t(%lf, %lf, %lf)\n", cv[j].xd, cv[j].yd, cv[j].zd);
    }
#endif

    // construct the Jacobian matrix
    // is central differencing ok ???
    Jac[0][0] = (cv[0].xd - cv[1].xd) / (2 * vdist);
    Jac[0][1] = (cv[2].xd - cv[3].xd) / (2 * vdist);
    Jac[0][2] = (cv[4].xd - cv[5].xd) / (2 * vdist);
    
    Jac[1][0] = (cv[0].yd - cv[1].yd) / (2 * vdist);
    Jac[1][1] = (cv[2].yd - cv[3].yd) / (2 * vdist);
    Jac[1][2] = (cv[4].yd - cv[5].yd) / (2 * vdist);
    
    Jac[2][0] = (cv[0].zd - cv[1].zd) / (2 * vdist);
    Jac[2][1] = (cv[2].zd - cv[3].zd) / (2 * vdist);
    Jac[2][2] = (cv[4].zd - cv[5].zd) / (2 * vdist);
    
    // find the eigenvalues and eigenvectors of the Jacobian
    
    JAMA::Eigenvalue<double>		ev(Jac);
    ev.getD(EigVals);
    ev.getV(EigVects);
    
#ifdef TRACE
    printf("\tThe eigenvalues of Jac:\n");
    printf("\t\t%lf  %lf  %lf\n", EigVals[0][0], EigVals[0][1], EigVals[0][2]);
    printf("\t\t%lf  %lf  %lf\n", EigVals[1][0], EigVals[1][1], EigVals[1][2]);
    printf("\t\t%lf  %lf  %lf\n", EigVals[2][0], EigVals[2][1], EigVals[2][2]);
    
    printf("\tThe eigenvectors of Jac:\n");
    printf("\t\t%lf  %lf  %lf\n", EigVects[0][0], EigVects[0][1], EigVects[0][2]);
    printf("\t\t%lf  %lf  %lf\n", EigVects[1][0], EigVects[1][1], EigVects[1][2]);
    printf("\t\t%lf  %lf  %lf\n", EigVects[2][0], EigVects[2][1], EigVects[2][2]);
#endif

    // analyze the eigenvalues:

    // all real parts of the eigenvalues are negative: - attracting node
    if(	(EigVals[0][0] < 0.00) &&
	(EigVals[1][1] < 0.00) &&
	(EigVals[2][2] < 0.00))
      {
	(*CritPts)[i].type = CPT_ATTRACTING_NODE;
#ifdef TRACE
	printf("Attracting node.\n");
#endif
      }
    else {
      // all real parts of the eigenvalues are pozitive: - repelling node
      if(	(EigVals[0][0] > 0.00) &&
		(EigVals[1][1] > 0.00) &&
		(EigVals[2][2] > 0.00))
	{
	  (*CritPts)[i].type = CPT_REPELLING_NODE;
#ifdef TRACE
	  printf("Repelling node.\n");
#endif
	  
	}
      else {
	// two real parts of the eigenvalues are of one sign 
	//   and the other has the opposite sign
	//	- this is a saddle point
	(*CritPts)[i].type = CPT_SADDLE;
#ifdef TRACE
	printf("Saddle point.\n");
#endif
      }
    }
    

    (*CritPts)[i].eval[0] = EigVals[0][0];
    (*CritPts)[i].eval[1] = EigVals[1][1];
    (*CritPts)[i].eval[2] = EigVals[2][2];
    
    // normalize the eigenvectors
    vecLength = sqrt(EigVects[0][0]*EigVects[0][0] + 
		     EigVects[1][0]*EigVects[1][0] + 
		     EigVects[2][0]*EigVects[2][0]);
    if(vecLength == 0.00) {
      vecLength = 1.00;
    }
    (*CritPts)[i].evect[0].xd = EigVects[0][0] / vecLength;
    (*CritPts)[i].evect[0].yd = EigVects[1][0] / vecLength;
    (*CritPts)[i].evect[0].zd = EigVects[2][0] / vecLength;
    
    vecLength = sqrt(EigVects[0][1]*EigVects[0][1] + 
		     EigVects[1][1]*EigVects[1][1] + 
		     EigVects[2][1]*EigVects[2][1]);
    if(vecLength == 0.00) {
      vecLength = 1.00;
    }
    (*CritPts)[i].evect[1].xd = EigVects[0][1] / vecLength;
    (*CritPts)[i].evect[1].yd = EigVects[1][1] / vecLength;
    (*CritPts)[i].evect[1].zd = EigVects[2][1] / vecLength;
    
    vecLength = sqrt(EigVects[0][2]*EigVects[0][2] + 
		     EigVects[1][2]*EigVects[1][2] + 
		     EigVects[2][2]*EigVects[2][2]);
    if(vecLength == 0.00) {
      vecLength = 1.00;
    }
    (*CritPts)[i].evect[2].xd = EigVects[0][2] / vecLength;
    (*CritPts)[i].evect[2].yd = EigVects[1][2] / vecLength;
    (*CritPts)[i].evect[2].zd = EigVects[2][2] / vecLength;
  }


#ifdef _DEBUG
  PrintElapsedTime("\tCP-3: Characterizing critical points.");
#endif
  
  return true;
}


inline Vector interpolation(double x, double y, double z, 
			    int sizx, int sizy, int sizz, Vector *forcevec)
{
  float alpha, beta, gamma;
  Vector forceInt;
  long slsz;
  
  alpha=x-int(x);
  beta=y-int(y);
  gamma=z-int(z);
  slsz=sizy*sizx;
  
  forceInt.xd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*(1-gamma)
    +forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*alpha*beta*(1-gamma)
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*gamma
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*(alpha*beta*gamma);
  
  forceInt.yd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*(1-gamma)
    +forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*(1-gamma)
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*gamma
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*gamma;
  
  forceInt.zd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*(1-gamma)
    +forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*(1-gamma)
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*gamma
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*gamma;
  
  return(forceInt);
}



// Sign change notes:
//
// if one component is 0 in all vectors, it should be considered as a 
//   change in sign
//   or else, we might miss some critical points for which we have only
//   NR_OF_SIGN_CHANGES - 1 changes in sign and the last component is 0
//   all over the cell. Example (2D):
//     *->--<-*
//     *->--<-*
// In the above example, the arrows indicate the the vectors at the 4 
//    corners of a 2D cell. The X component changes sign, but the Y 
//    component is 0 in all 4 corners, and thus we might miss this 
//    critical point if we require exactly 2 sign changes. Considering
//    0 everywhere as a change in sign, will include this critical point.
//


bool ChangeInSign(Vector* forces, int numForces) {
  int s, i;
  int count;
  bool change;
  
  // SEE notes on sign change above

  // change in sign for at least one of the vector components leads to a 
  //    sort of medial surface
  //
  // change in sign for at least NR_OF_SIGN_CHANGES vector components
  //
  
  if(numForces > 0) {
    count = 0;
    // check the components of these vectors.
    
    // X
    change = false;
    s = SIGN(forces[0].xd);
    for(i=1; i<numForces; i++) {
      if(SIGN(forces[i].xd) != s) {
	count = count + 1;
	change = true;
	break;
      }
    }
    
    if((!change) && (s == 0)) {
      // no change in sign but sign is 0
      count = count + 1;
      // printf("Yes\n");
    }
    
    if(count >= NR_OF_SIGN_CHANGES) {
      // return if at least NR_OF_SIGN_CHANGES vector components change sign
      return true;
    }

    // Y
    change = false;
    s = SIGN(forces[0].yd);
    for(i=1; i<numForces; i++) {
      if(SIGN(forces[i].yd) != s) {
	count = count + 1;
	change = true;
	break;
      }
    }

    if((!change) && (s == 0)) {
      // no change in sign but sign is 0
      count = count + 1;
      //printf("Yes\n");
    }
    
    if(count >= NR_OF_SIGN_CHANGES) {
      // return if at least NR_OF_SIGN_CHANGES vector components change sign
      return true;
    }
    
    // Z
    change = false;
    s = SIGN(forces[0].zd);
    for(i=1; i<numForces; i++) {
      if(SIGN(forces[i].zd) != s) {
	count = count + 1;
	change = true;
	break;
      }
    }
    
    if((!change) && (s == 0)) {
      // no change in sign but sign is 0
      count = count + 1;
      //printf("Yes\n");
    }

    if(count >= NR_OF_SIGN_CHANGES) {
      // return if at least NR_OF_SIGN_CHANGES vector components change sign
      return true;
    }
  }
  
  return false;
}


bool ChangeInSign(Vector* forces, int *indices, int numIndices) {
  int s, i;
  int count;
  bool change;
  
  //
  // SEE notes on sign change above
  //
  // change in sign for at least NR_OF_SIGN_CHANGES vector components
  if(numIndices > 0) {
    count = 0;
    // check the components of these vectors.
    
    // X
    change = false;
    s = SIGN(forces[indices[0]].xd);
    for(i=1; i<numIndices; i++) {
      if(SIGN(forces[indices[i]].xd) != s) {
	count = count + 1;
	change = true;
	break;
      }
    }
    
    if((!change) && (s == 0)) {
      // no change in sign but sign is 0
      //printf("Yes\n");
      count = count + 1;
    }
    
    if(count >= NR_OF_SIGN_CHANGES) {
      // return if at least NR_OF_SIGN_CHANGES vector components change sign
      return true;
    }

    // Y
    change = false;
    s = SIGN(forces[indices[0]].yd);
    for(i=1; i<numIndices; i++) {
      if(SIGN(forces[indices[i]].yd) != s) {
	count = count + 1;
	change = true;
	break;
      }
    }
    
    if((!change) && (s == 0)) {
      // no change in sign but sign is 0
      count = count + 1;
      //printf("Yes\n");
    }

    if(count >= NR_OF_SIGN_CHANGES) {
      // return if at least NR_OF_SIGN_CHANGES vector components change sign
      return true;
    }
    
    // Z
    s = SIGN(forces[indices[0]].zd);
    change = false;;
    for(i=1; i<numIndices; i++) {
      if(SIGN(forces[indices[i]].zd) != s) {
	count = count + 1;
	change = true;
	break;
      }
    }
    
    if((!change) && (s == 0)) {
      // no change in sign but sign is 0
      count = count + 1;
      //printf("Yes\n");
    }

    if(count >= NR_OF_SIGN_CHANGES) {
      // return if at least NR_OF_SIGN_CHANGES vector components change sign
      return true;
    }
  }
  
  return false;
}

bool FindCriticalPointInIntCell(
	int x, int y, int z, int* inds,
	int sX, int sY, int sZ, 
	Vector* ForceField, VoxelPositionDouble *critPt)
{

  int kk, jj, ii;
#ifdef TRACE
  printf("Testing cell: %d, %d, %d\n", x, y, z);
  printf("\tForces:\n");
  for(ii=0; ii < 8; ii++) {
    printf("\t\t(%f,\t%f,\t%f)\n", 
	   ForceField[inds[ii]].xd, ForceField[inds[ii]].yd, 
	   ForceField[inds[ii]].zd);
  }
#endif

  // if first vertex vector (coresponding to (x, y, z) with index in inds[0]) 
  //    is (0, 0, 0)
  //	then return (x, y, z) as the critical point.
  if(	(IS_ZERO(ForceField[inds[0]].xd)) &&
	(IS_ZERO(ForceField[inds[0]].yd)) &&
	(IS_ZERO(ForceField[inds[0]].zd)))
    {
      // printf("\tFirst vertex vector = 0. Found a critical point.\n");
      
      (*critPt).x = x;
      (*critPt).y = y;
      (*critPt).z = z;
      return true;
    }
  
  // the cell is a candidate cell if there is a change of sign
  //	in one of the vector components among all eight vertices of the cell
  if(ChangeInSign(ForceField, inds, 8)) {
    // it is a candidate cell
    // divide the cell in 8 subcells
    // for each of those 8 subcells do the candidate test again
    // and try to find the critical point in one of the candidate subcells
    // printf("\tCandidate cell\n");
    
    for(kk=0; kk < 2; kk++) {
      for(jj=0; jj < 2; jj++) {
	for(ii=0; ii < 2; ii++) {
	  if(
	     FindCriticalPointInFloatCell(
	        x + ii/2.00, y + jj/2.00, z + kk/2.00, 0.50,
		sX, sY, sZ, ForceField, critPt))
	  {
	    // crtPt is already set
	    return true;
	  }
	}
      }
    }
  }

  return false;
}


bool FindCriticalPointInFloatCell(
	double x, double y, double z, double cellsize,
	int sX, int sY, int sZ, 
	Vector* ForceField, VoxelPositionDouble *critPt)
{
  int kk, jj, ii;
  Vector cv[8];
  
#ifdef TRACE	
  printf("\ttesting subcell: %f, %f, %f. Cell size: %f\n", x, y, z, cellsize);
#endif
  // interpolate vector values at each of the 8 vertices of the cell
  cv[0] = interpolation(x, 		y, 		z,
			sX, sY, sZ, ForceField);
  cv[1] = interpolation(x+cellsize, 	y, 		z,
			sX, sY, sZ, ForceField);
  cv[2] = interpolation(x,		y+cellsize,	z,
			sX, sY, sZ, ForceField);
  cv[3] = interpolation(x+cellsize,	y+cellsize,	z,
			sX, sY, sZ, ForceField);
  cv[4] = interpolation(x,		y,		z+cellsize,
			sX, sY, sZ, ForceField);
  cv[5] = interpolation(x+cellsize,	y,		z+cellsize,
			sX, sY, sZ, ForceField);
  cv[6] = interpolation(x,		y+cellsize,	z+cellsize,	
			sX, sY, sZ, ForceField);
  cv[7] = interpolation(x+cellsize,	y+cellsize,	z+cellsize,	
			sX, sY, sZ, ForceField);

#ifdef TRACE	
  printf("\t\tForces:\n");
  for(int ii=0; ii < 8; ii++) {
    printf("\t\t\t(%f,\t%f,\t%f)\n", cv[ii].xd, cv[ii].yd, cv[ii].zd);
  }
#endif

  
  // if first vertex vector (coresponding to (x, y, z)) is (0, 0, 0)
  //	then return (x, y, z) as the critical point.
  if(	(IS_ZERO(cv[0].xd)) &&
	(IS_ZERO(cv[0].yd)) &&
	(IS_ZERO(cv[0].zd)))
    {
      // printf("\tFirst vertex vector = 0. Found a critical point.\n");
      
      (*critPt).x = x;
      (*critPt).y = y;
      (*critPt).z = z;
      return true;
    }
  
  
  // the cell is a candidate cell if there is a change of sign
  //	in one of the vector components among all eight vertices of the cell
  if(ChangeInSign(cv, 8)) {
    
    // if cell size < 1/CELL_SUBDIVISION_FACTOR,
    //	stop here and assume critical point is in the center of the cell
    if(cellsize <= (1.00 / CELL_SUBDIVISION_FACTOR)) {
      
      (*critPt).x = x + (cellsize / 2.00);
      (*critPt).y = y + (cellsize / 2.00);
      (*critPt).z = z + (cellsize / 2.00);
#ifdef TRACE
      Vector a;
      double l;
      a = interpolation((*critPt).x, (*critPt).y, (*critPt).z, sX, sY, sZ, ForceField);
      l = sqrt(a.xd*a.xd + a.yd*a.yd + a.zd*a.zd);
      printf("\n\tCell size too small. Assume critical point is in the middle. \n");
      printf("\tVector lenght = %f\n", l);
      
      printf("\t\tForces:\n");
      for(int ii=0; ii < 8; ii++) {
	printf("\t\t\t(%lf,\t%lf,\t%lf)\n", cv[ii].xd, cv[ii].yd, cv[ii].zd);
      }
#endif
      
      return true;
    }
    
    
    // it is a candidate cell and it's not small enough
    // divide the cell in 8 subcells
    // and try to find the critical point in one of the subcells
    // printf("\tCandidate subcell\n");
    for(kk=0; kk < 2; kk++) {
      for(jj=0; jj < 2; jj++) {
	for(ii=0; ii < 2; ii++) {
	  if(
	     FindCriticalPointInFloatCell(
	        x + ii*cellsize/2.00, y + jj*cellsize/2.00, 
		z + kk*cellsize/2.00, cellsize/2.00,
		sX, sY, sZ, ForceField, critPt))
	    {
	      // crtPt is already set
	      return true;
	    }
	}
      }
    }
  }
  
  return false;
}


