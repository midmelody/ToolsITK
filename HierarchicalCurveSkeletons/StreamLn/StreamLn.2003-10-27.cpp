// Form streamlines from critical points and seed points
//
// --- Author: Xiaosong Yuan, Bala, Nicu D. Cornea - Vizlab, Rutgers University
// --- Date: 9/28/2002
// --- Original file name: skel_streamline_float.cpp
//
// --- Changed to StreamLn.cpp by Nicu D. Cornea, on Tuesday, July 22, 2003
//
//

/*
Is it possible that I have to split a segment even if I'm very close
to one of it's end points. In that case, one of the resulting segments will 
be too small to fit the skeleton segment pattern.
*/

// #define TRACE

#include "StreamLn.h"

#define MAX_NUM_SKEL	       80000
#define MAX_NUM_SKEL_SEGMENTS  1000

#define SP_SMALL_DISTANCE	0.5		// defines "close to ..."
	// in terms of Manhattan distance |x1-x2| + |y1-y2| + |z1-z2|
	//

#define HD_SMALL_DISTANCE	2.00	// defines how close a high divergence
        // point should be
	// to the skeleton, for it to be ignored

#define SKEL_SEGMENTS_INC      5


#define  MIN_SEG_LENGTH  5  // minimum segment length (number of points)
   // only segments containing at least that many points will be included 
   //   in the skeleton. However, if the segment originated in a critical 
   //   point, we should keep it regardless of the length. Removing it, 
   //   might disconnect the basic skeleton.

#define SP_CLOSE_TO_SEGMENT_ENDPOINT  5  // defines the maximum number of 
   // skeleton points that can separate an intersection with a 
   // skeleton segment from the segment's end points so that the 
   // intersection is not reported. 
   // When constructing a skeleton segment, at each step we test for
   // intersection with the other existing segments. If we are close to 
   // another skeleton segment, we also check if we are also close to
   // one of that segment's end points. If we are, we keep going even if
   // the points in the 2 segments overlap, in the hope that we can end the 
   // current segment at the same point as the segment we are intersecting,
   // thus reducing the number of joints. Of course this does not work if
   // we are close to an end point but we are actualy moving towards the
   // other wnd of the segment, but in that case, the current segment will
   // be terminated once we get far enough from the end point.


int AddNewSkelSegment(Skeleton *Skel, int left, int first);
bool SegmentIsLongEnough(Skeleton *Skel, int crtSegment);
int  GetSkelSegment(Skeleton *Skel, int point, bool *endpoint);
bool DeleteLastSkelSegment(Skeleton *Skel);

inline Vector interpolation(double x, double y, double z, int sizx, int sizy,
			    int sizz, Vector *forcevec);
inline void rk2(double x, double y, double z, int sizx, int sizy, int sizz,
		double steps,
		Vector *Force_ini, VoxelPositionDouble *nextPos);

bool FollowStreamlines(
        int originCP,                 
	int originSkel,
	bool lookInCP, 
	double x, double y, double z,
	int sX, int sY, int sZ, Vector* ForceField,
	CriticalPoint *CritPts, int numCritPts,
	int *critPtInSkeleton,
	Skeleton *Skel);

inline double Distance(	VoxelPositionDouble *pos1,
			VoxelPositionDouble *pos2);

int CloseToSkel( VoxelPositionDouble *pos, int originSkel,
		 Skeleton *Skel, int nrAlreadyInSkel, 
		 double maxDistance, int crtSegmentLength);

int CloseToCP ( VoxelPositionDouble *pos, int	originCP, 
		CriticalPoint *CritPts, int numCritPts, 
		double maxDistance);

int AddSkeletonPoint(Skeleton *Skel, VoxelPositionDouble *pos);



//////////////////////////////////////////////////////////////////////////////
// Get basic skeleton (level 1): connects only critical points
//////////////////////////////////////////////////////////////////////////////
bool GetLevel1Skeleton(
        Vector* ForceField, 		// [in] vector field
	int L, int M, int N,		// [in] vector field size (X, Y and Z)
	CriticalPoint *CritPts,         // [in] critical points array
	int numCritPts, 		// [in] number of critical points
	Skeleton *Skel                  // [out] skeleton points
) {

#ifdef TRACE
  printf("TRACE: Starting GetLevel1Skeleton function...\n");
#endif

  if(Skel == NULL) {
    printf("Argument 1 (*Skel) to GetLevel1Skeleton cannot be NULL! Abort.\n");
    exit(1);
  }


  long idx, iidx, slsz, prvidx;
  int cc;
  int ii,jj,kk;
  int i,j,k;
  slsz = L*M;		// slice size

  int *critPtInSkeleton;
  int posInSkel;

  //
  // critPtInSkeleton indicates if a critical point is already part of the 
  //   skeleton (is in the Skel->Points array). If -1 then, the point is not 
  //   yet in the Points array. Otherwise it indicates the point's position in 
  //   the skeleton where the point can be found.
  //

  if(numCritPts < 1) {
    return true;
  }

  if((critPtInSkeleton = new int[numCritPts]) == NULL) {
    printf("UPS! - Error allocating memory for working data structures. \
              Abort.\n");
    exit(1);
  }
  // initialize to -1.
  for(i=0; i < numCritPts; i++) {
    critPtInSkeleton[i] = -1;
  }
  
  
#ifdef TRACE
  printf("Critical points:\n");
  for(i=0; i < numCritPts; i++) {
    printf("%d: (%f %f %f)\n", i, CritPts[i].position.x,
	   CritPts[i].position.y, CritPts[i].position.z);
  }
#endif

  //
  // follow the streamlines starting at saddles in the direction of the
  // positive eigenvector(s)
  //	until we are close to one of the skeleton points or a critical point,
  //    ignoring the points in the current segment.
  //

  for(i = 0; i < numCritPts; i++) {
     // start to follow force field at saddle points only
    if(CritPts[i].type != CPT_SADDLE) {
      continue;
    }

    idx = (int)CritPts[i].position.z * slsz + (int)CritPts[i].position.y *L
      + (int)CritPts[i].position.x;
    
#ifdef TRACE
    printf("Point: (%f, %f, %f): idx = %d, force here: (%f, %f, %f)\n",
	   CritPts[i].position.x, CritPts[i].position.y, CritPts[i].position.z,
	   idx,
	   ForceField[idx].xd, ForceField[idx].yd, ForceField[idx].zd);
#endif
    

#ifdef TRACE
    printf("Starting to follow the force field...\n");
#endif

    // the point is a saddle, so we should go in the direction pointed by
    // the pozitive eigenvectors
    for(j=0; j < 3; j++) {
      if(CritPts[i].eval[j] > 0) {
	

	//
	// if the critical point is not already in the skeleton, 
	//   add it as belonging to the current segment
	//
	posInSkel = critPtInSkeleton[i];
	if(posInSkel == -1) {
	  posInSkel = AddSkeletonPoint(Skel, &(CritPts[i].position));
	  //
	  // mark the point as being part of the skeleton already
	  //
	  critPtInSkeleton[i] = posInSkel;
	}

	// direction given by the eigenvector
	FollowStreamlines(
	  i, posInSkel, 
	  true, 
	  CritPts[i].position.x + (CritPts[i].evect[j].xd * SP_SMALL_DISTANCE),
	  CritPts[i].position.y + (CritPts[i].evect[j].yd * SP_SMALL_DISTANCE),
	  CritPts[i].position.z + (CritPts[i].evect[j].zd * SP_SMALL_DISTANCE),
	  L, M, N, ForceField,
	  CritPts, numCritPts,
	  critPtInSkeleton,
	  Skel);
	
	// the exact opposite direction of the eigenvector
	FollowStreamlines(
	  i, posInSkel, 
	  true,
	  CritPts[i].position.x - (CritPts[i].evect[j].xd * SP_SMALL_DISTANCE),
	  CritPts[i].position.y - (CritPts[i].evect[j].yd * SP_SMALL_DISTANCE),
	  CritPts[i].position.z - (CritPts[i].evect[j].zd * SP_SMALL_DISTANCE),
	  L, M, N, ForceField,
	  CritPts, numCritPts,
	  critPtInSkeleton,
	  Skel);
      }
    }
  }

#ifdef _DEBUG
  PrintElapsedTime("\tSL-2: Following streamlines starting at the critical points.");
#endif

  return true;
}


bool GetStreamLines(
	Vector* ForceField, int L, int M, int N,
	CriticalPoint *CritPts, int numCritPts,
	VoxelPosition *BdSeeds, int numBoundSeeds,
	VoxelPositionDouble *HDPoints, int numHDPoints,
        Skeleton *Skel)
{

#ifdef TRACE
  printf("TRACE: Starting GetStreamLines function...\n");
#endif

  if(Skel == NULL) {
    printf("Argument 1 (*Skel) to GetStreamlines cannot be NULL! Abort.\n");
    exit(1);
  }
  
  GetLevel1Skeleton( ForceField, L, M, N,
		     CritPts, numCritPts, 	
		     Skel);

  GetLevel2Skeleton( ForceField, L, M, N,
		     HDPoints, numHDPoints, 	
		     Skel);

  return true;
}







inline Vector interpolation(double x, double y, double z, int sizx, int sizy,
int sizz, Vector *forcevec)
    {
	float alpha, beta, gamma;
	Vector forceInt;
	long slsz;

	alpha=x-int(x);
	beta=y-int(y);
	gamma=z-int(z);
	slsz=sizy*sizx;

	forceInt.xd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-
	beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*
			(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*alpha*beta*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*
			beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*(alpha*
			beta*gamma);

	forceInt.yd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-
	beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*
			(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*
			beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*
			beta*gamma;

	forceInt.zd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-
	beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*
			(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-
			beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*
			(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*
			beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*
			beta*gamma;

	return(forceInt);
    }


inline void rk2(double x, double y, double z, int sizx, int sizy, int sizz,
double steps, Vector *Force_ini, VoxelPositionDouble *nextPos)
   {
	long slsz;
	Vector OutForce;
	float x1, y1, z1;
	slsz=sizy*sizx;

	OutForce=interpolation(x,y,z,sizx,sizy,sizz,Force_ini);

	x = x + OutForce.xd * steps;
	y = y + OutForce.yd * steps;
	z = z + OutForce.zd * steps;

	nextPos->x = x;
	nextPos->y = y;
	nextPos->z = z;

#ifdef TRACE
	//if(idxSeeds < numBoundSeeds) {

	printf("Next position: (%f, %f, %f) force here: (%f, %f, %f)\n",
		nextPos->x, nextPos->y, nextPos->z,
		OutForce.xd, OutForce.yd, OutForce.zd);

	//}
#endif
   }




//
// function FollowStreamlines
// creates an entire skeleton segment
//
bool FollowStreamlines(
	int originCP, 
	int originSkel,                
	bool lookInCP, 
	double x, double y, double z,
	int sX, int sY, int sZ, Vector* ForceField,
	CriticalPoint *CritPts, int numCritPts,
	int *critPtInSkeleton,
	Skeleton *Skel)
{

  long idx, prividx, slsz;
  int i, streamSteps;
  VoxelPositionDouble Startpos, Nextpos, LastAddedSkelPoint;
  int nrAlreadyInSkel;
  int step;
  bool stop, setfirst;
  int crtSegment, newSeg, segToBeSplit;
  bool endpoint;
  int sp, cp, s;
  int crtSegLength = 0;

  nrAlreadyInSkel = Skel->numPoints; 
     // used to limit the search in the skeleton
     // only to points not in the current segment
  Startpos.x = x;
  Startpos.y = y;
  Startpos.z = z;

  idx = (int)Startpos.z *slsz + (int)Startpos.y *sX + (int)Startpos.x;

  // check whether the start point is close to the already existing skeleton
  if(CloseToSkel( &Startpos, originSkel, Skel, nrAlreadyInSkel, 
		  SP_SMALL_DISTANCE, 0) != -1) 
  {
    // it's close to a skelton point - quit
#ifdef _DEBUG
    printf("Starting point is close to the existing skeleton. Break\n");
#endif
    return true;
  }
  // check if the start point is close to any critical point, except originCP
  // only if lookInCP is true
  if(lookInCP) {
    if(CloseToCP( &Startpos, originCP, CritPts, numCritPts, 
		  SP_SMALL_DISTANCE) != -1) 
      {
	// it's close to a critical point - quit
#ifdef _DEBUG
	printf("Starting point is close to a critical point. Break\n");
#endif
	return true;
      }
  }


  // add the start point to the skeleton
  s = AddSkeletonPoint(Skel, &Startpos);
  
  // add a new segment to the skeleton
  if(originSkel != -1) {
    crtSegment = AddNewSkelSegment(Skel, originSkel, s);
  }
  else {
    crtSegment = AddNewSkelSegment(Skel, s, s);
  }

  // !!
  // there is a new skeleton segment added after this point
  // make sure you set it corretly or delete it before leaving the function
  // !!
  LastAddedSkelPoint = Startpos;

  step = 0;
  stop = false;
  while(!stop)   {
    step++;

    if(step > 8000) {
      // stop if more than 8000 steps were performed
#ifdef TRACE
      printf("Step counter is %d. Break", step);
#endif
      // if the current segment is long enough, keep it, otherwise delete it.
      // however, if the segment originated in a critical point, we should 
      // keep it regardless of the length. Removing it, might disconnect the 
      // skeleton.
      if((originCP != -1) || SegmentIsLongEnough(Skel, crtSegment)) {
	Skel->Segments[crtSegment][SKEL_SEG_RIGHT] = 
	  Skel->Segments[crtSegment][SKEL_SEG_LAST];
      }
      else {
	// current segment is not long enough - delete it
	DeleteLastSkelSegment(Skel);
	// also delete the points added to skeleton during the processing of 
	//   this segment
	Skel->numPoints =  nrAlreadyInSkel;
	// #ifdef _DEBUG
	printf("Current segment is not long enough. Removed.\n");
	// #endif
      }
      return true;
    }
#ifdef TRACE
    printf("Step %d.", step);
#endif
    
    // move to the next position
    rk2( Startpos.x, Startpos.y, Startpos.z, sX, sY, sZ, 0.2,
	 ForceField, &Nextpos);

    // if Nextpos = Startpos then, we should stop
    if(	EQUAL(Nextpos.x, Startpos.x)	&&
	EQUAL(Nextpos.y, Startpos.y)	&&
	EQUAL(Nextpos.z, Startpos.z))
    {
#ifdef TRACE
      printf("New position is the same as start position. Break\n");
#endif
      // if the current segment is long enough, keep it, otherwise delete it.
      // however, if the segment originated in a critical point, we should 
      // keep it regardless of the length. Removing it, might disconnect the 
      // skeleton.
      if((originCP != -1) || SegmentIsLongEnough(Skel, crtSegment)) {
	Skel->Segments[crtSegment][SKEL_SEG_RIGHT] = 
	  Skel->Segments[crtSegment][SKEL_SEG_LAST];
      }
      else {
	// current segment is not long enough - delete it
	DeleteLastSkelSegment(Skel);
	// also delete the points added to skeleton during the processing of 
	//   this segment
	Skel->numPoints =  nrAlreadyInSkel;
	// #ifdef _DEBUG
	printf("Current segment is not long enough. Removed.\n");
	// #endif
      }
      return true;
    }

    // if this point is close to another point that is already part of the
    //   skeleton but not in the same segment,
    // or close to a critical point (if lookInCP is true), 
    //   but not the critical point that started this segment
    // we should stop
    // The points in the same segment are excluded from the search since 
    //   I'm only searching among the points that existed in the Skel array 
    //   when this function started.
    // But it is possible that a single point that belongs to the same 
    //   segment (the first one) was added before the call to this function 
    //   and that point is specified by the originSkel parameter.
    //

   
#ifdef TRACE
    printf("Checking if close to a skeleton point...\n");
    printf("\tcrtSegment = %d\n", crtSegment);
#endif    
    crtSegLength = Skel->Segments[crtSegment][SKEL_SEG_LAST] - 
      Skel->Segments[crtSegment][SKEL_SEG_FIRST];

#ifdef TRACE
    printf("crtSegLength = %d\n", crtSegLength);
#endif 

    sp = CloseToSkel( &Nextpos, originSkel, 
		      Skel, nrAlreadyInSkel, 
		      SP_SMALL_DISTANCE, crtSegLength);
    if(sp != -1) {
#ifdef TRACE
      printf("Position is close to a skeleton point. Break\n");
#endif
      //
      // we are close to a skeleton point, belonging to another segment.
      // if the current segment is not long enough, just delete the current
      // segment
      // however, if the segment originated in a critical point, we should 
      // keep it regardless of the length. Removing it, might disconnect the 
      // skeleton.
      if((originCP != -1) || SegmentIsLongEnough(Skel, crtSegment)) {
	// end the current segment
	Skel->Segments[crtSegment][SKEL_SEG_RIGHT] = sp;

	// find the segment that contains the skeleton point we are close to
	//   (that segment will be split into 2 if we are not close to one
	//    of it's end points)
	endpoint = false;
	segToBeSplit = GetSkelSegment(Skel, sp, &endpoint);
	if(!endpoint) {
	  // not close to one of the endpoints - the segment must be split 
	  // create a new skeleton segment
	  newSeg = AddNewSkelSegment(Skel, sp, sp+1);
	
	  // split the segment
	  Skel->Segments[newSeg][SKEL_SEG_RIGHT] = 
	    Skel->Segments[segToBeSplit][SKEL_SEG_RIGHT];
	  Skel->Segments[newSeg][SKEL_SEG_LAST] = 
	     Skel->Segments[segToBeSplit][SKEL_SEG_LAST];

	  // Skel->Segments[segToBeSplit][SKEL_SEG_LEFT] -> remains unchanged;
	  Skel->Segments[segToBeSplit][SKEL_SEG_RIGHT] = sp;
	  // Skel->Segments[segToBeSplit][SKEL_SEG_FIRST] -> remains unchanged;
	  Skel->Segments[segToBeSplit][SKEL_SEG_LAST] = sp - 1;
	}
	else {
	  // we are close to one of the endpoints of an existing segment
	  // we are happy - nothing to do
	}
      }
      else {
	// current segment is not long enough - delete it
	DeleteLastSkelSegment(Skel);
	// also delete the points added to skeleton during the processing of 
	//   this segment
	Skel->numPoints =  nrAlreadyInSkel;
	// #ifdef _DEBUG
	printf("Current segment is not long enough. Removed.\n");
	// #endif
      }

      return true;
    }
    else {
#ifdef TRACE
      printf("NO.\n Checking if close to a critical point...\n");
#endif     
      if(lookInCP) {
	cp = CloseToCP( &Startpos, originCP, 
			CritPts, numCritPts, 
			SP_SMALL_DISTANCE);
	if(cp != -1) {
#ifdef TRACE
	  printf("Position is close to a critical point. Break\n");
#endif
	  //
	  // we are close to a critical point that is not in the skeleton yet

	  // if the current segment is not long enough, just delete the current
	  // segment
	  // however, if the segment originated in a critical point, we should 
	  // keep it regardless of the length. Removing it, might disconnect 
	  // the skeleton.
	  if((originCP != -1) || SegmentIsLongEnough(Skel, crtSegment)) {

	    // add the critical point to the skeleton
	    //
	    sp = AddSkeletonPoint(Skel, &(CritPts[cp].position));
	    //
	    // mark the critical point as being part of the skeleton
	    //
	    critPtInSkeleton[cp] = sp;

	    // end the current segment
	    Skel->Segments[crtSegment][SKEL_SEG_RIGHT] = sp;
	    
	  }
	  else {
	    // current segment is not long enough - delete it
	    DeleteLastSkelSegment(Skel);
	    // also delete the points added to skeleton during the 
	    //   processing of this segment
	    Skel->numPoints =  nrAlreadyInSkel;
	    // #ifdef _DEBUG
	    printf("Current segment is not long enough. Removed.\n");
	    // #endif
	  }
	  
	  return true;
	}
	else {
#ifdef TRACE
	  printf("NO.\n");
#endif 
	}
      }
    }
    
    

    // if distance from last added position to the next is larger than ...
    //	add the next position to the skeleton.
    // if not, don't add the next position to the skeleton
    if(Distance(&LastAddedSkelPoint, &Nextpos) > SP_SMALL_DISTANCE) {
      // add next position to the skeleton
      s = AddSkeletonPoint(Skel, &Nextpos);

      // set the last point of the current segment to this point
      Skel->Segments[crtSegment][SKEL_SEG_LAST] = s;
      LastAddedSkelPoint = Nextpos;
    }

    // next position becomes current position
    Startpos = Nextpos;
  }

  return true;
}


inline double Distance(	VoxelPositionDouble *pos1,
			VoxelPositionDouble *pos2) 
{
  return ( fabs( pos1->x - pos2->x) + 
	   fabs( pos1->y - pos2->y) + 
	   fabs( pos1->z - pos2->z));
}


int CloseToSkel( VoxelPositionDouble *pos, int originSkel,
		 Skeleton *Skel, int nrAlreadyInSkel, 
		 double maxDistance, int crtSegmentLength)
{
  int i;
  double a, b, c;
  bool endpoint;
  int seg;
  int dToLeft, dToRight;

  // see if it's close to a skeleton point, ignoring the one point specified
  //   by originSkel
  for(i=0; i < nrAlreadyInSkel; i++) {
    // ignore the point specified by originSkel
    if(i == originSkel) continue;
    
    // incremental testing of the "close to" condition for higher speed
    if(fabs(Skel->Points[i].x - pos->x) < maxDistance) {
      a = fabs(Skel->Points[i].x - pos->x);
      if((a + fabs(Skel->Points[i].y - pos->y)) < maxDistance) {
	b = fabs(Skel->Points[i].y - pos->y);
	if((a + b + fabs(Skel->Points[i].z - pos->z)) < maxDistance) {
	  // the point is "close" to a skeleton point
	  // if the skeleton point we are close to is also close to one of 
	  //   the segment end points, (but it's not an end point)
	  //   do not report as being close to the point, wait to get to one
	  //   of the segment end points
	  endpoint = false;
	  seg = GetSkelSegment(Skel, i, &endpoint);
	  if(endpoint) {
	    // report the intersection if we are at an end point
	    return i;
	  }
	  else {
	    //
	    // see how close we are to the endpoints of seg
	    dToLeft  = i - Skel->Segments[seg][SKEL_SEG_FIRST];
	    dToRight = Skel->Segments[seg][SKEL_SEG_LAST] - i;
	    if( (dToLeft  > SP_CLOSE_TO_SEGMENT_ENDPOINT) &&
		(dToRight > SP_CLOSE_TO_SEGMENT_ENDPOINT))
	    {
	      // if we are NOT close to an end point, 
	      //    report the intersection
	      return i;
	    }		
	    // otherwise, we are close to a segment given by seg but, we
	    // are also just a few points from the seg's end point and we 
	    // hope to end the current segment at the same point where seg 
	    // ends. However if the current segment's length is comparable
	    // to SP_CLOSE_TO_SEGMENT_ENDPOINT then we should stop right 
	    // here because this won't look good

	    //  crtSegLength = Skel->Segments[crtSegment][SKEL_SEG_LAST] - 
	    //  Skel->Segments[crtSegment][SKEL_SEG_FIRST];
	    
	    // comparable means not more than twice
	    // SP_CLOSE_TO_SEGMENT_ENDPOINT
	    if(crtSegmentLength <= (2 * SP_CLOSE_TO_SEGMENT_ENDPOINT)) {
	      return i;
	    }
	       
	  }
	}
      }
    }
  }
  
  return -1;
}


int CloseToCP ( VoxelPositionDouble *pos, int	originCP, 
		 CriticalPoint *CritPts, int numCritPts, 
		 double maxDistance)
{
  int i;
  double a, b, c;
	
  // see if it's close to a critical point except the one specified in originCP
  for(i=0; i < numCritPts; i++) {
    // ignore the critical point that started this streamline
    if(i == originCP) continue;

    if(fabs(CritPts[i].position.x - pos->x) < maxDistance) {
      a = fabs(CritPts[i].position.x - pos->x);
      if((a + fabs(CritPts[i].position.y - pos->y)) < maxDistance) {
	b = fabs(CritPts[i].position.y - pos->y);
	if((a + b + fabs(CritPts[i].position.z - pos->z)) < maxDistance) {
	  // the point is "close" to a critical point
	  return i;
	}
      }
    }
  }
  
  return -1;
}

int AddSkeletonPoint(Skeleton *Skel, VoxelPositionDouble *pos) {

#ifdef _DEBUG
  // this shoudn't happen, but just to make sure
  if(Skel == NULL) {
    printf("Skel is NULL !!! Abort\n");
    exit(1);
  }
#endif

  int retval = Skel->numPoints;

  //
  // copy the point's position
  //
  Skel->Points[retval] = (*pos);
  
  // increase the number of skeleton points
  Skel->numPoints = Skel->numPoints + 1;
  if(Skel->numPoints >= Skel->sizePoints ) {
    printf("UPS! - Too many skeleton points to output ! - Abort.\n");
    exit(1);
  }

  return retval;
}

     


//////////////////////////////////////////////////////////////////////////////
// Get level 2 skeleton: adds segments starting at interior high divergence 
//    points
//////////////////////////////////////////////////////////////////////////////
bool GetLevel2Skeleton(
        Vector* ForceField, 		// [in] vector field
	int L, int M, int N,		// [in] vector field size (X, Y and Z)
	VoxelPositionDouble *HDPoints,  // [in] high divergence points array
	int numHDPoints,                // [in] number of high div. points
	Skeleton *Skel                  // [in, out] skeleton
	                                //   in - level 1 skeleton
	                                //   out - level 2 skeleton
) {

#ifdef TRACE
  printf("TRACE: GetLevel2Skeleton function...\n");
#endif
  
  if(Skel == NULL) {
    printf("UPS! - Skel is null!\n");
    exit(1);
  }

  long idx, iidx, slsz, prvidx;
  
  
  int cc;
  int ii,jj,kk;
  int i,j,k;
  slsz = L*M;		// slice size


  //
  // follow the streamlines starting at each high divergence point
  //	until we are close to one of the already existing skeleton points 
  //    ignoring the points in the current segment.
  //

  // start at high divergence points and follow streamlines
  if(HDPoints != NULL) {
    for(i=0; i < numHDPoints; i++) {

#ifdef _DEBUG
      // printf("Starting at high div point %d: (%f, %f, %f)...\n", i, HDPoints[i].x, HDPoints[i].y, HDPoints[i].z);
#endif

      FollowStreamlines(
	  -1, -1, 
	  false, 
	  HDPoints[i].x, HDPoints[i].y, HDPoints[i].z,
	  L, M, N, ForceField,
	  NULL, 0,
	  NULL,
	  Skel);
    }
  }


#ifdef _DEBUG
  PrintElapsedTime("\tSL-2: Following streamlines starting at the high divergence points.");
#endif

  //////
  return true;
}


//
// Allocates the Skel data structure if it's not allocated already
//
bool AllocateSkeleton(Skeleton **Skel, int numPoints, int numSegments) {
  int i;

#ifdef TRACE
  printf("Entering function AllocateSkeleton...\n");
#endif

  if(Skel == NULL) {
    printf("Argument 1 (**Skel) to AllocateSkeleton cannot be NULL ! \
Abort.\n");
    exit(1);
  }

  if((*Skel) == NULL) {
    (*Skel) = new Skeleton;
    if((*Skel) == NULL) {
      printf("UPS! - Error allocating memory for working data structures. \
Abort.\n");
    }
    (*Skel)->sizePoints = 0;
    (*Skel)->Points = NULL;
    (*Skel)->sizeSegments = 0;
    (*Skel)->Segments = NULL;
  }

  (*Skel)->numPoints = 0;
  (*Skel)->numSegments = 0;
  
  if(((*Skel)->Points == NULL) || ((*Skel)->sizePoints == 0)) {
    if(((*Skel)->Points = new VoxelPositionDouble[numPoints]) == NULL) {
      printf("UPS! - Error allocating memory for the working data structure. \
Abort.\n");
      exit(1);
    }
    (*Skel)->sizePoints = numPoints;
  }
  
  if(((*Skel)->Segments == NULL) || ((*Skel)->sizeSegments == 0)) {
    if(((*Skel)->Segments = new (int*)[numSegments]) == NULL) {
      printf("UPS! - Error allocating memory for the working data structure. \
Abort.\n");
      exit(1);
    }
    
    (*Skel)->sizeSegments = numSegments;

    for(i=0; i < (*Skel)->sizeSegments; i++) {
      (*Skel)->Segments[i] = NULL;
    }
    
  }  
  return true;
}


int AddNewSkelSegment(Skeleton *Skel, int left, int first) {
  int retval = -1;

  if(Skel->numSegments >= Skel->sizeSegments) {
    printf("Too many skeleton segments detected (%d). Maximum allowed: %d. \
Abort\n", Skel->numSegments + 1, Skel->sizeSegments);
    exit(1);
  }
  
  if(Skel->Segments[Skel->numSegments] == NULL) {
    if((Skel->Segments[Skel->numSegments] = new int[4]) == NULL) {
      printf("Error allocating memory for the working data structures. \
Abort\n");
      exit(1);
    }
  }
  
  Skel->Segments[Skel->numSegments][SKEL_SEG_LEFT] = left;   
  Skel->Segments[Skel->numSegments][SKEL_SEG_RIGHT] = -1;
  Skel->Segments[Skel->numSegments][SKEL_SEG_FIRST] = first;
  Skel->Segments[Skel->numSegments][SKEL_SEG_LAST] = first; 

  retval = Skel->numSegments;

  Skel->numSegments = Skel->numSegments + 1;
  return retval;
}

bool SegmentIsLongEnough(Skeleton *Skel, int crtSegment) {
  // a segment is long enough if it contains minimum MIN_SEG_LENGTH points
  // between it's left and the right end points
  if(abs(Skel->Segments[crtSegment][SKEL_SEG_LAST] - 
	 Skel->Segments[crtSegment][SKEL_SEG_FIRST]) < MIN_SEG_LENGTH) {
    return false;
  }
  return true;
}


// Function GetSkelSegment
// given a skeleton point, returns the segment it belongs to
//
int  GetSkelSegment(Skeleton *Skel, int point, bool *endpoint) {
  int i;
  for(i=0; i < Skel->numSegments; i++) {
    if( (point == Skel->Segments[i][SKEL_SEG_LEFT]) || 
	(point == Skel->Segments[i][SKEL_SEG_RIGHT]))
    {
      (*endpoint) = true;
      return i;
    }
    
    if( (point >= Skel->Segments[i][SKEL_SEG_FIRST]) && 
	(point <= Skel->Segments[i][SKEL_SEG_LAST]))
    {
      (*endpoint) = false;
      return i;
    }
  }
  
  printf("I shouldn't be here. Abort.\n");
  exit(1);
}

bool DeleteLastSkelSegment(Skeleton *Skel) {
  if(Skel->numSegments > 0) {
    Skel->numSegments = Skel->numSegments - 1;
  }
  return true;
}


//
// Frees a Skeleton data structure
//
bool FreeSkeleton(Skeleton **Skel) {
  int i;

  if(Skel == NULL) {
    return true;
  }

  if((*Skel) == NULL) {
    return true;
  }

  if((*Skel)->Points != NULL) {
    delete [] (*Skel)->Points;
  }

  if((*Skel)->Segments != NULL) {
    for(i=0; i < (*Skel)->sizeSegments; i++) {
      if((*Skel)->Segments[i] != NULL) {
	delete [] (*Skel)->Segments[i];
      }
    }
    delete [] (*Skel)->Segments;
  }

  delete (*Skel);
  (*Skel) = NULL;

  return true;
}
