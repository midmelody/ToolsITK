// Form streamlines from critical points and seed points
//
// --- Author: Xiaosong Yuan, Bala, Nicu D. Cornea - Vizlab, Rutgers University
// --- Date: 9/28/2002
// --- Original file name: skel_streamline_float.cpp
//
// --- Changed to StreamLn.cpp by Nicu D. Cornea, on Tuesday, July 22, 2003
//
//


// #define TRACE

#include "StreamLn.h"

#define MAX_NUM_SKEL	       80000
#define MAX_NUM_SKEL_SEGMENTS  1000

#define STEP_SIZE               0.2   // this is the size of the step
   // used to advance to the next position when following the vector field
   // 0.2 seems to be the best value.

#define PROGRESS_CHECK_INTERVAL  1000   // at each ... th step, we check 
   // whether we made same progress in the last ... steps (some more points
   // were added to the skeleton. If yes, we go on for another ... steps.
   // If not, we stop.
#define MAX_NUM_STEP_INTERVALS   50    // once we made 
// MAX_NUM_STEP_INTERVALS * PROGRESS_CHECK_INTERVAL, we should stop
// following the vector field. We are probably stuck near an attacting node
//

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


typedef enum{
  FSL_EXS_ERROR = 0,
  FSL_EXS_TOO_MANY_STEPS,
  FSL_EXS_CLOSE_TO_SKEL,
  FSL_EXS_CLOSE_TO_CP,
  FSL_EXS_STUCK, 
  FSL_EXS_OUTSIDE
} 
FSL_ExitStatus;

int AddNewSkelSegment(Skeleton *Skel, int point);
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
	bool lookInCP, 
	VoxelPositionDouble *startPt,
	Vector* whereTo,	
	int sX, int sY, int sZ, 
	Vector* ForceField,
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

  Vector whereTo;

  int *critPtInSkeleton;

  //
  // critPtInSkeleton indicates if a critical point is already part of the 
  //   skeleton (is in the Skel->Points array). If -1 then, the point is not 
  //   yet in the Points array. Otherwise it indicates the point's position in 
  //   the skeleton.
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

    printf("L1: Processed %d critical points out of %d.\r", i, numCritPts);
    fflush(stdout);
    
#ifdef TRACE
    idx = (int)CritPts[i].position.z * slsz + (int)CritPts[i].position.y *L
      + (int)CritPts[i].position.x;
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
	
	// direction given by the eigenvector
	whereTo.xd = CritPts[i].evect[j].xd;
	whereTo.yd = CritPts[i].evect[j].yd;
	whereTo.zd = CritPts[i].evect[j].zd;
	FollowStreamlines(
	  i, 
	  true, 
	  &(CritPts[i].position),
	  &whereTo,
	  L, M, N, ForceField,
	  CritPts, numCritPts,
	  critPtInSkeleton,
	  Skel);

	// the exact opposite direction of the eigenvector
	whereTo.xd = - CritPts[i].evect[j].xd;
	whereTo.yd = - CritPts[i].evect[j].yd;
	whereTo.zd = - CritPts[i].evect[j].zd;
	
	FollowStreamlines(
	  i,
	  true,
	  &(CritPts[i].position),
	  &whereTo,
	  L, M, N, ForceField,
	  CritPts, numCritPts,
	  critPtInSkeleton,
	  Skel);

      }
    }

  }

  printf("Done.\n");
#ifdef _DEBUG
  PrintElapsedTime("\tSL-2: Following streamlines starting at the critical points.");
#endif

  delete [] critPtInSkeleton;

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

	if((x > sizx) || (x < 0) || 
	   (y > sizy) || (y < 0) || 
	   (z > sizz) || (z < 0)) 
	{
	  forceInt.xd = 0.00;
	  forceInt.yd = 0.00;
	  forceInt.zd = 0.00;
	  return forceInt;
	}
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
     // long slsz;
     Vector OutForce;
     // float x1, y1, z1;
     // slsz=sizy*sizx;
     double len;

	OutForce=interpolation(x,y,z,sizx,sizy,sizz,Force_ini);

	//
	// normalize
	//
	len = sqrt((OutForce.xd * OutForce.xd) + 
		   (OutForce.yd * OutForce.yd) + 
		   (OutForce.zd * OutForce.zd));

	if(len > 0.00) {
	  OutForce.xd = OutForce.xd / len;
	  OutForce.yd = OutForce.yd / len;
	  OutForce.zd = OutForce.zd / len;
	}

	x = x + OutForce.xd * steps;
	y = y + OutForce.yd * steps;
	z = z + OutForce.zd * steps;

	nextPos->x = x;
	nextPos->y = y;
	nextPos->z = z;

#ifdef TRACE
	printf("Next position: (%f, %f, %f) force here: (%f, %f, %f)\n",
		nextPos->x, nextPos->y, nextPos->z,
		OutForce.xd, OutForce.yd, OutForce.zd);

#endif
   }




//
// function FollowStreamlines
// creates an entire skeleton segment
//
bool FollowStreamlines(
	int originCP, 
	bool lookInCP, 
	VoxelPositionDouble *startPt,
	Vector* whereTo,	
	int sX, int sY, int sZ, 
	Vector* ForceField,
	CriticalPoint *CritPts, int numCritPts,
	int *critPtInSkeleton,
	Skeleton *Skel)
{

#ifdef TRACE
  printf("Starting FollowStreamlines function\n\
\toriginCP = %d\n", originCP);
#endif

  int i;
  VoxelPositionDouble Startpos, Nextpos, LastAddedSkelPoint;
  int nrAlreadyInSkel;
  int step;
  bool stop;
  int crtSegment, newSeg, segToBeSplit;
  bool endpoint;
  int sp, cp, s, oCPSkelPos;
  int crtSegLength = 0;
  int saveCrtSegLen;

  FSL_ExitStatus exit_status = FSL_EXS_ERROR;
  int exit_data = -1;
  bool first_iteration = true;
  div_t res;

  nrAlreadyInSkel = Skel->numPoints; 
     // used to limit the search in the skeleton
     // only to points not in the current segment

  Startpos.x = startPt->x;
  Startpos.y = startPt->y;
  Startpos.z = startPt->z;
  
  oCPSkelPos = -1;  // the position of the original critical point
                    // in the Skeleton array if it was already there
                    // at the start of this function

  //
  // add the start point to the skeleton if not already there
  // the point can already be in the skeleton if it's a critical point
  // that was added as part of a previous segment (example: a saddle point
  // belongs to at least 2 segments)
  // the critPtInSkeleton array tells us if a critical point is already in the
  // skeleton or not and if it is, it tells us it's location.
  //
  
  s = -1;
  if(originCP != -1) {
    s = critPtInSkeleton[originCP];
    oCPSkelPos = s;
    if (s == -1) {
      s = AddSkeletonPoint(Skel, &Startpos);  
      // update the critPtInSkeleton array 
      critPtInSkeleton[originCP] = s;
    }
  }
  else {
    s = AddSkeletonPoint(Skel, &Startpos);
  }
  LastAddedSkelPoint = Startpos;
  
  // add a new segment to the skeleton starting and ending at this point
  crtSegment = AddNewSkelSegment(Skel, s);
  
  // !!
  // there is a new skeleton segment added after this point
  // make sure you set it correctly or delete it before leaving this function
  // !!
 

  exit_status = FSL_EXS_ERROR;
  exit_data = 0;

  first_iteration = true;
  step = 0;
  stop = false;
  saveCrtSegLen = crtSegLength;
  

  while(!stop)   {

#ifdef TRACE
    printf("Current position is: (%f, %f, %f)\n",
	   Startpos.x, Startpos.y, Startpos.z);
#endif

    // if Startpos is on the bounding box or outside, remove this segment
    if( (Startpos.x <= 0) || (Startpos.x >= (sX-1)) || 
	(Startpos.y <= 0) || (Startpos.y >= (sY-1)) || 
	(Startpos.z <= 0) || (Startpos.z >= (sZ-1))) 
    {
      exit_status = FSL_EXS_OUTSIDE;
      stop = true;
      break;
    }

    step++;
    // at every multiple of PROGRESS_CHECK_INTERVAL steps, 
    // take a moment to see if we are moving at all.
    res = div(step, PROGRESS_CHECK_INTERVAL);
    if(res.rem == 0) {
#ifdef TRACE
      printf("Step counter is %d(%d times %d). Check progress...\n", 
	     step, res.quot, PROGRESS_CHECK_INTERVAL);
#endif
      if(crtSegLength > saveCrtSegLen) {
	//
	// yes we are moving - continue
	//
#ifdef TRACE
      printf("Progress: added %d points in the last %d steps.\n", 
	     crtSegLength - saveCrtSegLen, PROGRESS_CHECK_INTERVAL);
#endif
	saveCrtSegLen = crtSegLength;
	//
	// but not forever. 
	// MAX_NUM_STEP_INTERVALS times PROGRESS_CHECK_INTERVAL steps 
	//   should be enough 
	//
	if(res.quot > MAX_NUM_STEP_INTERVALS) {	  
	  // stop if more than ... steps were performed. Most likely we are 
	  // chasing our own tail.
#ifdef TRACE
	  printf("Step counter is %d. Break", step);
#endif
	  exit_status = FSL_EXS_TOO_MANY_STEPS;
	  stop = true;
	  break;
	}
      }
      else {
	// we didn't add any skeleton point in the last 
	//   PROGRESS_CHECK_INTERVAL steps
	// we should stop;
#ifdef TRACE
	printf("No progress. Step counter is %d. Break\n", step);
#endif
	exit_status = FSL_EXS_TOO_MANY_STEPS;
	stop = true;
	break;
      }
    }

#ifdef TRACE
    printf("Step %d.", step);
#endif

    //
    // compute current segment length so far
    //
    crtSegLength = Skel->Segments[crtSegment][SKEL_SEG_LAST] - 
      Skel->Segments[crtSegment][SKEL_SEG_FIRST] + 1;

    // if this point is close to another point that is already part of the
    //   skeleton but not in the same segment,
    // we should stop
    // The points in the same segment are excluded from the search since 
    //   I'm only searching among the points that existed in the Skel array 
    //   when this function started.
   
#ifdef TRACE
    printf("\tcrtSegment = %d\n", crtSegment);
    printf("crtSegLength = %d\n", crtSegLength);

    printf("Checking if close to a skeleton point...\n");    
#endif 

    // 
    // if we are starting from a critical point and we are at the first 
    //   iteration, don't check if we are close to a skeleton point
    // We may be close to a skeleton point now, but we may be moving away from 
    //   it at the next step.
    //
    if((originCP != -1) && (first_iteration)) {
      // don't check 
#ifdef TRACE
	printf("First iteration starting from a critical point. Will not check if close to a skeleton point.\n");
#endif      
    }
    else {
      sp = CloseToSkel( &Startpos, oCPSkelPos,
			Skel, nrAlreadyInSkel, 
			SP_SMALL_DISTANCE, crtSegLength);
      if(sp != -1) {

#ifdef TRACE
	printf("Position is close to a skeleton point (index = %d). Break\n", 
	       sp);
#endif
	exit_status = FSL_EXS_CLOSE_TO_SKEL;
	exit_data = sp;
	stop = true;
	break;
      }
    }

    //
    // check if we are close to a critical point (if lookInCP is true)
    //

#ifdef TRACE
    printf("Checking if close to a critical point...\n");
#endif
    // 
    // if we are starting from a critical point and we are at the first 
    //   iteration, don't check if we are close to a critical point
    // We may be close to a critical point now, but we may be moving away from 
    //   it at the next step.
    //
    if((originCP != -1) && (first_iteration)) {
      // don't check 
#ifdef TRACE
	printf("First iteration starting from a critical point. Will not check if close to a critical point.\n");
#endif      
    }
    else {
      
      if(lookInCP) {
	cp = CloseToCP( &Startpos, originCP, 
			CritPts, numCritPts, 
			SP_SMALL_DISTANCE);
	if(cp != -1) {	
#ifdef TRACE
	  printf("Position is close to a critical point (index = %d). Break\n", 
		 cp);
#endif
	  exit_status = FSL_EXS_CLOSE_TO_CP;
	  exit_data = cp;
	  stop = true;
	  break;
	}
      }
    }

    //
    // add Startpos to the skeleton if not too close to last added point
    // for the first iteration, Startpos is equal to LastAddedSkelPint so
    // the starting position is not added twice.
    //
    if(Distance(&LastAddedSkelPoint, &Startpos) > SP_SMALL_DISTANCE) {

#ifdef TRACE
      printf("Current point added to Skel\n");
#endif
      // add current position to the skeleton
      s = AddSkeletonPoint(Skel, &Startpos);
      
      //
      // update current segment
      //
      if(Skel->Segments[crtSegment][SKEL_SEG_LEFT] ==
	 Skel->Segments[crtSegment][SKEL_SEG_FIRST]) 
      {
	// if LEFT and FIRST are the same, make FIRST equal to this new 
	// position
	Skel->Segments[crtSegment][SKEL_SEG_FIRST] = s;
      }
      Skel->Segments[crtSegment][SKEL_SEG_LAST] = s;
      Skel->Segments[crtSegment][SKEL_SEG_RIGHT] = s;
      
      // update LastAddedSkelPoint 
      LastAddedSkelPoint = Startpos;
    }

 
    // move to the next position
    // for the first iteration, the next position is given by the 
    // whereTo vector, if it's not NULL. If it is NULL, then we take
    // the force field value at the start position
    //
#ifdef TRACE
    printf("Moving to the next position...\n");
#endif 
    if(first_iteration) {
      first_iteration = false;
      if(whereTo != NULL) {
	Nextpos.x = Startpos.x + (whereTo->xd * STEP_SIZE);
	Nextpos.y = Startpos.y + (whereTo->yd * STEP_SIZE);
	Nextpos.z = Startpos.z + (whereTo->zd * STEP_SIZE);
      }
      else {
	rk2( Startpos.x, Startpos.y, Startpos.z, sX, sY, sZ, STEP_SIZE,
	     ForceField, &Nextpos);
      }
    }
    else {
      // for the subsequent iterations, we use the vector field value
      // at the current position
      rk2( Startpos.x, Startpos.y, Startpos.z, sX, sY, sZ, STEP_SIZE,
	   ForceField, &Nextpos);
    }


    //
    // if Nextpos = Startpos we are stuck and we should stop
    //
    if(	EQUAL(Nextpos.x, Startpos.x)	&&
	EQUAL(Nextpos.y, Startpos.y)	&&
	EQUAL(Nextpos.z, Startpos.z))
    {
#ifdef TRACE
      printf("New position is the same as start position. Break\n");
#endif
      exit_status = FSL_EXS_STUCK;
      stop = true;
      break;
    }

    //
    // next position becomes current position
    //
    Startpos = Nextpos;
    first_iteration = false;

  } // end while loop

  
  // usually, if the current segment is not long enough 
  // we are not going to keep it.
  // however, if the segment originated in a critical point, we should 
  // keep it regardless of the length. Removing it, might disconnect the 
  // skeleton.

  //
  // check segment length only if not originated at a critical point
  //
  if(originCP == -1) {
    if(!SegmentIsLongEnough(Skel, crtSegment)) {
      // the segment is not long enough - will be removed
      DeleteLastSkelSegment(Skel);
      // also delete the points added to skeleton during the processing of 
      //   this segment
      Skel->numPoints =  nrAlreadyInSkel;
      // #ifdef _DEBUG
      printf("Current segment is not long enough. Removed.\n");
      // #endif
      
      return true;
    }
  }

  //
  // we are going to keep the segment, let;s finish the job
  //
  switch(exit_status) {
  case FSL_EXS_TOO_MANY_STEPS:
    // nothing to do :)
    break;
  case FSL_EXS_CLOSE_TO_SKEL:
    // we are close to a skeleton point, belonging to another segment.
    
    // end the current segment at the intersection point
    Skel->Segments[crtSegment][SKEL_SEG_RIGHT] = exit_data;

    // find the segment that contains the skeleton point we are close to
    //   (that segment will be split into 2 if we are not close to one
    //    of it's end points)
    endpoint = false;
    segToBeSplit = GetSkelSegment(Skel, exit_data, &endpoint);
    if(!endpoint) {
      // not close to one of the endpoints - the segment must be split 
      // create a new skeleton segment starting at the intersection point
      // and ending where the original segment ended.
      newSeg = AddNewSkelSegment(Skel, exit_data);
      
      // LEFT and FIRST are set to exit_data by AddNewSkelSegment
      // Skel->Segments[newSeg][SKEL_SEG_LEFT] = exit_data;
      // Skel->Segments[newSeg][SKEL_SEG_FIRST] = exit_data;
      Skel->Segments[newSeg][SKEL_SEG_RIGHT] = 
	Skel->Segments[segToBeSplit][SKEL_SEG_RIGHT];
      Skel->Segments[newSeg][SKEL_SEG_LAST] = 
	Skel->Segments[segToBeSplit][SKEL_SEG_LAST];

      // the old segment now has to terminate at the intersection point
      // Skel->Segments[segToBeSplit][SKEL_SEG_LEFT] -> remains unchanged;
      // Skel->Segments[segToBeSplit][SKEL_SEG_FIRST] -> remains unchanged;
      Skel->Segments[segToBeSplit][SKEL_SEG_RIGHT] = exit_data;
      Skel->Segments[segToBeSplit][SKEL_SEG_LAST] = exit_data;
    }
    else {
      // we are close to one of the endpoints of an existing segment
      // we are happy - the original segment doesn't have to be split
    }
    break;

  case FSL_EXS_CLOSE_TO_CP:
    // we are close to a critical point that is not in the skeleton yet
    // (we know it's not in the skeleton yet because we first check if
    //  we are close to a skeleton point, and only of we are not close
    //  to any skeleton point we check whether we are close to a critical 
    //  point)
    
    // add the critical point to the skeleton
    sp = AddSkeletonPoint(Skel, &(CritPts[exit_data].position));
    //
    // mark the critical point as being part of the skeleton
    //
    critPtInSkeleton[exit_data] = sp;

    // end the current segment at the critical point
    Skel->Segments[crtSegment][SKEL_SEG_RIGHT] = sp;

    break;
   
  case FSL_EXS_STUCK:
    // we are stuck in a place - nothing to do, just exit
    break;
  case FSL_EXS_OUTSIDE:
    
    // the segment touched the bounding box
    // the segment will be removed
    DeleteLastSkelSegment(Skel);
    // also delete the points added to skeleton during the processing of 
    //   this segment
    Skel->numPoints =  nrAlreadyInSkel;
    // #ifdef _DEBUG
    printf("Current segment touched the bounding box. Removed.\n");
    // #endif
    
    break;
  case FSL_EXS_ERROR:
    printf("I shouldn't be here !!! Abort\n");
    exit(1);
    break;
  }
  
  /////
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
	    
	    // comparable means not more than twice as long as
	    //    SP_CLOSE_TO_SEGMENT_ENDPOINT
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
#ifdef _DEBUG
    // output the points we got so far
    /*
    int i;
    printf("Skeleton points so far:\n");
    for(i=0; i < Skel->numPoints; i++) {
      printf("%lf %lf %lf 1 0.5\n",
	     Skel->Points[i].x, 
	     Skel->Points[i].y, 
	     Skel->Points[i].z);
    }
    */
#endif
    exit(1);
  }

#ifdef TRACE
  printf("New skelton point added (index = %d)\n", retval);
#endif

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
	  -1,
	  false, 
	  &(HDPoints[i]),
	  NULL, 
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
    if(((*Skel)->Segments = new int*[numSegments]) == NULL) {
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


int AddNewSkelSegment(Skeleton *Skel, int point) {
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
  
  Skel->Segments[Skel->numSegments][SKEL_SEG_LEFT] = point;   
  Skel->Segments[Skel->numSegments][SKEL_SEG_RIGHT] = point;
  Skel->Segments[Skel->numSegments][SKEL_SEG_FIRST] = point;
  Skel->Segments[Skel->numSegments][SKEL_SEG_LAST] = point; 

  retval = Skel->numSegments;

#ifdef TRACE
  printf("New skelton segment added (index = %d)\n", retval);
#endif

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
