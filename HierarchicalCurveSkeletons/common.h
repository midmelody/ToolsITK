//////////////////////////////////////////////////////////
// Common include file for the skeletonization project.
//
// Nicu D. Cornea - Wed May 28 16:20:04 EDT 2003
//////////////////////////////////////////////////////////

#ifndef NCD_SKEL_COMMON_DEFINED
#define NCD_SKEL_COMMON_DEFINED

// Includes

#ifdef WIN32
	#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef WIN32
	#include <sys/time.h>
#else
	#include <time.h>
#endif


// Macros

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))

#define EPSILON				0.000001
#define IS_ZERO(nr)			(((nr) < EPSILON) && ((nr) > -EPSILON))
#define EQUAL(nr1, nr2)		(IS_ZERO((nr1) - (nr2)))


// Types

typedef struct {
	short x;
	short y;
	short z;
} VoxelPosition;

typedef struct {
	double xd;   // For large datasets, using double will cause memory shortage
	double yd;
	double zd;
} Vector;

enum CriticalPointType {
	CPT_SADDLE = 1,
	CPT_ATTRACTING_NODE,
	CPT_REPELLING_NODE,
	CPT_UNKNOWN
};

struct  VoxelPositionDouble
{
	double x;
	double y;
	double z;
};

struct CriticalPoint {
	VoxelPositionDouble		position;
	CriticalPointType		type;
        Vector 				evect[3];
	double				eval[3];
};

struct Skeleton {
  VoxelPositionDouble *Points;
  int sizePoints;
  int numPoints;
  int **Segments;
  int sizeSegments;
  int numSegments;
};

/*
struct SkeletonPoint {
  VoxelPositionDouble  position;
  unsigned short       *segments;
  unsigned short       segSize;
  unsigned short       segNumber;
};
*/

// Constants

#define SURF 			100		// surface voxel
#define BOUNDARY		110		// boundary voxel - participates in potential field calculation
#define INTERIOR 		200		// interior voxel
#define PADDING_MIN		210		// added voxels in order to thick the object
#define NR_PAD_VALUES	 40		// are in this range: PADDING_MIN to PADDING_MIN + NR_PAD_VALUES
#define EXTERIOR		  0		// background (exterior to the object) voxel (air)

// some constants
#define SKEL_SEG_LEFT   0   // left end point of the segment
#define SKEL_SEG_RIGHT  1   // right end point of the se
#define SKEL_SEG_FIRST  2   // first point of the segment 
                            //    excluding the left end point
#define SKEL_SEG_LAST   3   // last point of the segment, 
                            //   excluding the right end point

// Functions

void SetStartTime();
void PrintElapsedTime(const char* message);

///////////////////////////////////////////////////////////////////////////////
// reads volume data from a given file
///////////////////////////////////////////////////////////////////////////////
bool ReadVolume(char *filename, int L, int M, int N, unsigned char **vol);

///////////////////////////////////////////////////////////////////////////////
// checks that the volume has 1 plane of empty voxels sepparating the object 
//   from the bounding box in each direction
///////////////////////////////////////////////////////////////////////////////
bool CheckVolumePadding(
        unsigned char *vol, 	 // [in] volume to be checked
	int L, int M, int N   // [in] volume size (X, Y and Z). 
);

///////////////////////////////////////////////////////////////////////////////
// for volume vol, sets the flags array to SURF, INTERIOR or EXTERIOR for 
//     every voxel
// Creates the flags array
///////////////////////////////////////////////////////////////////////////////
bool SetFlags(unsigned char *vol, int L, int M, int N, unsigned char **flags);

///////////////////////////////////////////////////////////////////////////////
// for volume vol, sets the voxel values to SURF, INTERIOR or EXTERIOR
// in place - does not make a copy of the volume
///////////////////////////////////////////////////////////////////////////////
bool FlagVolume(unsigned char *vol, int L, int M, int N);


///////////////////////////////////////////////////////////////////////////////
// tests if the line between the 2 points crosses the boundary of the object
//	that is, if the line passes through a "blank" voxel
///////////////////////////////////////////////////////////////////////////////
bool IsLineCrossingBoundary(
	short x1, short y1, short z1,
	short x2, short y2, short z2,
	int sX, int sY, int sZ,
	unsigned char* flags);

///////////////////////////////////////////////////////////////////////////////
// function SaveSkeleton - saves the skeleton to a file
///////////////////////////////////////////////////////////////////////////////
bool SaveSkeleton(Skeleton *Skel, char *file, char mode = 0);


///////////////////////////////////////////////////////////////////////////////
// Function ReadVectorField
//   Reads a vector field from a file into a Vector array.
///////////////////////////////////////////////////////////////////////////////
bool ReadVectorField(Vector *field, int L, int M, int N, char *fileName);


///////////////////////////////////////////////////////////////////////////////
// Function SaveVectorField
//   Saves a vector field to a file
///////////////////////////////////////////////////////////////////////////////
bool SaveVectorField(Vector *field, int L, int M, int N, char *fileName);

#endif // NCD_SKEL_COMMON_DEFINED
