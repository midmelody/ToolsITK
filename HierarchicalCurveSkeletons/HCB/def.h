// type definitions
#ifndef NCD_BOUND_DEF
#define NCD_BOUND_DEF

#include "../common.h"

typedef struct {
	VoxelPosition	point;
	unsigned int* 	neighbors;
	unsigned char 	nrOfNeighbors;
	Vector			normal;
	double			curvature;
	double 			curvk2;
} bdElement;

#endif // NCD_BOUND_DEF
