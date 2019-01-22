#include "../common.h"
#include "tnt_array2d.h"
#include "jama_eig.h"

bool GetCriticalPoints(
	Vector* ForceField,         // [in] vector field
	int L, int M, int N,        // [in] volume size
	unsigned char *flags,       // [in] volume flags
	CriticalPoint **CritPts,    // [out] detected critical points
	int *numCritPts,            // [out] critical points count
	bool inOut = false          // [in] flag specifying if we know the 
                                    //    inside of the object from the outside
                                    //    false - can distinguish the inside
	                            //    true - cannot distinguish the inside
	                            //    DEFAULF: false
);

