//////////////////////////////////////////////////////////////////////////////
// ----  Computes the skeleton of a volume, using the potential field method.
//
// Input:	volume file name
//		size of volume
//		distance from object boundary where to place the charges (>=0)
//		potential field strenght (1 .. 10)
//		output file name
//
// Output:	skeleton points
//
// Last change: by Nicu D. Cornea
//
//////////////////////////////////////////////////////////////////////////////

// version 
#define pfSkel_Version "1.2.1.1"

char lineSkelOut = 0;

#include "pfSkel.h"

// #define TRACE

// callback functions for the pfSkel module
// for interactive mode
bool ChgParams(Skeleton* Skel, pfSkelCommand *cmd, void* other);  
// for non interactive mode
bool NotInteractive(Skeleton* Skel, pfSkelCommand *cmd, void* other);

int main(int argc, char *argv[]) {
  FILE *fskelout;
  int i, j, L,M,N;         // Sizes in x,y,z dimensions

  int distCharges, fieldStrenght, percHDPts;
  
  Skeleton *Skel = NULL;
  bool interactive, vfout, vfin;
  char *vectorfieldinputfile, *vectorfieldoutputfile;
 

  bool badOption = true;

  cbChangeParameters chgParamsFunction = NULL;
  void* chgParamsArg = NULL;


  /*  // test
  pfSkelCommand cmd;
  cmd.newHDP = 3.00;
  char file[50] = "cow.100x46x68";
  ChgParams(Skel, &cmd, (void*)file);
  exit(1);
  // end test
  */

  // print out verion information
  printf("\n%s - version %s\n", argv[0], pfSkel_Version);

  SetStartTime();

  if (argc < 9) {
    printf(" \
Usage: \n%s <volfile> <xs> <ys> <zs> <distCharges> <fieldStrenght> \n\
            <percHDPts> <SkelOutFile> [options].\n\n\
<volfile> - volume data file name\n\
<xs>, <ys>, <zs> - volume data size\n\
<distCharges> - distance from object boundary where to place the charges \n\
                   ( >=0 )\n\
<fieldStrenght> - potential field strenght (1 .. 10)\n\
<percHDPts> - percentage of highest divergence points to use in skeleton \n\
                 construction (0 .. 100).\n\
<SkelOutFile> - output the skeleton points to the specified file\n\
[options] - options. Valid options are:\n\
   -interactive | -int   = interactive mode. Asks for a new high \n\
      divergence percentage.\n\
   -vectorfieldinputfile | -vfin <file>  = read the vector field from the \n\
      specified file instead of computing it. If this option is present \n\
      the <distCharges> and <fieldStrength> parameters will be ignored. \n\
   -vectorfieldoutputfile | -vfout <file> = dump vector field to this  \n\
      file. \n\
   -lineskeleton | -ls = outputs a straight line skeleton instead of a \n\
      curve. This skeleton can be used for animation.\n\
", argv[0]);

/* Options to add next

  -noeuclideandistace | -ned = do not use euclidean distance. uses \n\
      Manhattan distance instead - faster, but less accurate. \n\
      IGNORED\n\
   -halfboundarypoints | -hbp = use only half of the boundary points to \n\
      place charges - faster but less accurate. IGNORED\n\

*/

    exit(1);
  }

  L 		= atoi(argv[2]);
  M 		= atoi(argv[3]);
  N 		= atoi(argv[4]); 
  distCharges 	= atoi(argv[5]);
  fieldStrenght = atoi(argv[6]);
  percHDPts	= atoi(argv[7]);

  vfin = false;
  vfout = false;
  vectorfieldinputfile = NULL;
  vectorfieldoutputfile = NULL;
  interactive = false;
  lineSkelOut = 0;

  if(argc > 9) {
    i = 9;
    while(i < argc) {
      badOption = true;
      if((strcmp(argv[i], "-interactive") == 0) || 
	 (strcmp(argv[i], "-int") == 0)) 
      {
	interactive = true;
	badOption = false;
      }
      else {
	if((strcmp(argv[i], "-vectorfieldinputfile") == 0) || 
	   (strcmp(argv[i], "-vfin") == 0)) 
	{
	   badOption = false;
	  if(argc > i+1) {
	    vfin = true;
	    vectorfieldinputfile = argv[i+1];
	    i++;
	  }
	  else {
	    vfin = false;
	    vectorfieldinputfile = NULL;
	    printf("** Need a file name after %s option. Ignored.\n", argv[i]);
	  }
	}
	else {
	  if((strcmp(argv[i], "-vectorfieldoutputfile") == 0) || 
	   (strcmp(argv[i], "-vfout") == 0)) 
	  {
	    badOption = false;
	    if(argc > i+1) {
	      vfout = true;
	      vectorfieldoutputfile = argv[i+1];
	      i++;
	    }
	    else {
	      vfout = false;
	      vectorfieldoutputfile = NULL;
	      printf("** Need a file name after %s option. Ignored.\n", 
		     argv[i]);
	    }
	  }
	  else {
	     if((strcmp(argv[i], "-lineskeleton") == 0) || 
		(strcmp(argv[i], "-ls") == 0)) 
	     {
	       badOption = false;
	       // output line skeleton
	       lineSkelOut = 1;
	     }
	  }
	}
      }

      if(badOption) {
	printf(" ** Unrecognized option: %s\n", argv[i]);
      }
      
      i++;
    }
  }

  //
  // setup interactive mode 
  //
  if(interactive) {
#ifdef _DEBUG
    printf("Interactive mode.\n");
#endif
    chgParamsFunction = &ChgParams;
    chgParamsArg = (void*)argv[8];
  }
  else {
    // not interactive
    chgParamsFunction = &NotInteractive;
    chgParamsArg = (void*)argv[8];
  }
  
  //
  // Allocate the Skeleton data structure
  //
  Skel = NULL;
  AllocateSkeleton(&Skel, 50000, 5000);

  //
  // compute the skeleton
  //
  pfSkel(argv[1], L, M, N, distCharges, fieldStrenght, percHDPts, 
	 Skel, 
	 chgParamsFunction, chgParamsArg, 
	 vectorfieldinputfile, vectorfieldoutputfile);
  

  //
  // free the skeleton data structure
  //
  FreeSkeleton(&Skel);


  printf("Done.\n");
  PrintElapsedTime("");
  return 0;
}


//
// callback function for the pfSkel module for interactive mode
//   This function is called from inside pfSkel function
//      when the computation is done and it allows changing of the parameters

bool ChgParams(Skeleton *Skel, pfSkelCommand *cmd, void *other) {
  char c[2];
  float newHD, newHC;

#ifdef TRACE
  printf("Start ChgParams\n");
  printf("Skel = %p, cmd = %p, other = %p\n", Skel, cmd, other);
#endif

  if(cmd == NULL) {
    return false;
  }
  if(other == NULL) {
    return false;
  } 

  // first thing: save the current skeleton.
  // cmd->newHDP contains the current value for the high divergence percentage
  // other should be a pointer to the filename
  //
  newHD = cmd->HDP;

#ifdef TRACE
  printf("Base output file name: %s\n", (char*)other);
#endif
  
  char *outfile;
  int len = strlen((char*)other) + 20;

#ifdef TRACE
  printf("Allocating %d chars for new file name.\n", len);
#endif

  if((outfile = new char[len]) == NULL) {
    printf("Error allocating memory for the working data structures. Abort.\n\
");
    return false;
  }
  

  sprintf(outfile, "%s-hd%3.2f.skel", (char*)other, newHD);
  printf("Saving skeleton to file %s ...\n", outfile);
  SaveSkeleton(Skel, outfile, lineSkelOut);

#ifdef TRACE
  printf("Deleting outfile....\n");
#endif
  delete [] outfile;

  printf("done.\n");

  
  //
  // present a menu
  //
  printf("Available commands: \n");
  printf("\
\tq - quit\n\
\tp - change parameters (you will be prompted for new values)\n");
  printf("Command: ");
  
  while(true) {
    scanf("%1c", &c);
    switch(c[0]) {
    case 'q':
      // quit
      cmd->cmdCode[0] = c[0];
      return false;
    case 'p':      
      cmd->cmdCode[0] = c[0];
      printf("New high divergence percentage (current value: %f): ", 
	     cmd->HDP);
      scanf("%f", &newHD);
      cmd->HDP = newHD;
      getchar();
      return true;
    default:
      printf("Invalid command: %c. See available commands above\nCommand: ",c);
    }
  }
  
  return false;
}



//
// callback function for the pfSkel module for non interactive mode
//   This function is called from inside pfSkel function
//      when the computation is done and it allows saving the skeleton

bool NotInteractive(Skeleton *Skel, pfSkelCommand *cmd, void *other) {
  char c[2];
  float newHD, newHC;

#ifdef TRACE
  printf("Start Not Interactive\n");
  printf("Skel = %p, cmd = %p, other = %p\n", Skel, cmd, other);
#endif

  if(cmd == NULL) {
    return false;
  }
  if(other == NULL) {
    return false;
  } 

  // first thing: save the current skeleton.
  // cmd->newHDP contains the current value for the high divergence percentage
  // other should be a pointer to the filename
  //
  newHD = cmd->HDP;

#ifdef TRACE
  printf("Base output file name: %s\n", (char*)other);
#endif
  
  char *outfile;
  int len = strlen((char*)other) + 20;

#ifdef TRACE
  printf("Allocating %d chars for new file name.\n", len);
#endif

  if((outfile = new char[len]) == NULL) {
    printf("Error allocating memory for the working data structures. Abort.\n\
");
    return false;
  }
  

  sprintf(outfile, "%s-hd%3.2f.skel", (char*)other, newHD);
  printf("Saving skeleton to file %s ...\n", outfile);
  SaveSkeleton(Skel, outfile, lineSkelOut);

#ifdef TRACE
  //  printf("Deleting outfile....\n");
#endif
  delete [] outfile;

  printf("done.\n");

  
  //
  // no menu - just give quit command
  //
  cmd->cmdCode[0] = 'q';
  return false;
}

