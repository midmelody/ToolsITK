#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <list>
#include "assert.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include <queue>

//ITK settings
const unsigned int Dimension = 3;
typedef unsigned char PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

//Voxel struct
struct VoxelType
{
  int x;
  int y;
  int z;
};

//Branch struct
struct BranchType
{
  int parentLabel;
  int voxelCt;
  std::list<VoxelType> voxelIndex;
  int maxLengthFromEnd;
  int maxLengthChildLabel;
  int childCt;
  std::list<int> childLabel;
  int generation;
  int subTree;
};

//Remove extension
char *remove_ext (char* mystr, char dot, char sep) 
{
  char *retstr, *lastdot, *lastsep;
  // Error checks and allocate string.
  if (mystr == NULL)
    return NULL;
  if ((retstr = (char*) malloc (strlen (mystr) + 1)) == NULL)
    return NULL;
  // Make a copy and find the relevant characters.
  strcpy (retstr, mystr);
  lastdot = strrchr (retstr, dot);
  lastsep = (sep == 0) ? NULL : strrchr (retstr, sep);
  // If it has an extension separator.
  if (lastdot != NULL) 
    {
      // and it's before the extenstion separator.
      if (lastsep != NULL) 
	{
	  if (lastsep < lastdot) 
	    {
	      // then remove it.
	      *lastdot = '\0';
            }
        } 
      else 
	{
	  // Has extension separator with no path separator.
	  *lastdot = '\0';
        }
    }
  // Return the modified string.
  return retstr;
}

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile outputTreeGenFilePre mode(0: whole tree; 1: two sub-trees)" << std::endl;
      return EXIT_FAILURE;
    }
 
  //Read skeleton
  VoxelType index;
  int label, parLabel, generation;
  //find total branch number
  FILE *treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  int branchCt = 0;
  int curLabel = 0;
  while(fscanf(treeFile, "%d %d %d %d %d", &index.x, &index.y, &index.z, &label, &parLabel) == 5)
    {
      if(label > curLabel)
	{
	  branchCt++;
	  curLabel = label;
	}
    }
  fclose(treeFile);

  //get information
  //index in branchList: branchLabel-1
  BranchType * branchList = new BranchType[branchCt];
  int a;
  for(a=0;a<branchCt;a++)
    {
      branchList[a].parentLabel = 0;
      branchList[a].voxelCt = 0;
      branchList[a].maxLengthFromEnd = 0;
      branchList[a].maxLengthChildLabel = 0;
      branchList[a].childCt = 0;
      branchList[a].generation = 1;
      if(a==0)
	branchList[a].subTree = 1;
      else if(a==1)
	branchList[a].subTree = 2;
      else if(a==2)
	branchList[a].subTree = 3;
      else
	branchList[a].subTree = 0;	
    }
 
  treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(treeFile, "%d %d %d %d %d", &index.x, &index.y, &index.z, &label, &parLabel) == 5)
    {
      branchList[label-1].parentLabel = parLabel;
      branchList[label-1].voxelCt++;
      branchList[label-1].voxelIndex.push_back(index);
    }
  fclose(treeFile);
 
  for(a=0;a<branchCt;a++)
    {
      int par = branchList[a].parentLabel-1;
      if(par>=0)
	{
	  branchList[par].childCt++;
	  branchList[par].childLabel.push_back(a+1);
	  branchList[a].generation = branchList[par].generation+1;
	  if(a>2)
	    branchList[a].subTree = branchList[par].subTree;
	}
    }
 
  //Write generation Tree
  int outMode = atoi(argv[3]);
  char outFileName[100]; 
  char outFileName1[100]; 
  char outFileName2[100]; 
  char outFileName3[100]; 
  char outFileName4[100];
  if(outMode==0)
    {
      strcpy(outFileName,argv[2]);
      strcat(outFileName,".txt");
      FILE *treeNewFile = fopen(outFileName, "w");
      for(a=0;a<branchCt;a++)
	{
	  std::list<VoxelType>::iterator voxelIt;
	  for(voxelIt=branchList[a].voxelIndex.begin(); voxelIt!=branchList[a].voxelIndex.end(); voxelIt++)
	    fprintf(treeNewFile, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, a+1, branchList[a].generation);
	}
      fclose(treeNewFile);
    }
  else
    {
      strcpy(outFileName1,argv[2]);
      strcat(outFileName1,"_1.txt");
      FILE *treeNewFile1 = fopen(outFileName1, "w");
      strcpy(outFileName2,argv[2]);
      strcat(outFileName2,"_2.txt");
      FILE *treeNewFile2 = fopen(outFileName2, "w");

      strcpy(outFileName3,remove_ext(argv[1], '.', '/'));
      strcat(outFileName3,"_1.txt");
      FILE *treeOriNewFile1 = fopen(outFileName3, "w");
      strcpy(outFileName4,remove_ext(argv[1], '.', '/'));
      strcat(outFileName4,"_2.txt");
      FILE *treeOriNewFile2 = fopen(outFileName4, "w");

      for(a=0;a<branchCt;a++)
	{
	  std::list<VoxelType>::iterator voxelIt;
	  for(voxelIt=branchList[a].voxelIndex.begin(); voxelIt!=branchList[a].voxelIndex.end(); voxelIt++)
	    {
	      if(branchList[a].subTree==2)
		{
		  fprintf(treeNewFile1, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, a+1, branchList[a].generation-1);
		  fprintf(treeOriNewFile1, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, a+1, branchList[a].parentLabel);
		}
	      if(branchList[a].subTree==3)
		{
		  fprintf(treeNewFile2, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, a+1, branchList[a].generation-1);	
		  fprintf(treeOriNewFile2, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, a+1, branchList[a].parentLabel);	
		}      
	    }
	}
      fclose(treeNewFile1); 
      fclose(treeNewFile2); 
      fclose(treeOriNewFile1); 
      fclose(treeOriNewFile2); 

      /*
      //Write sub Tree
      //FILE *treeSubFile = fopen(argv[3], "w");
      //for(a=0;a<branchCt;a++)
	{
	  std::list<VoxelType>::iterator voxelIt;
	  for(voxelIt=branchList[a].voxelIndex.begin(); voxelIt!=branchList[a].voxelIndex.end(); voxelIt++)
	    fprintf(treeSubFile, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, a+1, branchList[a].subTree);
	}
      fclose(treeSubFile);
      */
    }
  return EXIT_SUCCESS;
}
