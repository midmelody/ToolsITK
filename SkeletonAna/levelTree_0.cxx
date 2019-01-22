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


int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile outputTreeGenFile outputSubTreeFile" << std::endl;
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
  FILE *treeNewFile = fopen(argv[2], "w");
  for(a=0;a<branchCt;a++)
    {
      std::list<VoxelType>::iterator voxelIt;
      for(voxelIt=branchList[a].voxelIndex.begin(); voxelIt!=branchList[a].voxelIndex.end(); voxelIt++)
	fprintf(treeNewFile, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, a+1, branchList[a].generation);
    }
  fclose(treeNewFile);
  
  //Write sub Tree
  FILE *treeSubFile = fopen(argv[3], "w");
  for(a=0;a<branchCt;a++)
    {
      std::list<VoxelType>::iterator voxelIt;
      for(voxelIt=branchList[a].voxelIndex.begin(); voxelIt!=branchList[a].voxelIndex.end(); voxelIt++)
	fprintf(treeSubFile, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, a+1, branchList[a].subTree);
    }
  fclose(treeSubFile);

  return EXIT_SUCCESS;
}
