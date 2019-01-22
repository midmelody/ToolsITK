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
};


int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputImageFile outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }
 

  //Read skeleton
  VoxelType index;
  int label, parLabel;
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
      //Record parent child relationship
      int par = branchList[a].parentLabel;
      if(par>=1)
	{
	  branchList[par-1].childCt++;
	  branchList[par-1].childLabel.push_back(a+1);
	}
      //Initialize maxLength as voxel count
      branchList[a].maxLengthFromEnd = branchList[a].voxelCt;
    }
 
  //Dynamic programming for tree reconstruction
  for(a=branchCt-1;a>0;a--)
    {
      //Last branch layer
      if(branchList[a].childCt==0)
	{
	  branchList[a].maxLengthChildLabel = 0;
	}
      //update parent
      int par = branchList[a].parentLabel - 1;
      if(par>=0)
	{
	  int curLength = branchList[par].maxLengthFromEnd - branchList[par].voxelCt;
	  int childLength = branchList[a].maxLengthFromEnd;
	  if(curLength<childLength)
	    {
	      branchList[par].maxLengthFromEnd = childLength + branchList[par].voxelCt;
	      branchList[par].maxLengthChildLabel = a+1;
	    }
	}
    }

  //Print current tree info
  /*
  for(a=0;a<branchCt;a++)
    {
      std::cout<<"Branch #"<<a+1<<":"<<std::endl;
      std::cout<<"Parent: "<<branchList[a].parentLabel<<std::endl;
      std::cout<<"Total Voxel: "<<branchList[a].voxelCt<<std::endl;
      std::cout<<"Total Children: "<<branchList[a].childCt<<std::endl;
      std::cout<<"MaxLengthFromEnd: "<<branchList[a].maxLengthFromEnd<<std::endl;
      std::cout<<"MaxLengthChildLabel: "<<branchList[a].maxLengthChildLabel<<std::endl;
      std::cout<<"Children List:"<<std::endl;
      std::list<int>::iterator childrenIt;
      for(childrenIt=branchList[a].childLabel.begin(); childrenIt!=branchList[a].childLabel.end(); childrenIt++)
	std::cout<<(*childrenIt)<<std::endl;
      std::cout<<"-------------------------------------------------------"<<std::endl;
    }
  */

  //Disable branch with small absolute length or relative length (length from end point) 
  bool disableFlag[branchCt];
  int lengthThreshold = 5;
  float ratioThreshold = 0.2;
  for(a=0;a<branchCt;a++)
    {
      disableFlag[a] = false;
      int curLength = branchList[a].maxLengthFromEnd;
      //Absolute length
      if(curLength<=lengthThreshold)
	disableFlag[a] = true;
      //Relative length
      int par = branchList[a].parentLabel - 1; 
      if(par>=0)
	{
	  float minRatio = 100;
	  std::list<int>::iterator childrenIt;
	  for(childrenIt=branchList[par].childLabel.begin(); childrenIt!=branchList[par].childLabel.end(); childrenIt++)
	    {
	      int bro = (*childrenIt) - 1;
	      int broLength = branchList[bro].maxLengthFromEnd;
	      float ratio = float(curLength)/float(broLength);
	      if(minRatio>ratio)
		minRatio = ratio;
	    }
	  if(minRatio<ratioThreshold)
	    disableFlag[a] = true;
	}
    }
  //If disabled, disable whold tree
  for(a=0;a<branchCt;a++)
    {
      if(disableFlag[a])
	{
	  std::list<int>::iterator childrenIt;
	  for(childrenIt=branchList[a].childLabel.begin(); childrenIt!=branchList[a].childLabel.end(); childrenIt++)
	    {
	      int childIndex = (*childrenIt) - 1;
	      disableFlag[childIndex] = true;
	    }
	}
    }

  //Output new tree
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  //Parameters
  reader->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
  //Pipeline
  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate image for output
  ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions( region );
  outImage->SetSpacing( spacing );
  outImage->SetOrigin( origin );
  outImage->SetDirection( direction );
  outImage->Allocate();  
  //Initialize
  ImageType::IndexType pixelIndex;
  int i, j, k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  outImage->SetPixel(pixelIndex, 0);
	}
  int b;
  for(a=0;a<branchCt;a++)
    {
      if(!disableFlag[a])
	{
	  for(b=0;b<branchList[a].voxelCt;b++)
	    {
	      pixelIndex[0] = branchList[a].voxelIndex.front().x;
	      pixelIndex[1] = branchList[a].voxelIndex.front().y;
	      pixelIndex[2] = branchList[a].voxelIndex.front().z;
	      branchList[a].voxelIndex.pop_front();
	      outImage->SetPixel(pixelIndex,1);
	    }
	}
    }

  //Write output  
  writer->SetInput( outImage );
  try
    {
      writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }

  /*
  //Relabel: Prune and Merge
  int newLabel[branchCt];
  bool visitedFlag[branchCt];
  for(a=0;a<branchCt;a++)
    {
      newLabel[a] = 0;
      visitedFlag[a] = false;
    }
  int curReLabel = 1;
  for(a=0;a<branchCt;a++)
    {
      //not visited and not marked
      if((!visitedFlag[a])&&(newLabel[a]==0))
	{
	  //if disabled, label map to 0
	  if(disableFlag[a])
	    newLabel[a] = 0;
	  else
	    {
	      newLabel[a] = curReLabel;
	      //check children count
	      int childrenCt = 0;
	      std::list<int>::iterator childrenIt;
	      for(childrenIt=branchList[a].childLabel.begin(); childrenIt!=branchList[a].childLabel.end(); childrenIt++)
		{
		  int childIndex = (*childrenIt) - 1;
		  //count the number of not disabled
		  if(!disableFlag[childIndex])
		    childrenCt++;
		}
	      //if more than 1 child or leaf node, update current label
	      if(childrenCt>1)
		curReLabel++;
	      if(childrenCt==0)
		curReLabel++;
	      //if only 1 child, merge
	      if(childrenCt==1)
		{
		  for(childrenIt=branchList[a].childLabel.begin(); childrenIt!=branchList[a].childLabel.end(); childrenIt++)
		    {
		      int childIndex = (*childrenIt) - 1;
		      if(!disableFlag[childIndex])
			{
			  newLabel[childIndex] = newLabel[a];
			  branchList[childIndex].parentLabel = branchList[a].parentLabel;
			  //visitedFlag[childIndex] = true;
			}
		    }	 
		}
	    }
	}
    }

  //Print the info to check
  for(a=0;a<branchCt;a++)
    {
      std::cout<<"Branch #"<<a+1<<":"<<std::endl;
      std::cout<<"Parent: "<<branchList[a].parentLabel<<std::endl;
      std::cout<<"Total Voxel: "<<branchList[a].voxelCt<<std::endl;
      std::cout<<"Total Children: "<<branchList[a].childCt<<std::endl;
      std::cout<<"MaxLengthFromEnd: "<<branchList[a].maxLengthFromEnd<<std::endl;
      std::cout<<"MaxLengthChildLabel: "<<branchList[a].maxLengthChildLabel<<std::endl;
      std::cout<<"DisabledFlag: "<<disableFlag[a]<<std::endl;
      std::cout<<"NewLabel: "<<newLabel[a]<<std::endl;
      std::cout<<"Voxel List:"<<std::endl;
      std::list<VoxelType>::iterator voxelIt;
      for(voxelIt=branchList[a].voxelIndex.begin(); voxelIt!=branchList[a].voxelIndex.end(); voxelIt++)
	std::cout<<(*voxelIt).x<<" "<<(*voxelIt).y<<" "<<(*voxelIt).z<<std::endl;
      std::cout<<"Children List:"<<std::endl;
      std::list<int>::iterator childrenIt;
      for(childrenIt=branchList[a].childLabel.begin(); childrenIt!=branchList[a].childLabel.end(); childrenIt++)
	std::cout<<(*childrenIt)<<std::endl;
      std::cout<<"-------------------------------------------------------"<<std::endl;
    }
  
  //Write Tree
  FILE *treeNewFile = fopen(argv[2], "w");
  for(a=0;a<branchCt;a++)
    {
      if(newLabel[a]!=0)
	{
	  std::list<VoxelType>::iterator voxelIt;
	  for(voxelIt=branchList[a].voxelIndex.begin(); voxelIt!=branchList[a].voxelIndex.end(); voxelIt++)
	    {
	      if(branchList[a].parentLabel!=0)
		fprintf(treeNewFile, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, newLabel[a], newLabel[branchList[a].parentLabel-1]);
	      else
		fprintf(treeNewFile, "%d %d %d %d %d\n", (*voxelIt).x, (*voxelIt).y, (*voxelIt).z, newLabel[a], 0);
	    }
	}
    }
  fclose(treeNewFile);
   */

  return EXIT_SUCCESS;
}
