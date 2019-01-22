#include "itkImage.h"
#include "itkRecursiveGaussianImageFilter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

struct VoxelType
{
  float x;
  float y;
  float z;
};

//Branch struct
struct BranchType
{
  int voxelCt;
  int label;
  int parLabel;
  VoxelType * voxelIndex;
};

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile outputTreeFile " << std::endl;
      return EXIT_FAILURE;
    }

  //Read tree
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
  while(fscanf(treeFile, "%f %f %f %d %d", &index.x, &index.y, &index.z, &label, &parLabel) == 5)
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
    branchList[a].voxelCt = 0;
  treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(treeFile, "%f %f %f %d %d", &index.x, &index.y, &index.z, &label, &parLabel) == 5)
    {
      branchList[label-1].voxelCt++;
      branchList[label-1].label = label;
      branchList[label-1].parLabel = parLabel;    
    }
  fclose(treeFile);

  //Record point index
  int curVoxel[branchCt];
  for(a=0;a<branchCt;a++)
    {
      branchList[a].voxelIndex = new VoxelType[branchList[a].voxelCt];
      curVoxel[a] = 0;
    }
  treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(treeFile, "%f %f %f %d %d", &index.x, &index.y, &index.z, &label, &parLabel) == 5)
    {
      branchList[label-1].voxelIndex[curVoxel[label-1]] = index;
      curVoxel[label-1]++;
    }
  fclose(treeFile);

  /*
  //Branch voxel count
  for(a=0;a<branchCt;a++)
    {
      std::cout<<"("<<a+1<<": "<<branchList[a].voxelCt<<" "<<branchList[a].label<<" "<<branchList[a].parLabel<<"); ";
    }
  std::cout<<std::endl;
  */
  //ITK smooth
  const unsigned int Dimension = 1;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >  SmoothingFilterType;

  FILE *outFile = fopen(argv[2], "w+");
  for(a=0;a<branchCt;a++)
    {
      ImageType::SizeType size;
      size[0] = branchList[a].voxelCt;
      float smX[branchList[a].voxelCt]; 
      float smY[branchList[a].voxelCt];
      float smZ[branchList[a].voxelCt];
      int lineIdx;
      if(branchList[a].voxelCt>10)
	{
	  ImageType::RegionType region;
	  region.SetSize( size );
	  ImageType::Pointer lineImageX = ImageType::New();
	  ImageType::Pointer lineImageY = ImageType::New();
	  ImageType::Pointer lineImageZ = ImageType::New();
	  lineImageX->SetRegions( region );
	  lineImageX->Allocate(); 
	  lineImageY->SetRegions( region );
	  lineImageY->Allocate(); 
	  lineImageZ->SetRegions( region );
	  lineImageZ->Allocate(); 
	  ImageType::IndexType pixelIndex;

	  SmoothingFilterType::Pointer smoothingFilterX = SmoothingFilterType::New();
	  SmoothingFilterType::Pointer smoothingFilterY = SmoothingFilterType::New();
	  SmoothingFilterType::Pointer smoothingFilterZ = SmoothingFilterType::New();
	  smoothingFilterX->SetSigma( 3 ); 
	  smoothingFilterY->SetSigma( 3 ); 
	  smoothingFilterZ->SetSigma( 3 ); 
	  
	  for(lineIdx = 0; lineIdx<branchList[a].voxelCt; lineIdx++)
	    {
	      pixelIndex[0] = lineIdx;
	      lineImageX->SetPixel(pixelIndex, branchList[a].voxelIndex[lineIdx].x);
	      lineImageY->SetPixel(pixelIndex, branchList[a].voxelIndex[lineIdx].y);
	      lineImageZ->SetPixel(pixelIndex, branchList[a].voxelIndex[lineIdx].z);
	    }
	  smoothingFilterX->SetInput( lineImageX );
	  smoothingFilterX->Update();
	  smoothingFilterY->SetInput( lineImageY );
	  smoothingFilterY->Update();
	  smoothingFilterZ->SetInput( lineImageZ );
	  smoothingFilterZ->Update();

	  for(lineIdx = 0; lineIdx<branchList[a].voxelCt; lineIdx++)
	    {
	      pixelIndex[0] = lineIdx;
	      smX[lineIdx] = smoothingFilterX->GetOutput()->GetPixel(pixelIndex);
	      smY[lineIdx] = smoothingFilterY->GetOutput()->GetPixel(pixelIndex);
	      smZ[lineIdx] = smoothingFilterZ->GetOutput()->GetPixel(pixelIndex);
	    }   
	}
      else
	{
	  for(lineIdx = 0; lineIdx<branchList[a].voxelCt; lineIdx++)
	    {
	      smX[lineIdx] = branchList[a].voxelIndex[lineIdx].x;
	      smY[lineIdx] = branchList[a].voxelIndex[lineIdx].y;
	      smZ[lineIdx] = branchList[a].voxelIndex[lineIdx].z;
	    }   
	}

      //write out
      for(lineIdx = 0; lineIdx<branchList[a].voxelCt; lineIdx++)
	fprintf(outFile, "%f %f %f %d %d\n", smX[lineIdx], smY[lineIdx], smZ[lineIdx], branchList[a].label, branchList[a].parLabel);
      
      //std::cout<<"("<<a+1<<": "<<branchList[a].voxelCt<<" "<<branchList[a].label<<" "<<branchList[a].parLabel<<"); ";
    }
  //std::cout<<std::endl;

  fclose(outFile);

  return EXIT_SUCCESS;
}
