#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

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
  VoxelType * voxelIndex;
  VoxelType * voxelGradient;
};

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputImageFile outputImageFile " << std::endl;
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
    }
  fclose(treeFile);

  //Record point index
  int curVoxel[branchCt];
  for(a=0;a<branchCt;a++)
    {
      branchList[a].voxelIndex = new VoxelType[branchList[a].voxelCt];
      branchList[a].voxelGradient = new VoxelType[branchList[a].voxelCt];
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
      branchList[label-1].voxelGradient[curVoxel[label-1]].x = 0;
      branchList[label-1].voxelGradient[curVoxel[label-1]].y = 0;
      branchList[label-1].voxelGradient[curVoxel[label-1]].z = 0;
      curVoxel[label-1]++;
    }
  fclose(treeFile);

  //Compute tangent at every location using B-spline model
  float BSplineCoe[4];
  float t = 0;
  int b;
  BSplineCoe[0] = (-3+6*t-3*t*t)/6;
  BSplineCoe[1] = (-12*t+9*t*t)/6;
  BSplineCoe[2] = (3+6*t-9*t*t)/6;
  BSplineCoe[3] = (3*t*t)/6;
  for(a=0;a<branchCt;a++)
    {
      for(b=1;b<branchList[a].voxelCt-2;b++)
	{
	  branchList[a].voxelGradient[b].x = branchList[a].voxelGradient[b].x+BSplineCoe[0]*branchList[a].voxelIndex[b-1].x;
	  branchList[a].voxelGradient[b].x = branchList[a].voxelGradient[b].x+BSplineCoe[1]*branchList[a].voxelIndex[b].x;
	  branchList[a].voxelGradient[b].x = branchList[a].voxelGradient[b].x+BSplineCoe[2]*branchList[a].voxelIndex[b+1].x;
	  branchList[a].voxelGradient[b].x = branchList[a].voxelGradient[b].x+BSplineCoe[3]*branchList[a].voxelIndex[b+2].x;
	  branchList[a].voxelGradient[b].y = branchList[a].voxelGradient[b].y+BSplineCoe[0]*branchList[a].voxelIndex[b-1].y;
	  branchList[a].voxelGradient[b].y = branchList[a].voxelGradient[b].y+BSplineCoe[1]*branchList[a].voxelIndex[b].y;
	  branchList[a].voxelGradient[b].y = branchList[a].voxelGradient[b].y+BSplineCoe[2]*branchList[a].voxelIndex[b+1].y;
	  branchList[a].voxelGradient[b].y = branchList[a].voxelGradient[b].y+BSplineCoe[3]*branchList[a].voxelIndex[b+2].y;
	  branchList[a].voxelGradient[b].z = branchList[a].voxelGradient[b].z+BSplineCoe[0]*branchList[a].voxelIndex[b-1].z;
	  branchList[a].voxelGradient[b].z = branchList[a].voxelGradient[b].z+BSplineCoe[1]*branchList[a].voxelIndex[b].z;
	  branchList[a].voxelGradient[b].z = branchList[a].voxelGradient[b].z+BSplineCoe[2]*branchList[a].voxelIndex[b+1].z;
	  branchList[a].voxelGradient[b].z = branchList[a].voxelGradient[b].z+BSplineCoe[3]*branchList[a].voxelIndex[b+2].z;
	}
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
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
	  outImage->SetPixel(pixelIndex, reader->GetOutput()->GetPixel(pixelIndex));
	}

  //Write tree and corresponding gradient vector
  int c,d;
  for(a=0;a<branchCt;a++)
    {
      for(b=0;b<branchList[a].voxelCt;b++)
	{
	  pixelIndex[0]=int(branchList[a].voxelIndex[b].x);
	  pixelIndex[1]=int(branchList[a].voxelIndex[b].y);
	  pixelIndex[2]=int(branchList[a].voxelIndex[b].z);    
	  outImage->SetPixel(pixelIndex,2);
	}
    }

  a=4;
  for(c=1;c<branchList[a].voxelCt-2;c=c+10)
    {
      for(d=1;d<=10;d++)
	{
	  pixelIndex[0]=int(branchList[a].voxelIndex[c].x + d*branchList[a].voxelGradient[c].x);
	  pixelIndex[1]=int(branchList[a].voxelIndex[c].y + d*branchList[a].voxelGradient[c].y);
	  pixelIndex[2]=int(branchList[a].voxelIndex[c].z + d*branchList[a].voxelGradient[c].z);  
	  if((pixelIndex[0]>=0)&&(pixelIndex[0]<pRow)&&(pixelIndex[1]>=0)&&(pixelIndex[1]<pCol)&&(pixelIndex[2]>=0)&&(pixelIndex[2]<pSli))  
	    outImage->SetPixel(pixelIndex, 3);
	}

      //Show plane
      for(int x=branchList[a].voxelIndex[c].x-5; x<=branchList[a].voxelIndex[c].x+5; x++)
	for(int y=branchList[a].voxelIndex[c].y-5; y<=branchList[a].voxelIndex[c].y+5; y++)
	  {
	    int z = branchList[a].voxelIndex[c].z - (branchList[a].voxelGradient[c].x*(float(x)-branchList[a].voxelIndex[c].x)+branchList[a].voxelGradient[c].y*(float(y)-branchList[a].voxelIndex[c].y))/branchList[a].voxelGradient[c].z;
	    pixelIndex[0]=x;
	    pixelIndex[1]=y;
	    pixelIndex[2]=z;
	    if((pixelIndex[0]>=0)&&(pixelIndex[0]<pRow)&&(pixelIndex[1]>=0)&&(pixelIndex[1]<pCol)&&(pixelIndex[2]>=0)&&(pixelIndex[2]<pSli))  
	      outImage->SetPixel(pixelIndex, 4);
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

  return EXIT_SUCCESS;
}
