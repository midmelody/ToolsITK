#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

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
  VoxelType * voxelTangent;
};

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputImageFile outputFileFolder " << std::endl;
      return EXIT_FAILURE;
    }

  //Read tree
  VoxelType index;
  VoxelType tangent;
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
  while(fscanf(treeFile, "%f %f %f %d %d %f %f %f", &index.x, &index.y, &index.z, &label, &parLabel, &tangent.x, &tangent.y, &tangent.z) == 8)
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

  while(fscanf(treeFile, "%f %f %f %d %d %f %f %f", &index.x, &index.y, &index.z, &label, &parLabel, &tangent.x, &tangent.y, &tangent.z) == 8)
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
      branchList[a].voxelTangent = new VoxelType[branchList[a].voxelCt];
      curVoxel[a] = 0;
    }
  treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(treeFile, "%f %f %f %d %d %f %f %f", &index.x, &index.y, &index.z, &label, &parLabel, &tangent.x, &tangent.y, &tangent.z) == 8)
    {
      branchList[label-1].voxelIndex[curVoxel[label-1]] = index;
      branchList[label-1].voxelTangent[curVoxel[label-1]] = tangent;
      curVoxel[label-1]++;
    }
  fclose(treeFile);
  
  //Resample original image at a point along centerline on orthogonal plane to tangent vector
  //ITK settings
  //Read: 3D original image
  //Write: sample 2D resampled image
  const unsigned int Dimension3D = 3;
  typedef float PixelType3D;
  typedef itk::Image< PixelType3D, Dimension3D > ImageType3D;
  typedef itk::ImageFileReader< ImageType3D > ReaderType;
  const unsigned int Dimension2D = 2;
  typedef unsigned char PixelType2D;
  typedef itk::Image< PixelType2D, Dimension2D > ImageType2D;
  typedef itk::ImageFileWriter< ImageType2D > WriterType;
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  //Parameters
  reader->SetFileName( argv[2] );

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
  //Get in image specs
  ImageType3D::SizeType  size3D = reader->GetOutput()->GetRequestedRegion().GetSize();

  //Set out image specs
  ImageType2D::SizeType  size;
  size[0] = 51;
  size[1] = 51;
  int pRow, pCol;
  pRow = size[0];
  pCol = size[1];
  ImageType2D::RegionType region;
  region.SetSize( size );
  //Allocate image for output
  ImageType2D::Pointer outImage = ImageType2D::New();
  outImage->SetRegions( region );
  outImage->Allocate();  
  //Initialize
  float * temp = new float[size[0]*size[1]];
  ImageType2D::IndexType pixelIndex2D;
  int i, j;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex2D[0]=i;
	pixelIndex2D[1]=j;
	temp[j*pRow+i] = 0;
	outImage->SetPixel(pixelIndex2D, 0);
      }

  //Take a sample
  ImageType3D::IndexType pixelIndex3D;
  float xf, yf, zf;
  int b;
  for(a=0;a<branchCt;a++)
    //for(a=0;a<1;a++)
    {
      for(b=0;b<branchList[a].voxelCt;b++)
	//for(b=0;b<1;b++)
	{ 
	  char fileName[100];
	  sprintf(fileName,"%s/%d_%d.png",argv[3],a,b);
	  writer->SetFileName( fileName );

	  index = branchList[a].voxelIndex[b];
	  tangent = branchList[a].voxelTangent[b]; 
  
	  //Compute Rotation system
	  float sin1, sin2, cos1, cos2;
	  if((tangent.y!=0)||(tangent.z!=0))
	    {
	      sin1 = tangent.y/sqrt(tangent.y*tangent.y+tangent.z*tangent.z);
	      cos1 = tangent.z/sqrt(tangent.y*tangent.y+tangent.z*tangent.z);
	      sin2 = tangent.x/sqrt(tangent.x*tangent.x+tangent.y*tangent.y+tangent.z*tangent.z);
	      cos2 = sqrt(tangent.y*tangent.y+tangent.z*tangent.z)/sqrt(tangent.x*tangent.x+tangent.y*tangent.y+tangent.z*tangent.z);
	    }
	  else
	    {
	      sin1 = 1;
	      cos1 = 0;
	      sin2 = 1;
	      cos2 = 0;
	    }
	  
	  //Plane orthogonal to tangent and through index
	  for(int x=-(pRow-1)/2; x<=(pRow-1)/2; x++)
	    for(int y=-(pCol-1)/2; y<=(pCol-1)/2; y++)
	      {
		//Get coordinate after rotation
		float xG = cos2*float(x);
		float yG = -sin1*sin2*float(x)+cos1*float(y);
		float zG = -cos1*sin2*float(x)-sin1*float(y);
		xG = xG + index.x;
		yG = yG + index.y;
		zG = zG + index.z;
		if(xG<0) xG = 0;
		if(xG>size3D[0]-1) xG = size3D[0]-1; 
		if(yG<0) yG = 0;
		if(yG>size3D[1]-1) yG = size3D[1]-1; 
		if(zG<0) zG = 0;
		if(zG>size3D[2]-1) zG = size3D[2]-1; 

		pixelIndex3D[0] = int(xG);
		pixelIndex3D[1] = int(yG);
		pixelIndex3D[2] = int(zG);

		pixelIndex2D[0] = x + (size[0]-1)/2;
		pixelIndex2D[1] = y + (size[1]-1)/2;
		temp[pixelIndex2D[1]*pRow+pixelIndex2D[0]] = reader->GetOutput()->GetPixel(pixelIndex3D);
	      }	

	  float min=100000, max=0;
	  float valueF;
	  for(int t=0; t<pRow*pCol; t++)
	    {
	      valueF = temp[t];
	      if(min>valueF)
		min = valueF;
	      if(max<valueF)
		max = valueF;
	    }

	  unsigned char value;
	  for(int x=0; x<size[0]; x++)
	    for(int y=0; y<size[1]; y++)
	      {
		pixelIndex2D[0] = x;
		pixelIndex2D[1] = y;
		valueF = (temp[y*pRow+x]-min)/(max-min)*255;
		value = round(valueF);
		outImage->SetPixel(pixelIndex2D,value);
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

	  outImage->DisconnectPipeline();


    
	}
    }


  return EXIT_SUCCESS;
}
