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
  int label;
  int generation; 
};

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputMaskFile outputLabelFile outputGenerationFile" << std::endl;
      return EXIT_FAILURE;
    }

  //Read tree
  VoxelType voxel;
  int centerlineCt = 0;
  FILE *treeFile = fopen(argv[1], "r");
  while(fscanf(treeFile, "%f %f %f %d %d", &voxel.x, &voxel.y, &voxel.z, &voxel.label, &voxel.generation) == 5)
    {
      centerlineCt++;
    }
  fclose(treeFile);

  VoxelType *centerLine = new VoxelType[centerlineCt];
  int a=0;
  treeFile = fopen(argv[1], "r");
  while(fscanf(treeFile, "%f %f %f %d %d", &voxel.x, &voxel.y, &voxel.z, &voxel.label, &voxel.generation) == 5)
    {
      centerLine[a].x = voxel.x;
      centerLine[a].y = voxel.y;
      centerLine[a].z = voxel.z;
      centerLine[a].label = voxel.label;
      centerLine[a].generation = voxel.generation;
      a++;
    }
  fclose(treeFile);
  
  //Read: 3D airway image
  //Write: 3D label image, 3D generation image
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writerLabel = WriterType::New();
  WriterType::Pointer writerGeneration = WriterType::New();
  //Parameters
  reader->SetFileName( argv[2] );
  writerLabel->SetFileName( argv[3] );
  writerGeneration->SetFileName( argv[4] );

  //Get image specs
  reader->Update();
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

  //Allocate output image
  ImageType::Pointer imageLabel = ImageType::New();
  imageLabel->SetRegions( region );
  imageLabel->SetSpacing( spacing );
  imageLabel->SetOrigin( origin );
  imageLabel->SetDirection( direction );
  imageLabel->Allocate();  
  ImageType::Pointer imageGeneration = ImageType::New();
  imageGeneration->SetRegions( region );
  imageGeneration->SetSpacing( spacing );
  imageGeneration->SetOrigin( origin );
  imageGeneration->SetDirection( direction );
  imageGeneration->Allocate();  

  ImageType::IndexType pixelIndex;
  int i, j, k, value;
  float iC, jC, kC;
  float generation, label, distance;
  float minDistance;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = reader->GetOutput()->GetPixel(pixelIndex);

	  if(value>0)
	    {
	      minDistance = 100000000;
	      for(a=0;a<centerlineCt;a++)
		{
		  iC = centerLine[a].x; 
		  jC = centerLine[a].y; 
		  kC = centerLine[a].z; 
		  distance = (i-iC)*(i-iC)+(j-jC)*(j-jC)+(k-kC)*(k-kC);
		  if(minDistance > distance)
		    {
		      minDistance = distance;
		      label = centerLine[a].label;
		      generation = centerLine[a].generation;
		    }
		}
	      imageLabel->SetPixel(pixelIndex, label);
	      imageGeneration->SetPixel(pixelIndex, generation);	      
	    }
	  else
	    {
	      imageLabel->SetPixel(pixelIndex, 0);
	      imageGeneration->SetPixel(pixelIndex, 0);
	    }
	}	  

  writerLabel->SetInput( imageLabel );
  writerGeneration->SetInput( imageGeneration );

  writerLabel->Update();
  writerGeneration->Update(); 

  return EXIT_SUCCESS;
}
