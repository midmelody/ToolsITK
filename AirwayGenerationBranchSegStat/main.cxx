#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkPoint.h"
#include <iostream>

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
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputMaskImageFile " << std::endl;
      return EXIT_FAILURE;
    }
  
  //Three vectors record: generation, volume, and length for each airway branch segment (label)
  int genRec[1000];
  int volRec[1000];
  float lenRec[1000];
  int i, j, k;
  for(i=0; i<1000; i++)
    {
      genRec[i] = 0;
      volRec[i] = 0;
      lenRec[i] = 0.0;
    }
  
  //Read tree and assign generation record
  VoxelType voxel;
  int generationCt = 0;
  FILE *treeFile = fopen(argv[1], "r");
  while(fscanf(treeFile, "%f %f %f %d %d", &voxel.x, &voxel.y, &voxel.z, &voxel.label, &voxel.generation) == 5)
    {
      //Record total generation
      if(generationCt<voxel.generation)
	generationCt = voxel.generation;
      //Record corresponding generation for each label
      genRec[voxel.label] = voxel.generation;
    }
  fclose(treeFile);

  //Read image and assign volume record
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  //Filters
  ImageType::Pointer image = ImageType::New();
  ReaderType::Pointer readerLabel = ReaderType::New();
  //Parameters
  readerLabel->SetFileName( argv[2] );
  //Get image size
  readerLabel->Update();
  image = readerLabel->GetOutput();
  ImageType::SizeType size = image->GetRequestedRegion().GetSize();
  ImageType::SpacingType spacing = image->GetSpacing(); 
  float volUnit = spacing[0]*spacing[1]*spacing[2];
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::IndexType pixelIndex;
  int label;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  label = readerLabel->GetOutput()->GetPixel(pixelIndex);
	  if(label>0)
	    volRec[label] = volRec[label]+1;
	}	  

  //Read tree and image and assign length record
  ImageType::IndexType preIdx, curIdx, firstIdx, lastIdx;
  ImageType::PointType prePhy, curPhy, firstPhy, lastPhy;
  preIdx.Fill(0);
  curIdx.Fill(0);
  firstIdx.Fill(0);
  lastIdx.Fill(0);
  image->TransformIndexToPhysicalPoint(preIdx,prePhy);
  image->TransformIndexToPhysicalPoint(curIdx,curPhy);
  image->TransformIndexToPhysicalPoint(firstIdx,firstPhy);
  image->TransformIndexToPhysicalPoint(lastIdx,lastPhy);
  int preVoxelLab = 0;
  int preVoxelGen = 0;
  float tempLength;
  treeFile = fopen(argv[1], "r");
  while(fscanf(treeFile, "%f %f %f %d %d", &voxel.x, &voxel.y, &voxel.z, &voxel.label, &voxel.generation) == 5)
    {     
      //Set current point
      curIdx[0] = voxel.x;
      curIdx[1] = voxel.y;
      curIdx[2] = voxel.z;
      //Initialize length stat at beginning of each different branch segment (label)
      if(voxel.label!=preVoxelLab)
	{
	  //Record preIdx as lastIdx for previous branch if there is a "previous" branch, i.e. preVoxelGen>0
	  if(preVoxelGen>0)
	    {
	      lastIdx[0] = preIdx[0];  	  
	      lastIdx[1] = preIdx[1]; 
	      lastIdx[2] = preIdx[2]; 
	    }
	  //Set new preIdx to current index
	  preIdx[0] = voxel.x;
	  preIdx[1] = voxel.y;
	  preIdx[2] = voxel.z;
	  //Set new firstIdx to current index
	  firstIdx[0] = voxel.x;
	  firstIdx[1] = voxel.y;
	  firstIdx[2] = voxel.z;
	  //Set new preVoxelGen and preVoxelLab
	  preVoxelGen = voxel.generation;
	  preVoxelLab = voxel.label;
	}
      //Compute current distance in physical space
      image->TransformIndexToPhysicalPoint(preIdx,prePhy);
      image->TransformIndexToPhysicalPoint(curIdx,curPhy);
      tempLength = curPhy.EuclideanDistanceTo(prePhy);
      lenRec[voxel.label] = lenRec[voxel.label] + tempLength;
      //Set previous point
      preIdx[0] = curIdx[0];
      preIdx[1] = curIdx[1];
      preIdx[2] = curIdx[2];      
    }
  fclose(treeFile);

  //Output Stats
  for(i=0;i<1000;i++)
    {
      if( (genRec[i]>0)&&(genRec[i]<10) )
	{
	  std::cout<<"Generation "<<genRec[i];
	  std::cout<<" Volume "<<volRec[i]*volUnit; 
	  std::cout<<" Length "<<lenRec[i]<<std::endl;
	}
    }

  return EXIT_SUCCESS;
}
