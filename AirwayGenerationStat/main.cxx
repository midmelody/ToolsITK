#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
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

struct StatType
{
  float volume;
  float length;
};

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputMaskImageFile [outputGenImageFile]" << std::endl;
      return EXIT_FAILURE;
    }

  //Read tree
  VoxelType voxel;
  int generationCt = 0;
  FILE *treeFile = fopen(argv[1], "r");
  int genRec[1000]; 
  for(int i=0;i<1000;i++)
    genRec[i] = 0;
  while(fscanf(treeFile, "%f %f %f %d %d", &voxel.x, &voxel.y, &voxel.z, &voxel.label, &voxel.generation) == 5)
    {
      //Record total generation
      if(generationCt<voxel.generation)
	generationCt = voxel.generation;
      //Record corresponding generation for each label
      genRec[voxel.label] = voxel.generation;
    }
  fclose(treeFile);

  StatType *statGeneration = new StatType[generationCt];
  //Initialize
  int g;
  for(g=0;g<generationCt;g++)
    {
      statGeneration[g].volume = 0;
      statGeneration[g].length = 0;
    }

  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ImageType::Pointer image = ImageType::New();
  ImageType::Pointer imageOut = ImageType::New();
  ReaderType::Pointer readerLabel = ReaderType::New();
  WriterType::Pointer writerGen = WriterType::New();

  //Parameters
  readerLabel->SetFileName( argv[2] );
  if(argc==4)
    writerGen->SetFileName( argv[3] );
  
  //Get image size
  readerLabel->Update();
  image = readerLabel->GetOutput();
  ImageType::SizeType size = image->GetRequestedRegion().GetSize();
  ImageType::SpacingType spacing = image->GetSpacing(); 
  ImageType::PointType origin = image->GetOrigin(); 
  ImageType::DirectionType direction = image->GetDirection();
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image if needed
  if(argc==4)
    {
      imageOut->SetRegions( region );
      imageOut->SetSpacing( spacing );
      imageOut->SetOrigin( origin );
      imageOut->SetDirection( direction );
      imageOut->Allocate();  
      imageOut->FillBuffer(0);
    }

  float volUnit = spacing[0]*spacing[1]*spacing[2];
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 

  //Run stats
  //Length count and length
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
      statGeneration[voxel.generation-1].length = statGeneration[voxel.generation-1].length+tempLength;
      
      //Set previous point
      preIdx[0] = curIdx[0];
      preIdx[1] = curIdx[1];
      preIdx[2] = curIdx[2];      
    }
  fclose(treeFile);

  ImageType::IndexType pixelIndex;
  int i, j, k, value, label;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  label = readerLabel->GetOutput()->GetPixel(pixelIndex);
	  if(label>0)
	    {
	      int curGen = genRec[label]-1;
	      if(curGen>=0)
		statGeneration[genRec[label]-1].volume = statGeneration[genRec[label]-1].volume + 1;
	      if(argc==4)
		imageOut->SetPixel(pixelIndex,genRec[label]);
	    }
	}	  

  //Output Stats
  for(g=0;g<generationCt;g++)
    {
      std::cout<<"----------"<<std::endl;
      std::cout<<"Generation "<<g+1<<std::endl;
      std::cout<<"Volume "<< statGeneration[g].volume*volUnit <<std::endl;
      std::cout<<"Length "<< statGeneration[g].length <<std::endl; 
    }

  //Write output
  if(argc==4)
    {
      writerGen->SetInput(imageOut);
      writerGen->Update();
    }
  return EXIT_SUCCESS;
}
