#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

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
  float meanInt;
  float varInt;
  float maxInt;
  float minInt;
  float volume;
  float lengthCT;
  float length;
};

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputMaskFile inputImageFile maxValue" << std::endl;
      return EXIT_FAILURE;
    }

  //Read tree
  VoxelType voxel;
  int generationCt = 0;
  FILE *treeFile = fopen(argv[1], "r");
  while(fscanf(treeFile, "%f %f %f %d %d", &voxel.x, &voxel.y, &voxel.z, &voxel.label, &voxel.generation) == 5)
    {
      if(generationCt<voxel.generation)
	generationCt = voxel.generation;
    }
  fclose(treeFile);
  
  std::cout<<"GenerationCt "<<generationCt<<std::endl;

  StatType *statGeneration = new StatType[generationCt];
  //Initialize
  int g;
  for(g=0;g<generationCt;g++)
    {
      statGeneration[g].meanInt = 0;
      statGeneration[g].varInt = 0;
      statGeneration[g].maxInt = -100000000;
      statGeneration[g].minInt = 100000000;
      statGeneration[g].volume = 0;
      statGeneration[g].lengthCT = 0;
      statGeneration[g].length = 0;
    }

  treeFile = fopen(argv[1], "r");
  
  //Read: 3D labeled airway/wall image; 3D original image
  //Output: statistics
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  //Filters
  ReaderType::Pointer readerImage = ReaderType::New();
  ReaderType::Pointer readerLabel = ReaderType::New();
  
  //Parameters
  readerImage->SetFileName( argv[3] );
  readerLabel->SetFileName( argv[2] );
  int maxValue = atoi(argv[4]);
  //Get image size
  readerImage->Update();
  readerLabel->Update();
  ImageType::SizeType size = readerImage->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::SpacingType spacing = readerImage->GetOutput()->GetSpacing(); 
  float volUnit = spacing[0]*spacing[1]*spacing[2];
  float lengthUnit = spacing[0]*spacing[0]+spacing[1]*spacing[1]+spacing[2]*spacing[2];
  lengthUnit = sqrt(lengthUnit);
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
 
  //Run stats
  //Length count and length
  int voxX = 0;
  int voxY = 0;
  int voxZ = 0;
  int voxelLabel = 1;
  float tempLength;
  while(fscanf(treeFile, "%f %f %f %d %d", &voxel.x, &voxel.y, &voxel.z, &voxel.label, &voxel.generation) == 5)
    {     
      statGeneration[voxel.generation-1].lengthCT++;
      //Initialize length stat at beginning of each different branch segment (label)
      if(voxel.label!=voxelLabel)
	{
	  voxX = voxel.x;
	  voxY = voxel.y;
	  voxZ = voxel.z;
	  voxelLabel = voxel.label;
	}

      tempLength = (voxel.x-voxX)*spacing[0]*(voxel.x-voxX)*spacing[0];
      tempLength = tempLength + (voxel.y-voxY)*spacing[1]*(voxel.y-voxY)*spacing[1];
      tempLength = tempLength + (voxel.z-voxZ)*spacing[2]*(voxel.z-voxZ)*spacing[2];
      tempLength = sqrt(tempLength);
      statGeneration[voxel.generation-1].length = statGeneration[voxel.generation-1].length+tempLength;
      
      voxX = voxel.x;
      voxY = voxel.y;
      voxZ = voxel.z;      
    }
  fclose(treeFile);

  //std::cout<<"^^^^^^^^^^^^"<<std::endl;
  //for(g=0;g<generationCt;g++)
  //  std::cout<<statGeneration[g].lengthCT<<" "<<statGeneration[g].length<<std::cout<<std::endl;

  ImageType::IndexType pixelIndex;
  int i, j, k, value, label;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = readerImage->GetOutput()->GetPixel(pixelIndex);
	  label = readerLabel->GetOutput()->GetPixel(pixelIndex);
	  if((label>0)&&(label<=generationCt)&&(value<maxValue))
	    {
	      /*
	      statGeneration[label-1].meanInt = statGeneration[label-1].meanInt + value;
	      if(statGeneration[label-1].maxInt < value) 
	        statGeneration[label-1].maxInt = value;
	      if(statGeneration[label-1].minInt > value)
		statGeneration[label-1].minInt = value;
	      */
	      statGeneration[label-1].volume = statGeneration[label-1].volume + 1;
	      
	      //for(g=0;g<generationCt;g++)
	      //std::cout<<statGeneration[g].volume<<" ";
	      //std::cout<<std::endl;
	    }
	}	  



  /*
  for(g=0;g<generationCt;g++)
    {
      statGeneration[g].meanInt = statGeneration[g].meanInt/statGeneration[g].volume;
    }

  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = readerImage->GetOutput()->GetPixel(pixelIndex);
	  label = readerLabel->GetOutput()->GetPixel(pixelIndex);
	  if(label>0)
	    {
	      statGeneration[label-1].varInt = statGeneration[label-1].varInt+(value-statGeneration[label-1].meanInt)*(value-statGeneration[label-1].meanInt);
	    }
	}	 
  */
  //Output Stats
  for(g=0;g<generationCt;g++)
    {
      std::cout<<"----------"<<std::endl;
      std::cout<<"Generation "<<g+1<<std::endl;
      //std::cout<<"MeanIntensity "<< statGeneration[g].meanInt<<std::endl;
      //std::cout<<"VarIntensity "<< statGeneration[g].varInt/statGeneration[g].volume<<std::endl;
      //std::cout<<"MaxIntensity "<< statGeneration[g].maxInt <<std::endl;
      //std::cout<<"MinIntensity "<< statGeneration[g].minInt <<std::endl;
      std::cout<<"Volume "<< statGeneration[g].volume*volUnit <<std::endl;
      std::cout<<"LengthCount "<< statGeneration[g].lengthCT*lengthUnit <<std::endl;
      std::cout<<"Length "<< statGeneration[g].length <<std::endl; 
   }

  return EXIT_SUCCESS;
}
