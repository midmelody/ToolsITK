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
  float meanInt;
  float varInt;
  float maxInt;
  float minInt;
  float volume;
  float lengthCount;
  float length;
  float lengthPhy;
  float lengthDirect;
  float lengthDirectPhy;
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
      statGeneration[g].lengthCount = 0;
      statGeneration[g].length = 0;
      statGeneration[g].lengthPhy = 0;
      statGeneration[g].lengthDirect = 0;
      statGeneration[g].lengthDirectPhy = 0;
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
  ImageType::Pointer image = ImageType::New();
  ReaderType::Pointer readerLabel = ReaderType::New();
  
  //Parameters
  readerImage->SetFileName( argv[3] );
  readerLabel->SetFileName( argv[2] );
  int maxValue = atoi(argv[4]);
  //Get image size
  readerImage->Update();
  readerLabel->Update();
  image = readerImage->GetOutput();
  ImageType::SizeType size = image->GetRequestedRegion().GetSize();
  ImageType::SpacingType spacing = image->GetSpacing(); 
  float volUnit = spacing[0]*spacing[1]*spacing[2];
  float lengthUnit = spacing[0]*spacing[0]+spacing[1]*spacing[1]+spacing[2]*spacing[2];
  lengthUnit = sqrt(lengthUnit);
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
  while(fscanf(treeFile, "%f %f %f %d %d", &voxel.x, &voxel.y, &voxel.z, &voxel.label, &voxel.generation) == 5)
    {     
      //length by counting with max spacing setting
      statGeneration[voxel.generation-1].lengthCount++;

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
	      std::cout<<"First idx for branch "<<preVoxelLab<<": "<<firstIdx<<std::endl;	 
	      lastIdx[0] = preIdx[0];  	  
	      lastIdx[1] = preIdx[1]; 
	      lastIdx[2] = preIdx[2]; 
	      std::cout<<"Last idx for branch "<<preVoxelLab<<": "<<lastIdx<<std::endl;	      
	      //Compute direct length for last segment
	      tempLength = (firstIdx[0]-lastIdx[0])*spacing[0]*(firstIdx[0]-lastIdx[0])*spacing[0];
	      tempLength = tempLength + (firstIdx[1]-lastIdx[1])*spacing[1]*(firstIdx[1]-lastIdx[1])*spacing[1];
	      tempLength = tempLength + (firstIdx[2]-lastIdx[2])*spacing[2]*(firstIdx[2]-lastIdx[2])*spacing[2];
	      tempLength = sqrt(tempLength);
	      std::cout<<"Straight length: "<<tempLength<<std::endl;
	      statGeneration[preVoxelGen-1].lengthDirect = statGeneration[preVoxelGen-1].lengthDirect+tempLength;
	      //Compute in physical space
	      image->TransformIndexToPhysicalPoint(firstIdx,firstPhy);
	      image->TransformIndexToPhysicalPoint(lastIdx,lastPhy);
	      tempLength = firstPhy.EuclideanDistanceTo(lastPhy);
	      std::cout<<"Straight length in physical space: "<<tempLength<<std::endl;
	      statGeneration[preVoxelGen-1].lengthDirectPhy = statGeneration[preVoxelGen-1].lengthDirectPhy+tempLength;
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
      //Compute current distance
      tempLength = (curIdx[0]-preIdx[0])*spacing[0]*(curIdx[0]-preIdx[0])*spacing[0];
      tempLength = tempLength + (curIdx[1]-preIdx[1])*spacing[1]*(curIdx[1]-preIdx[1])*spacing[1];
      tempLength = tempLength + (curIdx[2]-preIdx[2])*spacing[2]*(curIdx[2]-preIdx[2])*spacing[2];
      tempLength = sqrt(tempLength);
      statGeneration[voxel.generation-1].length = statGeneration[voxel.generation-1].length+tempLength;
      //Compute in physical space
      image->TransformIndexToPhysicalPoint(preIdx,prePhy);
      image->TransformIndexToPhysicalPoint(curIdx,curPhy);
      tempLength = curPhy.EuclideanDistanceTo(prePhy);
      statGeneration[voxel.generation-1].lengthPhy = statGeneration[voxel.generation-1].lengthPhy+tempLength;
      
      //Set previous point
      preIdx[0] = curIdx[0];
      preIdx[1] = curIdx[1];
      preIdx[2] = curIdx[2];      

      std::cout<<preIdx<<curIdx<<firstIdx<<lastIdx<<std::endl;
    }
  fclose(treeFile);

  //std::cout<<"^^^^^^^^^^^^"<<std::endl;
  //for(g=0;g<generationCt;g++)
  //  std::cout<<statGeneration[g].lengthCount<<" "<<statGeneration[g].length<<std::cout<<std::endl;


  /*
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
	      
	      //statGeneration[label-1].meanInt = statGeneration[label-1].meanInt + value;
	      //if(statGeneration[label-1].maxInt < value) 
	      //  statGeneration[label-1].maxInt = value;
	      //if(statGeneration[label-1].minInt > value)
		//statGeneration[label-1].minInt = value;
	      
	      statGeneration[label-1].volume = statGeneration[label-1].volume + 1;
	      
	      //for(g=0;g<generationCt;g++)
	      //std::cout<<statGeneration[g].volume<<" ";
	      //std::cout<<std::endl;
	    }
	}	  

*/

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
      //std::cout<<"Volume "<< statGeneration[g].volume*volUnit <<std::endl;
      std::cout<<"LengthCount "<< statGeneration[g].lengthCount*lengthUnit <<std::endl;
      std::cout<<"Length "<< statGeneration[g].length <<std::endl; 
      std::cout<<"LengthPhy "<<      statGeneration[g].lengthPhy <<std::endl; 
      std::cout<<"LengthDirect "<<    statGeneration[g].lengthDirect <<std::endl; 
      std::cout<<"LengthDirectPhy "<<    statGeneration[g].lengthDirectPhy <<std::endl; 
   }

  return EXIT_SUCCESS;
}
