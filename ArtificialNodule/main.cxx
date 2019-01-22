#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <cstdlib>
#include <ctime>
int main( int argc, char * argv[] )
{
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile maskDTImageFile outputImageFile outputMaskFile depth radius [randSeed]" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New();  
  WriterType::Pointer writer = WriterType::New();
  WriterType::Pointer writerM = WriterType::New(); 
  //Parameters
  reader->SetFileName( argv[1] );
  reader->Update();
  readerM->SetFileName( argv[2] );
  readerM->Update();  
  writer->SetFileName( argv[3] );
  writerM->SetFileName( argv[4] );  
  float depth = atof(argv[5]);
  float radius = atof(argv[6]);
  if(argc == 8)
    srand(atoi(argv[7]));
  else
    srand((unsigned)time(0));
 
  //Set image specs
  ImageType::SizeType size = reader->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
  ImageType::RegionType region;
  region.SetSize( size ); 
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( reader->GetOutput()->GetSpacing() );
  image->SetOrigin( reader->GetOutput()->GetOrigin() );
  image->SetDirection( reader->GetOutput()->GetDirection() );
  image->Allocate();
  image->FillBuffer( 0.0 );
  ImageType::Pointer imageM = ImageType::New();
  imageM->SetRegions( region );
  imageM->SetSpacing( reader->GetOutput()->GetSpacing() );
  imageM->SetOrigin( reader->GetOutput()->GetOrigin() );
  imageM->SetDirection( reader->GetOutput()->GetDirection() );
  imageM->Allocate();
  imageM->FillBuffer( 0.0 );
  //Set values
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::IndexType pixelIndex;
  int i, j, k, value, valueM;

  //Count total "shell" points
  int totalShell = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0] = i;
	  pixelIndex[1] = j;
	  pixelIndex[2] = k;	  
	  valueM = readerM->GetOutput()->GetPixel(pixelIndex);
	  if((-valueM>depth-1)&&(-valueM<depth+0.5))
	    totalShell++;
	}
  std::cout<<"Total shell points count:"<<totalShell<<std::endl;

  // Randomly pick the Nth mask point as center
  // N is a random number between 1 and total "shell" point count 
  int centerCtN = 1 + rand() % static_cast<int>(totalShell);
  std::cout<<"Select the "<<centerCtN<<"-th point"<<std::endl;

  // Find center point
  int centerCt = 0;
  int centerX, centerY, centerZ;
  bool stopFlag = false;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0] = i;
	  pixelIndex[1] = j;
	  pixelIndex[2] = k;
	  //Set initial output as input
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  image->SetPixel(pixelIndex, value);
	  //Find the N-th point
	  valueM = readerM->GetOutput()->GetPixel(pixelIndex);
	  if((-valueM>depth-1)&&(-valueM<depth+0.5)&&!stopFlag)
	    {
	      if(centerCt<centerCtN)
		centerCt++;
	      else
		{
		  centerX = i;
		  centerY = j;
		  centerZ = k;
		  stopFlag = true;
		}  
	    }
	}
  std::cout<<centerX<<" "<<centerY<<" "<<centerZ<<": "<<depth<<" "<<radius<<std::endl;
  
  //Assign value around selected center
  int minX, maxX, minY, maxY, minZ, maxZ;
  minX = centerX - radius/spacing[0];
  minY = centerY - radius/spacing[1];
  minZ = centerZ - radius/spacing[2];
  if(minX < 0)
    minX = 0;
  if(minY < 0)
    minY = 0;
  if(minZ < 0)
    minZ = 0;
  maxX = centerX + radius/spacing[0];
  maxY = centerY + radius/spacing[1];
  maxZ = centerZ + radius/spacing[2];
  if(maxX > pRow)
    maxX = pRow;
  if(maxY > pCol)
    maxY = pCol;
  if(maxZ > pSli)
    maxZ = pSli;
  
  double dist;
  for(i=minX;i<maxX;i++)
    for(j=minY;j<maxY;j++)
      for(k=minZ;k<maxZ;k++)
	{
	  dist = (double(i)-centerX)*(double(i)-centerX)*spacing[0]*spacing[0];
	  dist = dist+(double(j)-centerY)*(double(j)-centerY)*spacing[1]*spacing[1];
	  dist = dist+(double(k)-centerZ)*(double(k)-centerZ)*spacing[2]*spacing[2];
	  //std::cout<<dist<<std::endl;
	  if(dist<=radius*radius)
	    {
	      pixelIndex[0] = i;
	      pixelIndex[1] = j;
	      pixelIndex[2] = k;
	      valueM = readerM->GetOutput()->GetPixel(pixelIndex);
	      if(valueM<2)
		{
		  value = 30 + rand() % static_cast<int>(70);
		  image->SetPixel(pixelIndex, value);
		  imageM->SetPixel(pixelIndex, 1);
		}
	    }
	}
  
  writer->SetInput( image );
  writer->Update();
  writerM->SetInput( imageM );
  writerM->Update(); 

  return EXIT_SUCCESS;
}
