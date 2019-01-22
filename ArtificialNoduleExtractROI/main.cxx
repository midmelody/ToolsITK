#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <cstdlib>
#include <ctime>
int main( int argc, char * argv[] )
{
  if( argc < 8 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile maskImageFile maskDTImageFile outputImageFile(no ext) outputMaskFile(no ext) depth size [randSeed]" << std::endl;
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
  ReaderType::Pointer readerMDT = ReaderType::New();
  
  WriterType::Pointer writer = WriterType::New(); 
  WriterType::Pointer writerM = WriterType::New();
  //Parameters
  reader->SetFileName( argv[1] );
  reader->Update();
  readerM->SetFileName( argv[2] );
  readerM->Update();  
  readerMDT->SetFileName( argv[3] );
  readerMDT->Update();  
  float depth = atof(argv[6]);
  float size = atof(argv[7])*4;
  if(argc == 9)
    srand(atoi(argv[8])*(unsigned)time(0));
  else
    srand((unsigned)time(0));
 
  //Set output image specs - ROI image and a mask image with same size as input
  ImageType::SizeType sizeIn = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = sizeIn[0];
  int pCol = sizeIn[1];
  int pSli = sizeIn[2];
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
  float spaX = spacing[0];
  float spaY = spacing[1]; 
  float spaZ = spacing[2];
  //std::cout<<"Spacing: "<<spaX<<" "<<spaY<<" "<<spaZ<<std::endl;
  ImageType::SizeType sizeOut;
  sizeOut[0] = size/spaX;
  sizeOut[1] = size/spaY;  
  sizeOut[2] = size/spaZ+1;
  std::cout<<"Size: "<<sizeOut<<std::endl; 
  ImageType::RegionType region;
  region.SetSize( sizeOut );
  
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( reader->GetOutput()->GetSpacing() );
  image->SetOrigin( reader->GetOutput()->GetOrigin() );
  image->SetDirection( reader->GetOutput()->GetDirection() );
  image->Allocate();
  image->FillBuffer( -1000.0 );
  ImageType::Pointer imageM = ImageType::New();
  imageM->SetRegions( region );
  imageM->SetSpacing( reader->GetOutput()->GetSpacing() );
  imageM->SetOrigin( reader->GetOutput()->GetOrigin() );
  imageM->SetDirection( reader->GetOutput()->GetDirection() );
  imageM->Allocate();
  imageM->FillBuffer( 0.0 );
  
  //Set values
  ImageType::IndexType pixelIndex, pixelIndexOut;
  int i, j, k, value, valueM, valueMDT;

  //Count total "shell" points
  int totalShell = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0] = i;
	  pixelIndex[1] = j;
	  pixelIndex[2] = k;	  
	  valueMDT = readerMDT->GetOutput()->GetPixel(pixelIndex);
	  if((-valueMDT>depth-1)&&(-valueMDT<depth+0.5))
	    totalShell++;
	}
  //std::cout<<"Total shell points count:"<<totalShell<<std::endl;

  // Randomly pick the Nth mask point as center
  // N is a random number between 1 and total "shell" point count 
  int centerCtN = 1 + rand() % static_cast<int>(totalShell);
  //std::cout<<"Select the "<<centerCtN<<"-th point"<<std::endl;

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
	  //Find the N-th point
	  valueMDT = readerMDT->GetOutput()->GetPixel(pixelIndex);
	  if((-valueMDT>depth-1)&&(-valueMDT<depth+0.5)&&!stopFlag)
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
  std::cout<<"Center: "<<centerX<<" "<<centerY<<" "<<centerZ<<std::endl;
  std::cout<<"Depth: "<<depth<<std::endl;
  
  //Assign value around selected center
  int minX, maxX, minY, maxY, minZ, maxZ;
  minX = centerX - sizeOut[0]/2;
  minY = centerY - sizeOut[1]/2;
  minZ = centerZ - sizeOut[2]/2;
  if(minX < 0)
    minX = 0;
  if(minY < 0)
    minY = 0;
  if(minZ < 0)
    minZ = 0;
  maxX = minX + sizeOut[0];
  maxY = minY + sizeOut[1];
  maxZ = minZ + sizeOut[2];
  if(maxX > pRow)
    maxX = pRow;
  if(maxY > pCol)
    maxY = pCol;
  if(maxZ > pSli)
    maxZ = pSli;
  //std::cout<<"Min: "<<minX<<" "<<minY<<" "<<minZ<<std::endl;
  //std::cout<<"Max: "<<maxX<<" "<<maxY<<" "<<maxZ<<std::endl;
  for(i=minX;i<maxX;i++)
    for(j=minY;j<maxY;j++)
      for(k=minZ;k<maxZ;k++)
	{
	  pixelIndex[0] = i;
	  pixelIndex[1] = j;
	  pixelIndex[2] = k;
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  valueM = readerM->GetOutput()->GetPixel(pixelIndex);	  
	  pixelIndexOut[0] = i-minX;
	  pixelIndexOut[1] = j-minY;	  
	  pixelIndexOut[2] = k-minZ;
	  image->SetPixel(pixelIndexOut, value);
	  imageM->SetPixel(pixelIndexOut, valueM);	  
	}

  char filename[100];
  sprintf(filename, "%s_%d_%d_%d.hdr", argv[4],centerX,centerY,centerZ);
  char filenameM[100];
  sprintf(filenameM, "%s_%d_%d_%d.hdr", argv[5],centerX,centerY,centerZ);
  
  writer->SetFileName( filename ); 
  writer->SetInput( image );
  writer->Update();


  writerM->SetFileName( filenameM ); 
  writerM->SetInput( imageM );
  writerM->Update();
  return EXIT_SUCCESS;
}
