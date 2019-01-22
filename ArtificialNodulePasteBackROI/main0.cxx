#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile inputROIImageFile outputImageFile refROIName" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer readerIMG = ReaderType::New();
  ReaderType::Pointer readerROI = ReaderType::New();  
  WriterType::Pointer writer = WriterType::New(); 

  //Parameters
  readerIMG->SetFileName( argv[1] );
  readerIMG->Update();
  readerROI->SetFileName( argv[2] );
  readerROI->Update();
  writer->SetFileName( argv[3] );
  
  //Parse for center location
  char *filenameRef;
  int centerX, centerY, centerZ;
  int caseID, nodID, rad;
  char *lastdot;
  filenameRef = basename(argv[4]);
  lastdot = strrchr(filenameRef, '.');
  if (lastdot != NULL)
    *lastdot = '\0'; 
  sscanf(filenameRef, "%d_%d_%d_%d_%d_%d", &caseID, &nodID, &rad, &centerX, &centerY, &centerZ);
  std::cout<<"Center: "<<centerX<<" "<<centerY<<" "<<centerZ<<std::endl;

  //Get image specs
  ImageType::SizeType sizeIMG = readerIMG->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = sizeIMG[0];
  int pCol = sizeIMG[1];
  int pSli = sizeIMG[2];
  ImageType::SizeType sizeROI = readerROI->GetOutput()->GetRequestedRegion().GetSize();
  //Set output 
  ImageType::Pointer image = ImageType::New();
  image = readerIMG->GetOutput();
  
  //Replace around the center point
  int minX, maxX, minY, maxY, minZ, maxZ;
  minX = centerX - sizeROI[0]/2;
  minY = centerY - sizeROI[1]/2;
  minZ = centerZ - sizeROI[2]/2;
  if(minX < 0)
    minX = 0;
  if(minY < 0)
    minY = 0;
  if(minZ < 0)
    minZ = 0;
  maxX = minX + sizeROI[0];
  maxY = minY + sizeROI[1];
  maxZ = minZ + sizeROI[2];
  if(maxX > pRow)
    maxX = pRow;
  if(maxY > pCol)
    maxY = pCol;
  if(maxZ > pSli)
    maxZ = pSli;
 
  //Set output image value around selected center
  ImageType::IndexType pixelIndex, pixelIndexROI;
  int i, j, k, value;
  for(i=minX;i<maxX-1;i++)
    for(j=minY;j<maxY-1;j++)
      for(k=minZ;k<maxZ-1;k++)
	{
	  pixelIndexROI[0] = i-minX;
	  pixelIndexROI[1] = j-minY;	  
	  pixelIndexROI[2] = k-minZ;
	  value = readerROI->GetOutput()->GetPixel(pixelIndexROI);	  
	  pixelIndex[0] = i;
	  pixelIndex[1] = j;
	  pixelIndex[2] = k;
	  image->SetPixel(pixelIndex, value);
	}
 
  writer->SetInput( image );
  writer->Update();

  return EXIT_SUCCESS;
}
