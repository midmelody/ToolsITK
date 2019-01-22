#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile startX startY endX endY" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 2;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  
  //Parameters
  reader->SetFileName( argv[1] );
  reader->Update();
  writer->SetFileName( argv[2] );
  int startX = atoi(argv[3]);
  int startY = atoi(argv[4]);
  int endX = atoi(argv[5]);
  int endY = atoi(argv[6]);
  int sizeX = endX-startX+1;
  int sizeY = endY-startY+1;

  //Set image specs
  ImageType::SizeType  size;
  size[0] = sizeX;
  size[1] = sizeY;
  ImageType::RegionType region;
  region.SetSize( size ); 
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->Allocate(); 
  int pRow, pCol;
  pRow = size[0];
  pCol = size[1];
  
  //Read image 1 by 1
  ImageType::IndexType pixelInIndex, pixelOutIndex;
  int i, j, k, value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelInIndex[0] = i+startX;
	pixelInIndex[1] = j+startY;	    
	value = reader->GetOutput()->GetPixel(pixelInIndex);
	pixelOutIndex[0] = i;
	pixelOutIndex[1] = j;
	image->SetPixel(pixelOutIndex, value);
      }

  writer->SetInput( image );
  writer->Update();

  return EXIT_SUCCESS;
}
