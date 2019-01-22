#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile thre" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 2;
  typedef unsigned char PixelType;
  typedef itk::RGBPixel< PixelType > RGBPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;  
  typedef itk::Image< RGBPixelType, Dimension > RGBImageType;
  typedef itk::ImageFileReader< RGBImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  
  //Parameters
  reader->SetFileName( argv[1] );
  reader->Update();
  writer->SetFileName( argv[2] );
  int thre = atoi(argv[3]);
  //Set image specs
  ImageType::SizeType size;
  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize( size ); 
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->Allocate();
  
  int pRow, pCol;
  pRow = size[0];
  pCol = size[1];
  
  //Read image and output green channel
  ImageType::IndexType pixelIndex;
  int i, j, k;
  RGBPixelType value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex[0] = i;
	pixelIndex[1] = j;	    
	value = reader->GetOutput()->GetPixel(pixelIndex);
	//if((i==245)&&(j==228))
	// std::cout<<value<<std::endl;
	if((value[1]>thre)&&(value[0]<10)&&(value[2]<10))
	  image->SetPixel(pixelIndex, 255);
	else
	  image->SetPixel(pixelIndex, 0);	  
      }

  writer->SetInput( image );
  writer->Update();

  return EXIT_SUCCESS;
}
