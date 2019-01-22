#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileWriter.h"
#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef unsigned char OutPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< OutPixelType, Dimension > OutImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::BinaryThresholdImageFilter< ImageType, OutImageType > BinaryThresholdFilterType; 
  typedef itk::ImageFileWriter< OutImageType > WriterType;
    
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
  WriterType::Pointer writer = WriterType::New();
   
  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

  reader->Update();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::IndexType pixelIndex;
  pixelIndex[0] = (size[0]-1)/2;
  pixelIndex[1] = (size[1]-1)/2;
  pixelIndex[2] = (size[2]-1)/2;

  float threshold = reader->GetOutput()->GetPixel(pixelIndex);

  thresholder->SetInput(reader->GetOutput());
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetLowerThreshold( threshold );
  thresholder->SetUpperThreshold( threshold );

  writer->SetInput( thresholder->GetOutput() );
  try
    {
      writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }



  return EXIT_SUCCESS;
}
