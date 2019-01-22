#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkFlipImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 2;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::FlipImageFilter< ImageType > FlipImageFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  
  //Pipeline
  reader->Update();

  FlipImageFilterType::FlipAxesArrayType flipAxes;
  flipAxes[0] = true;
  flipAxes[1] = false;
  flipFilter->SetInput(reader->GetOutput());
  flipFilter->SetFlipAxes(flipAxes);

  writer->SetInput( flipFilter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
