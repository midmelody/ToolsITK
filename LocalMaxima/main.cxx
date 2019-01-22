#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkValuedRegionalMaximaImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>

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
  typedef double InputPixelType;
  typedef unsigned char OutputPixelType;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::ValuedRegionalMaximaImageFilter < InputImageType, InputImageType > ValuedRegionalMaximaImageFilter;
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ValuedRegionalMaximaImageFilter::Pointer filter = ValuedRegionalMaximaImageFilter::New ();
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters and Pipeline
  reader->SetFileName( argv[1] );
  filter->SetInput( reader->GetOutput() );
  rescaler->SetInput( filter->GetOutput() );
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 255 );
  writer->SetFileName( argv[2] );
  writer->SetInput( rescaler->GetOutput() );

  //Write output  
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
