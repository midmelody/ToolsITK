#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkN4BiasFieldCorrectionImageFilter.h"
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
  typedef unsigned int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::OtsuThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typedef itk::N4BiasFieldCorrectionImageFilter<ImageType, ImageType, ImageType> CorrecterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New();
  ThresholderType::Pointer thresholder = ThresholderType::New();
  CorrecterType::Pointer correcter = CorrecterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
 
  //Pipeline
  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }

  thresholder->SetInput( reader->GetOutput() );
  thresholder->SetNumberOfHistogramBins( 200 );
  thresholder->SetInsideValue( 0 );
  thresholder->SetOutsideValue( 1 );

  CorrecterType::VariableSizeArrayType  maximumNumberOfIterations( 4 );
  maximumNumberOfIterations.Fill( 10 );
  correcter->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  correcter->SetNumberOfFittingLevels( 4 );
  correcter->SetConvergenceThreshold( 0.000001 );
  correcter->SetInput( reader->GetOutput() );
  correcter->SetMaskImage( thresholder->GetOutput() );

  //Write output  
  writer->SetInput( correcter->GetOutput() );
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
