#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <iostream>

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile variance upperThreshold lowerThreshod" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef double InputPixelType;
  typedef unsigned char OutputPixelType;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::CannyEdgeDetectionImageFilter< InputImageType, InputImageType > EdgeFilterType;
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;


  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  EdgeFilterType::Pointer edger = EdgeFilterType::New();
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters and Pipeline
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  float variance = atof( argv[3] );
  float upperThreshold = atof( argv[4] );
  float lowerThreshold = atof( argv[5] );
  edger->SetInput( reader->GetOutput() );
  edger->SetVariance( variance );
  edger->SetUpperThreshold( upperThreshold );
  edger->SetLowerThreshold( lowerThreshold );
  rescaler->SetInput( edger->GetOutput() );
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 255 );
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
