#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " featureImageFile maskImageFile seedImageFile ratio" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > BinaryThresholdFilterType; 
  typedef itk::MultiplyImageFilter< ImageType, ImageType, ImageType > MultiplyFilterType; 
  typedef itk::StatisticsImageFilter< ImageType > StatFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer readerFI = ReaderType::New();
  ReaderType::Pointer readerMI = ReaderType::New();
  BinaryThresholdFilterType::Pointer maskThresholdFilter = BinaryThresholdFilterType::New();
  MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
  StatFilterType::Pointer statFilter = StatFilterType::New();
  BinaryThresholdFilterType::Pointer seedThresholdFilter = BinaryThresholdFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  readerFI->SetFileName( argv[1] );
  readerMI->SetFileName( argv[2] );
  writer->SetFileName( argv[3] ); 
  float ratio = atof( argv[4] );

  //Pipeline
  try
    {
      readerFI->Update();
      readerMI->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  maskThresholdFilter->SetInput( readerMI->GetOutput() );
  maskThresholdFilter->SetOutsideValue( 0 );
  maskThresholdFilter->SetInsideValue( 1 );
  maskThresholdFilter->SetLowerThreshold( 1 ) ;
  maskThresholdFilter->SetUpperThreshold( 10 );

  multiplyFilter->SetInput1( maskThresholdFilter->GetOutput() );
  multiplyFilter->SetInput2( readerFI->GetOutput() ); 
  
  statFilter->SetInput( multiplyFilter->GetOutput() ); 
  statFilter->Update();
 
  float thresholdL, thresholdU;
  thresholdU = statFilter->GetMaximum();
  std::cout<<"Max "<<thresholdU<<std::endl;
  thresholdL = thresholdU * ratio;
  seedThresholdFilter->SetInput( multiplyFilter->GetOutput() );
  seedThresholdFilter->SetOutsideValue( 0 );
  seedThresholdFilter->SetInsideValue( 1 );
  seedThresholdFilter->SetLowerThreshold( thresholdL );
  seedThresholdFilter->SetUpperThreshold( thresholdU );

  //Write output
  writer->SetInput( seedThresholdFilter->GetOutput() );
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

 
