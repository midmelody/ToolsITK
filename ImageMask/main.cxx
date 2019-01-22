#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile maskImageFile outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef double PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > BinaryThresholdFilterType; 
  typedef itk::MultiplyImageFilter <ImageType, ImageType > MultiplyImageFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer readerFI = ReaderType::New();
  ReaderType::Pointer readerMI = ReaderType::New();
  BinaryThresholdFilterType::Pointer maskThresholdFilter = BinaryThresholdFilterType::New();
  MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  readerFI->SetFileName( argv[1] );
  readerMI->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
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

  //Write output  
  writer->SetInput( multiplyFilter->GetOutput() );
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
