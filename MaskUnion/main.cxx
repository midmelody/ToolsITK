#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOrImageFilter.h"
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
      std::cerr << argv[0] << " inputImageFile1 inputImageFile2 outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::OrImageFilter< ImageType, ImageType, ImageType > OrFilterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > BinaryThresholdFilterType; 
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  OrFilterType::Pointer orFilter = OrFilterType::New();
  BinaryThresholdFilterType::Pointer maskThresholdFilter = BinaryThresholdFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
 
  //Pipeline
  try
    {
      reader1->Update();
      reader2->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  orFilter->SetInput1(reader1->GetOutput());
  orFilter->SetInput2(reader2->GetOutput());

  maskThresholdFilter->SetInput(orFilter->GetOutput());
  maskThresholdFilter->SetOutsideValue( 0 );
  maskThresholdFilter->SetInsideValue( 1 );
  maskThresholdFilter->SetLowerThreshold( 1 ) ;
  maskThresholdFilter->SetUpperThreshold( 10 );


  //Write output  
  writer->SetInput( maskThresholdFilter->GetOutput() );
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
