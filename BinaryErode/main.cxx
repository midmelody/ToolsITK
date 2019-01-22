#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program erodes a binary image with given radius

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile radius" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::BinaryBallStructuringElement< PixelType, Dimension > StructuringElementType;
  typedef itk::BinaryErodeImageFilter <ImageType, ImageType, StructuringElementType> BinaryErodeImageFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  unsigned int radius = atoi( argv[3] );

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
  
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  erodeFilter->SetInput(reader->GetOutput());
  erodeFilter->SetKernel(structuringElement);
  erodeFilter->SetBackgroundValue(0);
  erodeFilter->SetForegroundValue(1);
 
  //Write output  
  writer->SetInput( erodeFilter->GetOutput() );
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
