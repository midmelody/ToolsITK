#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMultiplyImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
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
      std::cerr << argv[0] << " inputImageFile1 inputImageFile2  outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef double PixelType;
  typedef int OutPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< OutPixelType, Dimension > OutImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType > InRescaleIntensityImageFilterType;
  typedef itk::MultiplyImageFilter <ImageType, ImageType > MultiplyImageFilterType;
  typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType > RescaleIntensityImageFilterType;
  typedef itk::ImageFileWriter< OutImageType > WriterType;

  //Filters
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New();
  RescaleIntensityImageFilterType::Pointer rescaler = RescaleIntensityImageFilterType::New();
  InRescaleIntensityImageFilterType::Pointer rescalerIn1 = InRescaleIntensityImageFilterType::New();
  InRescaleIntensityImageFilterType::Pointer rescalerIn2 = InRescaleIntensityImageFilterType::New();
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
  rescalerIn1->SetInput(reader1->GetOutput());
  rescalerIn1->SetOutputMaximum(10000);
  rescalerIn1->SetOutputMinimum(0);
  rescalerIn2->SetInput(reader2->GetOutput());
  rescalerIn2->SetOutputMaximum(10000);
  rescalerIn2->SetOutputMinimum(0);

  multiplyFilter->SetInput1(rescalerIn1->GetOutput());
  multiplyFilter->SetInput2(rescalerIn2->GetOutput());

  //Write output  
  rescaler->SetInput( multiplyFilter->GetOutput() );
  rescaler->SetOutputMaximum(10000);
  rescaler->SetOutputMinimum(0);
  writer->SetInput( rescaler->GetOutput() );
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
