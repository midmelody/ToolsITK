#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile1 inputImageFile2 label1 label2 outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;
  typedef itk::MaximumImageFilter< ImageType, ImageType, ImageType > MaximumFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  RescaleFilterType::Pointer rescaleFilter1 = RescaleFilterType::New();
  RescaleFilterType::Pointer rescaleFilter2 = RescaleFilterType::New();
  MaximumFilterType::Pointer maximumFilter = MaximumFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );
  int label1 = atoi( argv[3] );
  int label2 = atoi( argv[4] );
  writer->SetFileName( argv[5] );
 
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
  
  rescaleFilter1->SetInput(reader1->GetOutput());
  rescaleFilter1->SetOutputMinimum(0);
  rescaleFilter1->SetOutputMaximum(label1);

  rescaleFilter2->SetInput(reader2->GetOutput());
  rescaleFilter2->SetOutputMinimum(0);
  rescaleFilter2->SetOutputMaximum(label2);

  maximumFilter->SetInput1(rescaleFilter1->GetOutput());
  maximumFilter->SetInput2(rescaleFilter2->GetOutput());

  //Write output  
  writer->SetInput( maximumFilter->GetOutput() );
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
