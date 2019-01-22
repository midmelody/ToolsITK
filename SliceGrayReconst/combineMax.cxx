#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMaximumImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile1 inputImageFile2 inputImageFile3 outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3; 
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::MaximumImageFilter< ImageType, ImageType, ImageType> MaxFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  ReaderType::Pointer reader3 = ReaderType::New();
  MaxFilterType::Pointer maxFilter1 = MaxFilterType::New();
  MaxFilterType::Pointer maxFilter2 = MaxFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );
  reader3->SetFileName( argv[3] );
  writer->SetFileName( argv[4] );

  //Pipeline
  try
    {
      reader1->Update();
      reader2->Update();
      reader3->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  maxFilter1->SetInput1( reader1->GetOutput() );
  maxFilter1->SetInput2( reader2->GetOutput() );
  maxFilter2->SetInput1( reader3->GetOutput() );
  maxFilter2->SetInput2( maxFilter1->GetOutput() );

  writer->SetInput( maxFilter2->GetOutput() );

  try
    {
      writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout << "ExceptionObject caught !" << std::endl;
      std::cout << err << std::endl;
      return -1;
    }

  return EXIT_SUCCESS;
}
