#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMaximumImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile1 inputImageFile2 outputImageFile " << std::endl;
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
  MaxFilterType::Pointer maxFilter = MaxFilterType::New();
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
  
  maxFilter->SetInput1( reader1->GetOutput() );
  maxFilter->SetInput2( reader2->GetOutput() );
 
  writer->SetInput( maxFilter->GetOutput() );

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
