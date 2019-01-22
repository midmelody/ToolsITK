#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkContourMeanDistanceImageFilter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputMaskFile1 inputMaskFile2 " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ContourMeanDistanceImageFilter< ImageType, ImageType > DistFilterType;

  //Filters
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  DistFilterType::Pointer distFilter = DistFilterType::New();

  //Parameters
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );

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
  
  distFilter->SetInput1(reader1->GetOutput());
  distFilter->SetInput2(reader2->GetOutput());
  distFilter->SetUseImageSpacing(true);
  distFilter->Update();
  
  std::cout<<"HausdorffDistance: "<<distFilter->GetMeanDistance()<<std::endl;

  return EXIT_SUCCESS;
}
