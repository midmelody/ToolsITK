#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkFlipImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"

//The program resamples image to target spacing
int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile Axis(0/1/2)" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::FlipImageFilter<ImageType> FlipImageFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

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
  
  //Flip
  itk::FixedArray<bool, 3> flipAxes;
  flipAxes[0] = false;
  flipAxes[1] = false;
  flipAxes[2] = false;
  int axis = atoi(argv[3]);
  if(axis==0)
    flipAxes[0] = true;
  else if(axis==1)
    flipAxes[1] = true;
  else if(axis==2)
    flipAxes[2] = true;
  else
    return EXIT_FAILURE;

  flipFilter->SetInput(reader->GetOutput());
  flipFilter->SetFlipAxes(flipAxes);

  //Write output  
  writer->SetInput( flipFilter->GetOutput() );
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
