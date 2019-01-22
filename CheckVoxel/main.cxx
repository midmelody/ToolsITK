#include "itkImage.h"
#include "itkImageFileReader.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program counts the voxel in a binary image for volume

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile X Y Z" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef double PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();

  //Parameters
  reader->SetFileName( argv[1] );
 
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
  
  ImageType::IndexType index;
  index[0] = atoi(argv[2]);
  index[1] = atoi(argv[3]);
  index[2] = atoi(argv[4]);


  ImageType::PixelType pixel;
  pixel = reader->GetOutput()->GetPixel( index );

  std::cout<<"Intensity: "<< pixel <<std::endl;

  return EXIT_SUCCESS;
}
