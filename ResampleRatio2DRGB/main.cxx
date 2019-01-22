#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkVectorResampleImageFilter.h"
#include "itkIdentityTransform.h"
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
      std::cerr << argv[0] << " inputImageFile outputImageFile ratio" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const   unsigned int Dimension = 2;
  typedef unsigned char PixelType;
  typedef itk::RGBPixel< PixelType > RGBPixelType; 
  typedef itk::Image< RGBPixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::VectorResampleImageFilter< ImageType, ImageType > ResampleFilterType;
  typedef itk::IdentityTransform< double, Dimension >  TransformType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  TransformType::Pointer transform = TransformType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  resampleFilter->SetTransform( transform );
  float ratio = atof(argv[3]);
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
  
  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  std::cout<<spacing<<" "<<size<<" "<<origin<<std::endl<<direction<<std::endl;
  float targetSizeX = size[0] * ratio;
  float targetSizeY = size[1] * ratio;
  size[0] = targetSizeX;
  size[1] = targetSizeY;
  float targetSpacingX = spacing[0] / ratio;
  float targetSpacingY = spacing[1] / ratio;
  spacing[0] = targetSpacingX;
  spacing[1] = targetSpacingY;
  resampleFilter->SetOutputOrigin( origin );
  resampleFilter->SetOutputSpacing( spacing );
  resampleFilter->SetOutputDirection( direction );
  resampleFilter->SetSize( size );
  resampleFilter->SetInput( reader->GetOutput() );

  //Write output  
  writer->SetInput( resampleFilter->GetOutput() );
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
