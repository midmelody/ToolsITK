#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program resamples image to target spacing
int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile targetSizeX targetSizeY targetSizeZ" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ResampleImageFilter< ImageType, ImageType >  ResampleFilterType;
  typedef itk::FlipImageFilter <ImageType> FlipImageFilterType;
  typedef itk::IdentityTransform< double, Dimension >  TransformType;
  typedef itk::BSplineInterpolateImageFunction< ImageType, double >  InterpolatorType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  TransformType::Pointer transform = TransformType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  resampleFilter->SetTransform( transform );
  resampleFilter->SetInterpolator( interpolator );

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
  float targetSizeX = atof(argv[3]);
  float targetSizeY = atof(argv[4]);
  float targetSizeZ = atof(argv[5]);
  //spacing[0] = size[0] * spacing[0] / targetSizeX; 
  //spacing[1] = size[1] * spacing[1] / targetSizeY; 
  //spacing[2] = size[2] * spacing[2] / targetSizeZ;
  size[0] = targetSizeX;
  size[1] = targetSizeY;
  size[2] = targetSizeZ;
  resampleFilter->SetOutputOrigin( origin );
  resampleFilter->SetOutputSpacing( spacing );
  resampleFilter->SetOutputDirection( direction );
  resampleFilter->SetSize( size );
  resampleFilter->SetInput( reader->GetOutput() );
  std::cout<<spacing<<" "<<size<<" "<<origin<<std::endl<<direction<<std::endl;
 
  //Flip
  itk::FixedArray<bool, 3> flipAxes;
  flipAxes[0] = false;
  flipAxes[1] = false;
  flipAxes[2] = true;
  flipFilter->SetInput(resampleFilter->GetOutput());
  flipFilter->SetFlipAxes(flipAxes);

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
