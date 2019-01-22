#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile degree" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 2;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::AffineTransform< double, Dimension >  TransformType;
  typedef itk::LinearInterpolateImageFunction< ImageType, double >  InterpolatorType;
  typedef itk::ResampleImageFilter<ImageType, ImageType>  FilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();
  TransformType::Pointer transform = TransformType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  int degree = atoi( argv[3] ); 
  
  //Pipeline
  reader->Update();
  ImageType::Pointer inputImage = ImageType::New();
  inputImage = reader->GetOutput();

  ImageType::SpacingType spacing = inputImage->GetSpacing();
  ImageType::PointType origin = inputImage->GetOrigin();
  ImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
  filter->SetOutputOrigin( origin );
  filter->SetOutputSpacing( spacing );
  filter->SetOutputDirection( inputImage->GetDirection() );
  filter->SetSize( size );
  double imageCenterX = origin[0] + spacing[0] * size[0] / 2.0;
  double imageCenterY = origin[1] + spacing[1] * size[1] / 2.0;
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( 2 );
  filter->SetInput(inputImage);

  TransformType::OutputVectorType translation1;
  translation1[0] =   -imageCenterX;
  translation1[1] =   -imageCenterY;
  transform->Translate( translation1 );
  const double degreesToRadians = std::atan(1.0) / 45.0;
  const double angle = degree * degreesToRadians;
  transform->Rotate2D( -angle, false );
  TransformType::OutputVectorType translation2;
  translation2[0] =   imageCenterX;
  translation2[1] =   imageCenterY;
  transform->Translate( translation2, false );
  filter->SetTransform( transform );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
