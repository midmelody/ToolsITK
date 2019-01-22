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
  if( argc < 5 )
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
  typedef itk::FlipImageFilter<ImageType> FlipImageFilterType;
  typedef itk::AffineTransform< double, Dimension >  TransformType;
  typedef itk::LinearInterpolateImageFunction< ImageType, double >  InterpolatorType;
  typedef itk::ResampleImageFilter<ImageType, ImageType>  FilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
  FilterType::Pointer filter = FilterType::New();
  TransformType::Pointer transform = TransformType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( 0 );
  
  //Pipeline
  reader->Update();
  
  //Flip version
  FlipImageFilterType::FlipAxesArrayType flipAxes;
  flipAxes[0] = true;
  flipAxes[1] = false;
  flipFilter->SetInput(reader->GetOutput());
  flipFilter->SetFlipAxes(flipAxes);
  flipFilter->Update();

  int i,j;
  int degree[15] = {30,45,60,90,120,135,150,180,210,225,240,270,300,315,330};
  ImageType::Pointer inputImage = ImageType::New();
  ImageType::Pointer outputImage = ImageType::New();
  char filename[100];    
  for(i=0;i<2;i++)
    {
      if(i==0)
	inputImage = reader->GetOutput();
      else
	inputImage = flipFilter->GetOutput();
      
      ImageType::SpacingType spacing = inputImage->GetSpacing();
      ImageType::PointType origin = inputImage->GetOrigin();
      ImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
      filter->SetOutputOrigin( origin );
      filter->SetOutputSpacing( spacing );
      filter->SetOutputDirection( inputImage->GetDirection() );
      filter->SetSize( size );
      double imageCenterX = origin[0] + spacing[0] * size[0] / 2.0;
      double imageCenterY = origin[1] + spacing[1] * size[1] / 2.0;
      
      for(j=0;j<15;j++)
	{
	  //std::cout<<degree[j]<<""<<imageCenterX<<""<<imageCenterY<<std::endl;
	  filter->SetInput(inputImage);
	  TransformType::OutputVectorType translation1;
	  translation1[0] =   -imageCenterX;
	  translation1[1] =   -imageCenterY;
	  transform->Translate( translation1 );
	  const double degreesToRadians = std::atan(1.0) / 45.0;
	  const double angle = degree[j] * degreesToRadians;
	  transform->Rotate2D( -angle, false );
	  TransformType::OutputVectorType translation2;
	  translation2[0] =   imageCenterX;
	  translation2[1] =   imageCenterY;
	  transform->Translate( translation2, false );
	  filter->SetTransform( transform );
	  outputImage = filter->GetOutput();

	  sprintf(filename, "%s_%d_%d.hdr", argv[1], i, j);
          writer->SetFileName( filename );
	  writer->SetInput( outputImage );
	  writer->Update();
	  //outputImage->DisconnectPipeline();
	}
    }
 
  return EXIT_SUCCESS;
}
