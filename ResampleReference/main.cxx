#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program resamples image to target spacing
int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " InputImageFile RefImageFile OutputImageFile ForceFlag" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ResampleImageFilter< ImageType, ImageType >  ResampleFilterType;
  typedef itk::IdentityTransform< double, Dimension >  TransformType;
  typedef itk::BSplineInterpolateImageFunction< ImageType, double >  InterpolatorType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer readerIn = ReaderType::New();
  ReaderType::Pointer readerRef = ReaderType::New();
  ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  TransformType::Pointer transform = TransformType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  readerIn->SetFileName( argv[1] );
  readerRef->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
  bool forceFlag = atoi( argv[4] );
  resampleFilter->SetTransform( transform );
  resampleFilter->SetInterpolator( interpolator );

  //Pipeline
  try
    {
      readerIn->Update();
      readerRef->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  //Get image specs
  ImageType::SpacingType spacingRef = readerRef->GetOutput()->GetSpacing();
  ImageType::SpacingType spacingIn = readerIn->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = readerRef->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = readerRef->GetOutput()->GetDirection();
  ImageType::SizeType  sizeIn = readerIn->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::SizeType  sizeRef = readerRef->GetOutput()->GetRequestedRegion().GetSize();  
  ImageType::SizeType  size;
  
  size[0] = sizeIn[0]*spacingIn[0]/spacingRef[0];
  size[1] = sizeIn[1]*spacingIn[1]/spacingRef[1];
  size[2] = sizeIn[2]*spacingIn[2]/spacingRef[2];
  if(forceFlag)
    {
      size[0] = sizeRef[0];
      size[1] = sizeRef[1];
      size[2] = sizeRef[2];
    }

  
  resampleFilter->SetOutputOrigin( origin );
  resampleFilter->SetOutputSpacing( spacingRef );
  resampleFilter->SetOutputDirection( direction );
  resampleFilter->SetSize( size );
  resampleFilter->SetInput( readerIn->GetOutput() );
 
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
