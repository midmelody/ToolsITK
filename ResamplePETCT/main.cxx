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
#include "stdio.h" 

//The program resamples image to target spacing
int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " PETImageFile regPETImageFile CTImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef double PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ResampleImageFilter< ImageType, ImageType >  ResampleFilterType;
  typedef itk::FlipImageFilter <ImageType> FlipImageFilterType;
  typedef itk::IdentityTransform< double, Dimension >  TransformType;
  typedef itk::BSplineInterpolateImageFunction< ImageType, double >  InterpolatorType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer readerPET = ReaderType::New();
  ReaderType::Pointer readerCT = ReaderType::New();
  ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  TransformType::Pointer transform = TransformType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  readerPET->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  readerCT->SetFileName( argv[3] );
  resampleFilter->SetTransform( transform );
  resampleFilter->SetInterpolator( interpolator );

  //Pipeline
  try
    {
      readerPET->Update();
      readerCT->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  //Get image specs
  ImageType::SpacingType spacing = readerPET->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = readerPET->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = readerPET->GetOutput()->GetDirection();
  ImageType::SizeType  sizePET = readerPET->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::SizeType  sizeCT = readerCT->GetOutput()->GetRequestedRegion().GetSize();
  std::cout<<"sizeCT: "<<sizeCT[0]<<" "<<sizeCT[1]<<" "<<sizeCT[2]<<std::endl;
  std::cout<<"sizePET: "<<sizePET[0]<<" "<<sizePET[1]<<" "<<sizePET[2]<<std::endl;
  spacing[0] = spacing[0]*sizePET[0]/sizeCT[0];
  spacing[1] = spacing[1]*sizePET[1]/sizeCT[1];
  spacing[2] = spacing[2]*sizePET[2]/sizeCT[2];
  resampleFilter->SetOutputOrigin( origin );
  resampleFilter->SetOutputSpacing( spacing );
  resampleFilter->SetOutputDirection( direction );
  resampleFilter->SetSize( sizeCT );
  resampleFilter->SetInput( readerPET->GetOutput() );
 
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
