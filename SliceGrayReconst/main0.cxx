#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::SliceBySliceImageFilter< ImageType, ImageType > SliceFilterType;
  typedef SliceFilterType::InternalInputImageType SliceImageType;
  typedef itk::BinaryBallStructuringElement< PixelType, 2 > StructuringElementType;
  typedef itk::GrayscaleMorphologicalClosingImageFilter< SliceImageType, SliceImageType, StructuringElementType >  ClosingFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  SliceFilterType::Pointer slicer = SliceFilterType::New();
  ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  StructuringElementType  structuringElement;
  structuringElement.SetRadius( 1 );
  structuringElement.CreateStructuringElement();
  closingFilter->SetKernel( structuringElement );

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
  
  slicer->SetFilter( closingFilter );
  slicer->SetInput( reader->GetOutput() );

  //Write output  
  writer->SetInput( slicer->GetOutput() );
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
