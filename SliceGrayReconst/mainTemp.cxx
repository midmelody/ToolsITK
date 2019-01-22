#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkOrientImageFilter.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkMaximumImageFilter.h"
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
      std::cerr << argv[0] << " inputImageFile outputImageFile radius direction(1:axial 2:coronal 3:sagittal)" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::OrientImageFilter<ImageType,ImageType> OrientImageFilterType;
  typedef itk::SliceBySliceImageFilter< ImageType, ImageType > SliceFilterType;
  typedef SliceFilterType::InternalInputImageType SliceImageType;
  typedef itk::BinaryBallStructuringElement< PixelType, 2 > StructuringElementType;
  typedef itk::GrayscaleMorphologicalClosingImageFilter< SliceImageType, SliceImageType, StructuringElementType >  ClosingFilterType;
  typedef itk::GrayscaleErodeImageFilter< SliceImageType, SliceImageType, StructuringElementType> ErodeFilterType;
  typedef itk::MaximumImageFilter< SliceImageType, SliceImageType, SliceImageType> maxFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  OrientImageFilterType::Pointer orienterPre = OrientImageFilterType::New();
  OrientImageFilterType::Pointer orienterPost = OrientImageFilterType::New();
  SliceFilterType::Pointer slicer = SliceFilterType::New();
  ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
  SliceFilterType::Pointer slicerSecond = SliceFilterType::New();
  ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  StructuringElementType  structuringElement;
  structuringElement.SetRadius( atoi(argv[3]) );
  structuringElement.CreateStructuringElement();
  closingFilter->SetKernel( structuringElement );
  if(atoi(argv[4])==1)
    orienterPre->SetDesiredCoordinateOrientationToAxial();
  else if(atoi(argv[4])==2)
    orienterPre->SetDesiredCoordinateOrientationToCoronal();
  else if(atoi(argv[4])==3)
    orienterPre->SetDesiredCoordinateOrientationToSagittal();
  else
    {
      std::cerr << "Orientation selection whitin 1, 2, 3 " << std::endl;
      return EXIT_FAILURE; 
    }

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
  
  orienterPre->UseImageDirectionOn();
  orienterPre->SetInput( reader->GetOutput() );

  orienterPre->Update();

  slicer->SetFilter( closingFilter );
  slicer->SetInput( orienterPre->GetOutput() );


  StructuringElementType  structuringElementD;
  structuringElementD.SetRadius( 1 );
  structuringElementD.CreateStructuringElement();
  erodeFilter->SetKernel(structuringElementD);

  slicerSecond->SetInput(slicer->GetOutput());
  slicerSecond->SetFilter(erodeFilter);

  maxFilterType->SetInput1(erodeFilter->GetOutput());
  maxFilterType->SetInput2(slicerSecond->GetInput());








  orienterPost->UseImageDirectionOn();
  orienterPost->SetInput( slicerSecond->GetOutput() );
  orienterPost->SetDesiredCoordinateOrientationToAxial();
  orienterPost->Update();

  //Write output  
  writer->SetInput( orienterPost->GetOutput() );
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
