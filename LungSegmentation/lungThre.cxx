#include "itkImage.h"
#include "itkImageFileReader.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"

#include "itkOtsuThresholdImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile airwayImageFile outputImageFile thrshold" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::BinaryBallStructuringElement< PixelType, Dimension > StructuringElementType;
  typedef itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementType> ErodeFilterType;
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType> DilateFilterType;
  typedef itk::OtsuThresholdImageFilter< ImageType, ImageType > ThresholdFilterType;
  typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType > SubtractImageFilter;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > BinaryThresholdImageFilter;
  typedef itk::ConnectedComponentImageFilter< ImageType, ImageType > ConnectedComponentFilterType;
  typedef itk::RelabelComponentImageFilter< ImageType, ImageType >  RelabelFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerAirway = ReaderType::New();
  DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
  ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  SubtractImageFilter::Pointer subtractFilter = SubtractImageFilter::New();
  BinaryThresholdImageFilter::Pointer binaryThresholdFilter = BinaryThresholdImageFilter::New();
  BinaryThresholdImageFilter::Pointer binaryThresholdFilterF = BinaryThresholdImageFilter::New();
  ConnectedComponentFilterType::Pointer connectComponentFilter = ConnectedComponentFilterType::New();
  RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  WriterType::Pointer writer = WriterType::New();
  WriterType::Pointer writerT = WriterType::New();
  //Parameters
  reader->SetFileName( argv[1] );
  readerAirway->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
 
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

  int thre =  atoi(argv[4]);  
  //std::cout<<thre<<std::endl;
  binaryThresholdFilterF->SetInput( reader->GetOutput() );
  binaryThresholdFilterF->SetLowerThreshold( 0 );
  binaryThresholdFilterF->SetUpperThreshold( thre );
  binaryThresholdFilterF->SetOutsideValue( 0 );
  binaryThresholdFilterF->SetInsideValue( 1 );

  //writerT->SetFileName( "threTemp.hdr" );
  // writerT->SetInput(binaryThresholdFilterF->GetOutput());
  //writerT->Update();

  StructuringElementType  dilateSE;
  dilateSE.SetRadius( 5 );
  dilateSE.CreateStructuringElement();
  dilateFilter->SetKernel( dilateSE );
  dilateFilter->SetInput(readerAirway->GetOutput());
  dilateFilter->SetDilateValue( 1 );

  subtractFilter->SetInput1( binaryThresholdFilterF->GetOutput() );
  subtractFilter->SetInput2( dilateFilter->GetOutput() );

  binaryThresholdFilter->SetInput( subtractFilter->GetOutput() );
  binaryThresholdFilter->SetOutsideValue( 0 );
  binaryThresholdFilter->SetInsideValue( 1 );
  binaryThresholdFilter->SetLowerThreshold( 1 ) ;
  binaryThresholdFilter->SetUpperThreshold( 10 );
  /*
  StructuringElementType  erodeSE;
  erodeSE.SetRadius( 5 );
  erodeSE.CreateStructuringElement();
  erodeFilter->SetKernel( erodeSE );
  erodeFilter->SetErodeValue( 1 );
  erodeFilter->SetInput(binaryThresholdFilter->GetOutput());
  */
  connectComponentFilter->SetInput( binaryThresholdFilter->GetOutput() );
  connectComponentFilter->SetFullyConnected(0);

  relabelFilter->SetInput( connectComponentFilter->GetOutput() );
  relabelFilter->SetMinimumObjectSize(100000);

  //Write output  
  writer->SetInput( relabelFilter->GetOutput() );
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


  typedef std::vector< itk::SizeValueType > SizesInPixelsType;
  const SizesInPixelsType &  sizesInPixels = relabelFilter->GetSizeOfObjectsInPixels();
  SizesInPixelsType::const_iterator sizeItr = sizesInPixels.begin();
  SizesInPixelsType::const_iterator sizeEnd = sizesInPixels.end();
  std::cout << "Number of pixels per class " << std::endl;
  unsigned int kclass = 0;
  while( sizeItr != sizeEnd )
    {
      std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
      ++kclass;
      ++sizeItr;
    }








  return EXIT_SUCCESS;
}
