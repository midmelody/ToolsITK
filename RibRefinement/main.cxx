#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
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
      std::cerr << argv[0] << " inputImageFile outputImageFile minObjectSize maxObjectSize" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ConnectedComponentImageFilter< ImageType, ImageType > ConnectedComponentFilterType;
  typedef itk::RelabelComponentImageFilter< ImageType, ImageType >  RelabelFilterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType >  ThresholdFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ConnectedComponentFilterType::Pointer connectComponentFilter = ConnectedComponentFilterType::New();
  RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  ConnectedComponentFilterType::Pointer connectComponentFilter2 = ConnectedComponentFilterType::New();
  RelabelFilterType::Pointer relabelFilter2 = RelabelFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
 
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
  
  connectComponentFilter->SetInput( reader->GetOutput() );
  connectComponentFilter->SetFullyConnected(1);
  relabelFilter->SetInput( connectComponentFilter->GetOutput() );
  relabelFilter->SetMinimumObjectSize(atoi(argv[3]));
  relabelFilter->Update();

  typedef std::vector< itk::SizeValueType > SizesInPixelsType;
  const SizesInPixelsType &  sizesInPixels = relabelFilter->GetSizeOfObjectsInPixels();
  SizesInPixelsType::const_iterator sizeItr = sizesInPixels.begin();
  SizesInPixelsType::const_iterator sizeEnd = sizesInPixels.end();
  std::cout << "Number of pixels per class " << std::endl;
  unsigned int kclass = 1;
  int threLow = 1;
  while( sizeItr != sizeEnd )
    {
      std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
      if(*sizeItr>atoi(argv[4]))
	threLow++;
      ++kclass;
      ++sizeItr;
    }
  int threHigh = kclass;
  thresholdFilter->SetInput( relabelFilter->GetOutput());
  const PixelType outsideValue = 0;
  const PixelType insideValue  = 1;
  thresholdFilter->SetOutsideValue( outsideValue );
  thresholdFilter->SetInsideValue(  insideValue  );
  PixelType lowerThreshold = threLow;
  PixelType upperThreshold = threHigh-1;
  if( lowerThreshold > upperThreshold )
    lowerThreshold = upperThreshold;

  std::cout<<"Low: "<<lowerThreshold<<" High: "<<upperThreshold<<std::endl;

  thresholdFilter->SetLowerThreshold( lowerThreshold );
  thresholdFilter->SetUpperThreshold( upperThreshold );

  connectComponentFilter2->SetInput( thresholdFilter->GetOutput() );
  connectComponentFilter2->SetFullyConnected(1);
  relabelFilter2->SetInput( connectComponentFilter2->GetOutput() );
  relabelFilter2->SetMinimumObjectSize(atoi(argv[3]));

  //Write output  
  writer->SetInput( relabelFilter2->GetOutput() );
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
