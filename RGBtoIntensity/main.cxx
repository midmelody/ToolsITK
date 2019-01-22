#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{
  //Type Definitions
  const   unsigned int Dimension = 2;
  typedef unsigned char PixelType;
  typedef itk::RGBPixel< PixelType > RGBPixelType;

  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< RGBPixelType, Dimension > RGBImageType; 

  typedef itk::ImageFileReader< RGBImageType > ReaderType;
  typedef itk::RGBToLuminanceImageFilter< RGBImageType, ImageType > RGBConvertFilterType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;


  //ITK Filters
  ReaderType::Pointer reader = ReaderType::New();
  RGBConvertFilterType::Pointer RGBconverter = RGBConvertFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Filter parameters
  reader->SetFileName( argv[1] );
  RGBconverter->SetInput( reader->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->SetInput( RGBconverter->GetOutput() );
  writer->Update();
 
  return EXIT_SUCCESS;
}
