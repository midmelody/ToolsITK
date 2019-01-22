#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkConstantPadImageFilter.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile padSize" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ConstantPadImageFilter <ImageType, ImageType> ConstantPadImageFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ConstantPadImageFilterType::Pointer padFilter = ConstantPadImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  ImageType::SizeType extendRegion;
  extendRegion[0] = atoi( argv[3] );
  extendRegion[1] = atoi( argv[3] );
  extendRegion[2] = atoi( argv[3] );
  ImageType::PixelType constantPixel = 0;
  padFilter->SetInput(reader->GetOutput());
  padFilter->SetPadLowerBound(extendRegion);
  padFilter->SetPadUpperBound(extendRegion);
  padFilter->SetConstant(constantPixel);
  writer->SetInput( padFilter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
