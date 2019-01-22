#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkConstantPadImageFilter.h"
#include "itkCropImageFilter.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile Size_Pad/Trim" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ConstantPadImageFilter <ImageType, ImageType> ConstantPadImageFilterType;
  typedef itk::CropImageFilter <ImageType, ImageType> CropImageFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ConstantPadImageFilterType::Pointer padFilter = ConstantPadImageFilterType::New();
  CropImageFilterType::Pointer cropFilter = CropImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  int size = atoi( argv[3] );
  if(size >0)
    {
      ImageType::SizeType extendRegion;
      extendRegion[0] = size;
      extendRegion[1] = size;
      extendRegion[2] = size;
      ImageType::PixelType constantPixel = 0;
      padFilter->SetInput(reader->GetOutput());
      padFilter->SetPadLowerBound(extendRegion);
      padFilter->SetPadUpperBound(extendRegion);
      padFilter->SetConstant(constantPixel);
      writer->SetInput( padFilter->GetOutput() );
    }
  else
    {
      ImageType::SizeType cropSize;  
      cropSize[0] = -size;
      cropSize[1] = -size;
      cropSize[2] = -size;
      cropFilter->SetInput(reader->GetOutput());
      cropFilter->SetBoundaryCropSize(cropSize);
      writer->SetInput( cropFilter->GetOutput() );
    }

  writer->Update();

  return EXIT_SUCCESS;
}
