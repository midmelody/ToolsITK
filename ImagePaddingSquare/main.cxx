#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkConstantPadImageFilter.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
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

  reader->Update();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = size[0];
  int pCol = size[1];
  int diff;
  ImageType::SizeType extendRegion;
  extendRegion[2] = 0;
  if(pRow<pCol)
    {
      diff = pCol - pRow;
      extendRegion[0] = diff/2;
      extendRegion[1] = 0;
    }
  else
    {
      diff = pRow - pCol;
      extendRegion[0] = 0;      
      extendRegion[1] = diff/2;
    }

  ImageType::PixelType constantPixel = 0;
  padFilter->SetInput(reader->GetOutput());
  padFilter->SetPadLowerBound(extendRegion);
  padFilter->SetPadUpperBound(extendRegion);
  padFilter->SetConstant(constantPixel);
  writer->SetInput( padFilter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
