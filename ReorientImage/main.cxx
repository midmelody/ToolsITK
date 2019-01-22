#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkOrientImageFilter.h"
#include "itkImageFileWriter.h"

//The program re-orient the image to RAI

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
  typedef double PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::OrientImageFilter<ImageType,ImageType> OrientFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  OrientFilterType::Pointer orienter = OrientFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
 
  //Pipeline
  reader->Update();
  orienter->UseImageDirectionOn();
  orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
  orienter->SetInput( reader->GetOutput() );
  orienter->Update();
  writer->SetInput( orienter->GetOutput() );
  writer->Update();
  
  return EXIT_SUCCESS;
}
