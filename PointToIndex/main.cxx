#include "itkImage.h"
#include "itkImageFileReader.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage x y z";
      std::cerr << " " << std::endl;
      return EXIT_FAILURE;
    }

  typedef int PixelType;
  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::Point< double, ImageType::ImageDimension > PointType;
 
  ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName( argv[1] );
  reader->Update();

  PointType point;
  point[0] = atof( argv[2] );
  point[1] = atof( argv[3] );
  point[2] = atof( argv[4] );
  
  ImageType::IndexType pixelIndex;

  reader->GetOutput()->TransformPhysicalPointToIndex( point, pixelIndex );

  std::cout<<pixelIndex<<std::endl;


  return EXIT_SUCCESS;
}

