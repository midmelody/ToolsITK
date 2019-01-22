#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkVersorRigid3DTransform.h"

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << argv[0];
      std::cerr << " inputImageFile outputImagefile parameterFile refImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  const    unsigned int    Dimension = 3;
  typedef  float   PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerRef = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  readerRef->SetFileName( argv[4] );

  typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  typedef itk::VersorRigid3DTransform< double > TransformType;
  TransformType::Pointer transform = TransformType::New();

  typedef itk::LinearInterpolateImageFunction< ImageType, double >  InterpolatorType;
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();

  reader->Update();
  readerRef->Update();
  const ImageType * inputImage = reader->GetOutput();
  const ImageType * refImage = readerRef->GetOutput();
  const ImageType::SpacingType & spacing = refImage->GetSpacing();
  const ImageType::PointType & origin  = refImage->GetOrigin();
  const ImageType::SizeType size = refImage->GetLargestPossibleRegion().GetSize();
  
  resampler->SetInput( reader->GetOutput() );
  resampler->SetOutputOrigin( origin );
  resampler->SetOutputSpacing( spacing );
  resampler->SetOutputDirection( inputImage->GetDirection() );
  resampler->SetSize( size );
  resampler->SetInterpolator( interpolator );
  resampler->SetDefaultPixelValue( 0 );
  resampler->SetTransform( transform );
  writer->SetInput( resampler->GetOutput() );

  TransformType::MatrixType matrix;
  TransformType::OffsetType offset;
  FILE *paraFile = fopen(argv[3], "r");
  int i,j;
  float temp;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	if(fscanf(paraFile, "%f\n", &temp) == 1)
	  matrix[i][j] = temp;
      }
  i = 0;
  while (fscanf(paraFile, "%f\n", &temp) == 1)
    {
      offset[i] = temp;
      i++;
    }

  std::cout<<matrix<<std::endl;
  std::cout<<offset<<std::endl;

  fclose(paraFile);
  
  matrix[0][0] = 1;
  matrix[0][1] = 0;
  matrix[0][2] = 0;
  matrix[1][0] = 0;
  matrix[1][1] = 1;
  matrix[1][2] = 0;
  matrix[2][0] = 0;
  matrix[2][1] = 0;
  matrix[2][2] = 1;

  transform->SetMatrix(matrix);
  transform->SetOffset(offset);

  writer->Update();
  return EXIT_SUCCESS;
}
