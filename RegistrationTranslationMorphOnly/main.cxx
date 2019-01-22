#include "itkImage.h"
#include "itkTranslationTransform.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImageFile  outputImageFile ";
    std::cerr << " transX transY transZ" << std::endl;
    return EXIT_FAILURE;
    }
  
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType; 
  typedef itk::TranslationTransform< double, Dimension > TransformType;
  typedef itk::ResampleImageFilter< ImageType, ImageType> ResampleImageFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  TransformType::Pointer transform = TransformType::New();
  TransformType::OutputVectorType translation;
  translation[0] = atof( argv[3] );
  translation[1] = atof( argv[4] );
  translation[2] = atof( argv[5] ); 
  transform->Translate(translation);

  ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
  resampler->SetTransform( transform.GetPointer() );
  resampler->SetInput( reader->GetOutput() );
  resampler->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  reader->GetOutput()->GetOrigin() );
  resampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
  resampler->SetOutputDirection( reader->GetOutput()->GetDirection() );
  resampler->SetDefaultPixelValue( 100 );
  
  WriterType::Pointer writer =  WriterType::New();
  writer->SetFileName( argv[2] );  
  writer->SetInput( resampler->GetOutput() );
  writer->Update();
 
  return EXIT_SUCCESS;
}
