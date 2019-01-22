#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile outputImageFile sigma " << std::endl;
    return EXIT_FAILURE;
    }

  typedef float InputPixelType;
  typedef unsigned char OutputPixelType;
  typedef itk::Image< InputPixelType,  3 >   InputImageType;
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;
  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typedef itk::RecursiveGaussianImageFilter< InputImageType, InputImageType >  FilterType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  FilterType::Pointer filterX = FilterType::New();
  FilterType::Pointer filterY = FilterType::New();
  FilterType::Pointer filterZ = FilterType::New();

  filterX->SetDirection( 0 );   // 0 --> X direction
  filterY->SetDirection( 1 );   // 1 --> Y direction
  filterZ->SetDirection( 2 );   // 2 --> Z direction

  filterX->SetOrder( FilterType::ZeroOrder );
  filterY->SetOrder( FilterType::ZeroOrder );
  filterZ->SetOrder( FilterType::ZeroOrder );

  filterX->SetNormalizeAcrossScale( false );
  filterY->SetNormalizeAcrossScale( false );
  filterZ->SetNormalizeAcrossScale( false );

  filterX->SetInput( reader->GetOutput() );
  filterY->SetInput( filterX->GetOutput() );
  filterZ->SetInput( filterY->GetOutput() );

  const double sigma = atof( argv[3] );

  filterX->SetSigma( sigma );
  filterY->SetSigma( sigma );
  filterZ->SetSigma( sigma );

  filterZ->Update();
 
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  rescaler->SetInput( filterY->GetOutput() );
  writer->SetInput( rescaler->GetOutput() );
  writer->Update();


  return EXIT_SUCCESS;
}

