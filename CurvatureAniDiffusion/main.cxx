#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

int main( int argc, char * argv[] )
{
  if( argc < 6 ) 
    { 
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile  outputImageFile ";
    std::cerr << "numberOfIterations  timeStep  conductance" << std::endl;
    return EXIT_FAILURE;
    }
  
  typedef float    PixelType;
  typedef int OutputPixelType;
  typedef itk::Image< PixelType, 3 >   ImageType;
  typedef itk::Image< OutputPixelType, 3 > OutputImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::CurvatureAnisotropicDiffusionImageFilter< ImageType, ImageType >  FilterType;
  typedef itk::CastImageFilter< ImageType, OutputImageType > CastFilterType;


  FilterType::Pointer filter = FilterType::New();
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  filter->SetInput( reader->GetOutput() );

  const unsigned int numberOfIterations = atoi( argv[3] );
  const double       timeStep = atof( argv[4] );
  const double       conductance = atof( argv[5] );

  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );
  filter->Update();

  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(filter->GetOutput());

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput(  castFilter->GetOutput());
  writer->Update();

  return EXIT_SUCCESS;
}

