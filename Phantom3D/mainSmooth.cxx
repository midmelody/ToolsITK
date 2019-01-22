#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << "inputImage outputImage smoothSigma" << std::endl;
      return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef unsigned char OutPixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef itk::Image< OutPixelType, Dimension >  OutImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType > GaussianFilterType;
  typedef itk::RescaleIntensityImageFilter< ImageType, OutImageType> RescaleFilterType;
  typedef itk::ImageFileWriter< OutImageType > WriterType;


  //Read image 
  ReaderType::Pointer Reader = ReaderType::New();
  Reader->SetFileName( argv[1] );
  Reader->Update();
 
  float smoothRatio = atof(argv[3]);
  GaussianFilterType::Pointer GaussX = GaussianFilterType::New();
  GaussianFilterType::Pointer GaussY = GaussianFilterType::New();
  GaussianFilterType::Pointer GaussZ = GaussianFilterType::New();
  GaussX->SetDirection( 0 );
  GaussY->SetDirection( 1 );
  GaussZ->SetDirection( 2 );
  GaussX->SetSigma( Reader->GetOutput()->GetSpacing()[0] * smoothRatio);
  GaussY->SetSigma( Reader->GetOutput()->GetSpacing()[1] * smoothRatio);
  GaussZ->SetSigma( Reader->GetOutput()->GetSpacing()[2] * smoothRatio);
  GaussX->SetZeroOrder();
  GaussY->SetZeroOrder();
  GaussZ->SetZeroOrder();
  GaussX->SetInput( Reader->GetOutput() );
  GaussY->SetInput( GaussX->GetOutput() );
  GaussZ->SetInput( GaussY->GetOutput() );
  GaussZ->Update();

  //Rescale

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(GaussZ->GetOutput());
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);

  //Output
  WriterType::Pointer  Writer  = WriterType::New();
  Writer->SetInput(rescaleFilter->GetOutput() );
  Writer->SetFileName(argv[2]);
  Writer->Update();
  return EXIT_SUCCESS;
}
