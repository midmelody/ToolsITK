#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkImageRegistrationMethod.h"
#include "itkAffineTransform.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkGradientDescentOptimizer.h"
#include "itkNormalizeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"
 

 
int main( int argc, char *argv[] )
{
  //Type Definitions
  const   unsigned int Dimension = 2;

  typedef unsigned char PixelType;
  typedef itk::RGBPixel< PixelType > RGBPixelType;
  typedef float InternalPixelType;

  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< RGBPixelType, Dimension > RGBImageType; 
  typedef itk::Image< float, Dimension> InternalImageType;

  typedef itk::ImageFileReader< RGBImageType > ReaderType;
  typedef itk::RGBToLuminanceImageFilter< RGBImageType, ImageType > RGBConvertFilterType;
  typedef itk::NormalizeImageFilter< ImageType, InternalImageType > NormalizeFilterType;
  typedef itk::DiscreteGaussianImageFilter< InternalImageType,InternalImageType > GaussianFilterType;

  typedef itk::AffineTransform< double, Dimension > TransformType;
  typedef itk::GradientDescentOptimizer OptimizerType;
  typedef itk::LinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
  typedef itk::ImageRegistrationMethod< InternalImageType, InternalImageType > RegistrationType;
  typedef itk::MutualInformationImageToImageMetric< InternalImageType, InternalImageType > MetricType;
  typedef itk::VectorResampleImageFilter< RGBImageType, RGBImageType > VectorResampleFilterType;
  typedef itk::ImageFileWriter< RGBImageType >  WriterType;


  //ITK Filters
  ReaderType::Pointer readerF = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New(); 
  RGBConvertFilterType::Pointer RGBconverterF = RGBConvertFilterType::New();
  RGBConvertFilterType::Pointer RGBconverterM = RGBConvertFilterType::New();
  NormalizeFilterType::Pointer normalizerF = NormalizeFilterType::New();
  NormalizeFilterType::Pointer normalizerM = NormalizeFilterType::New();
  GaussianFilterType::Pointer smootherF  = GaussianFilterType::New();
  GaussianFilterType::Pointer smootherM = GaussianFilterType::New();
  TransformType::Pointer transform = TransformType::New();
  OptimizerType::Pointer optimizer = OptimizerType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  MetricType::Pointer metric = MetricType::New();
  RegistrationType::Pointer registration = RegistrationType::New();
  VectorResampleFilterType::Pointer RGBResampler = VectorResampleFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Filter parameters
  readerF->SetFileName( argv[1] );
  readerM->SetFileName( argv[2] );
  RGBconverterF->SetInput( readerF->GetOutput() );
  RGBconverterM->SetInput( readerM->GetOutput() );
  normalizerF->SetInput( RGBconverterF->GetOutput() );
  normalizerM->SetInput( RGBconverterM->GetOutput() );
  smootherF->SetVariance( 2.0 );
  smootherM->SetVariance( 2.0 );
  smootherF->SetInput( normalizerF->GetOutput() );
  smootherM->SetInput( normalizerM->GetOutput() );
  registration->SetOptimizer( optimizer );
  registration->SetTransform( transform );
  registration->SetInterpolator( interpolator );
  registration->SetMetric( metric );
  metric->SetFixedImageStandardDeviation(  0.4 );
  metric->SetMovingImageStandardDeviation( 0.4 );
  registration->SetFixedImage( smootherF->GetOutput() );
  registration->SetMovingImage( smootherM->GetOutput() );

  normalizerF->Update();
  ImageType::RegionType fixedImageRegion = normalizerF->GetOutput()->GetBufferedRegion();
  registration->SetFixedImageRegion( fixedImageRegion );
  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );
  initialParameters[0] = 1.0;
  initialParameters[1] = 0.0;
  initialParameters[2] = 0.0;
  initialParameters[3] = 1.0;
  initialParameters[4] = 0.0;
  initialParameters[5] = 0.0;
  registration->SetInitialTransformParameters( initialParameters );
  const unsigned int numberOfPixels = fixedImageRegion.GetNumberOfPixels();
  const unsigned int numberOfSamples = static_cast< unsigned int >( numberOfPixels * 0.01 );
  metric->SetNumberOfSpatialSamples( numberOfSamples );
  optimizer->SetLearningRate( 0.1 );
  optimizer->SetNumberOfIterations( 1000 );
  optimizer->MaximizeOn(); 
  registration->Update();
  ParametersType finalParameters = registration->GetLastTransformParameters();
  std::cout << "Final Parameters: " << finalParameters << std::endl;
  unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  double bestValue = optimizer->GetValue();
  std::cout << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;
  std::cout << " Numb. Samples = " << numberOfSamples    << std::endl;
  TransformType::Pointer finalTransform = TransformType::New();
  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );
  
  RGBResampler->SetInput( readerM->GetOutput() );
  RGBResampler->SetSize( readerF->GetOutput()->GetLargestPossibleRegion().GetSize() );
  RGBResampler->SetOutputOrigin(  readerF->GetOutput()->GetOrigin() );
  RGBResampler->SetOutputSpacing( readerF->GetOutput()->GetSpacing() );
  RGBResampler->SetOutputDirection( readerF->GetOutput()->GetDirection() );
  RGBResampler->SetTransform( finalTransform );

  writer->SetInput( RGBResampler->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();
 
  return EXIT_SUCCESS;
}
