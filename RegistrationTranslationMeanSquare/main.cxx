#include "itkImageRegistrationMethodv4.h"
#include "itkTranslationTransform.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << " outputImagefile" << std::endl;
    return EXIT_FAILURE;
    }
  
  const unsigned int                          Dimension = 3;
  typedef  float                              PixelType;
  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;
  typedef itk::TranslationTransform< double, Dimension > TransformType;
  typedef itk::RegularStepGradientDescentOptimizerv4<double> OptimizerType;
  typedef itk::MeanSquaresImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;
  typedef itk::ImageRegistrationMethodv4< FixedImageType, MovingImageType, TransformType > RegistrationType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  TransformType::Pointer initialTransform   = TransformType::New();
  
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInitialTransform(     initialTransform     );
  
  typedef itk::LinearInterpolateImageFunction< FixedImageType, double > FixedLinearInterpolatorType;
  typedef itk::LinearInterpolateImageFunction< MovingImageType, double > MovingLinearInterpolatorType;
  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );  
  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );

  FixedLinearInterpolatorType::Pointer fixedInterpolator =   FixedLinearInterpolatorType::New();
  MovingLinearInterpolatorType::Pointer movingInterpolator = MovingLinearInterpolatorType::New();

  metric->SetFixedInterpolator(  fixedInterpolator  );
  metric->SetMovingInterpolator(  movingInterpolator  );

  TransformType::Pointer movingInitialTransform = TransformType::New();

  TransformType::ParametersType initialParameters( movingInitialTransform->GetNumberOfParameters() );
  initialParameters[0] = 0.0;
  initialParameters[1] = 0.0;
  initialParameters[2] = 0.0; 

  movingInitialTransform->SetParameters( initialParameters );

  registration->SetMovingInitialTransform( movingInitialTransform );

  TransformType::Pointer   identityTransform = TransformType::New();
  identityTransform->SetIdentity();

  registration->SetFixedInitialTransform( identityTransform );

  optimizer->SetLearningRate( 0.1 );
  optimizer->SetMinimumStepLength( 0.001 );
  optimizer->SetRelaxationFactor( 0.5 );
  optimizer->SetNumberOfIterations( 200 );

  const unsigned int numberOfLevels = 1;

  RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize( 1 );
  shrinkFactorsPerLevel[0] = 1;

  RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize( 1 );
  smoothingSigmasPerLevel[0] = 0;

  registration->SetNumberOfLevels ( numberOfLevels );
  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

  try
    {
      registration->Update();
      //std::cout << "Optimizer stop condition: "
      //	<< registration->GetOptimizer()->GetStopConditionDescription()
      //	<< std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }

  TransformType::ConstPointer transform = dynamic_cast< const TransformType *> ( registration->GetTransform() );

  TransformType::ParametersType finalParameters = transform->GetParameters();
  const double TranslationAlongX = finalParameters[0];
  const double TranslationAlongY = finalParameters[1];
  const double TranslationAlongZ = finalParameters[2];
  
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

  const double bestValue = optimizer->GetValue();

  //std::cout << "Result = " << std::endl;
  std::cout << " TranslationX " << TranslationAlongX  << "; ";
  std::cout << " TranslationY " << TranslationAlongY  << "; ";
  std::cout << " TranslationZ " << TranslationAlongZ  << std::endl; 
  //std::cout << " Iterations    = " << numberOfIterations << std::endl;
  //std::cout << " Metric value  = " << bestValue          << std::endl;

  typedef itk::CompositeTransform< double, Dimension > CompositeTransformType;
  CompositeTransformType::Pointer outputCompositeTransform = CompositeTransformType::New();
  outputCompositeTransform->AddTransform( movingInitialTransform );
  outputCompositeTransform->AddTransform( registration->GetModifiableTransform() );

  typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >    ResampleFilterType;

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( movingImageReader->GetOutput() );
  resampler->SetTransform( outputCompositeTransform );

  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 100 );

  typedef int                            OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[3] );

  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();
  
  return EXIT_SUCCESS;
}
