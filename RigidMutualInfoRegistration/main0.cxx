#include "itkImageRegistrationMethodv4.h"
#include "itkTranslationTransform.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkCommand.h"

#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:
  typedef itk::RegularStepGradientDescentOptimizerv4<double> OptimizerType;
  typedef   const OptimizerType *                            OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
    OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << argv[0];
      std::cerr << " fixedImageFile  movingImageFile outputImagefile " << std::endl;
      return EXIT_FAILURE;
    }

  const    unsigned int    Dimension = 3;
  typedef  float   PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;

  typedef itk::TranslationTransform< double, Dimension >         TransformType;
  typedef itk::RegularStepGradientDescentOptimizerv4<double>     OptimizerType;
  typedef itk::ImageRegistrationMethodv4< ImageType, ImageType, TransformType > RegistrationType;
  typedef itk::MattesMutualInformationImageToImageMetricv4< ImageType, ImageType > MetricType;
 
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  MetricType::Pointer metric = MetricType::New();
  registration->SetMetric( metric  );
  registration->SetOptimizer(     optimizer     );
  unsigned int numberOfBins = 24;
  metric->SetNumberOfHistogramBins( numberOfBins );
  metric->SetUseMovingImageGradientFilter( false );
  metric->SetUseFixedImageGradientFilter( false );

  
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer fixedImageReader  = ImageReaderType::New();
  ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );
  fixedImageReader->Update();
  movingImageReader->Update();

  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );

  optimizer->SetLearningRate( 8.00 );
  optimizer->SetMinimumStepLength( 0.001 );
  optimizer->SetNumberOfIterations( 200 );
  optimizer->ReturnBestParametersAndValueOn();

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

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

  RegistrationType::MetricSamplingStrategyType  samplingStrategy  =  RegistrationType::RANDOM;
  double samplingPercentage = 0.80;
  registration->SetMetricSamplingStrategy( samplingStrategy );
  registration->SetMetricSamplingPercentage( samplingPercentage );
  registration->MetricSamplingReinitializeSeed( 121213 );

  registration->Update();

  TransformType::ParametersType finalParameters = registration->GetOutput()->Get()->GetParameters();

  double TranslationAlongX = finalParameters[0];
  double TranslationAlongY = finalParameters[1];

  unsigned int numberOfIterations = optimizer->GetCurrentIteration();

  double bestValue = optimizer->GetValue();

  //
  std::cout << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " Translation X = " << TranslationAlongX  << std::endl;
  std::cout << " Translation Y = " << TranslationAlongY  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;
  std::cout << " Stop Condition  = " << optimizer->GetStopConditionDescription() << std::endl;

  typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType;

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( registration->GetTransform() );
  resample->SetInput( movingImageReader->GetOutput() );

  ImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  PixelType defaultPixelValue = 100;

  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( defaultPixelValue );

  typedef  unsigned char  OutputPixelType;

  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter< ImageType, ImageType > CastFilterType;

  typedef itk::ImageFileWriter< ImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[3] );

  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();
  
  return EXIT_SUCCESS;
}
