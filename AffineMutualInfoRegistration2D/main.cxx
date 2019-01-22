#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkImageRegistrationMethod.h"
#include "itkAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"

class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate(): m_CumulativeIterationIndex(0) {};
public:
  typedef  itk::RegularStepGradientDescentOptimizer  OptimizerType;
  typedef  const OptimizerType *                     OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer =
      dynamic_cast< OptimizerPointer >( object );
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
	return;
      }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << "  " <<
      m_CumulativeIterationIndex++ << std::endl;
  }
private:
  unsigned int m_CumulativeIterationIndex;
};

template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );
protected:
  RegistrationInterfaceCommand() {};
public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
  typedef   OptimizerType *                            OptimizerPointer;
  void Execute(itk::Object * object, const itk::EventObject & event)
  {
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
	return;
      }
    RegistrationPointer registration =
                        dynamic_cast<RegistrationPointer>( object );
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >(
                       registration->GetOptimizer() );

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : "
              << registration->GetCurrentLevel()  << std::endl;
    std::cout << std::endl;

    if ( registration->GetCurrentLevel() == 0 )
      {
      optimizer->SetMaximumStepLength( 4.00 );
      optimizer->SetMinimumStepLength(  0.05 );
      }
    else
      {
      optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / 4.0 );
      optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 10.0 );
      }
  }
  void Execute(const itk::Object * , const itk::EventObject & )
    { return; }
};

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
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::LinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
  typedef itk::MattesMutualInformationImageToImageMetric< InternalImageType, InternalImageType > MetricType;
  typedef OptimizerType::ScalesType       OptimizerScalesType;
  typedef itk::MultiResolutionImageRegistrationMethod< InternalImageType, InternalImageType > RegistrationType;
  typedef itk::RecursiveMultiResolutionPyramidImageFilter< InternalImageType, InternalImageType > FixedImagePyramidType;
  typedef itk::RecursiveMultiResolutionPyramidImageFilter< InternalImageType, InternalImageType > MovingImagePyramidType;

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
  registration->SetFixedImage( smootherF->GetOutput() );
  registration->SetMovingImage( smootherM->GetOutput() );
  smootherF->Update();
  registration->SetFixedImageRegion( smootherF->GetOutput()->GetBufferedRegion() );

  typedef itk::CenteredTransformInitializer< TransformType, ImageType, ImageType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  RGBconverterF->GetOutput() );
  initializer->SetMovingImage( RGBconverterM->GetOutput() );
  initializer->MomentsOn();
  initializer->InitializeTransform();
  registration->SetInitialTransformParameters( transform->GetParameters() );

  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = 1.0;
  optimizerScales[4] = 1.0 / 1e7;
  optimizerScales[5] = 1.0 / 1e7;

  optimizer->SetScales( optimizerScales );

  metric->SetNumberOfHistogramBins( 128 );
  metric->SetNumberOfSpatialSamples( 5000000 );
  metric->ReinitializeSeed( 76926294 );
  optimizer->SetNumberOfIterations(  200  );
  optimizer->SetRelaxationFactor( 0.8 );

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
  CommandType::Pointer command = CommandType::New();
  registration->AddObserver( itk::IterationEvent(), command );
  registration->SetNumberOfLevels( 3 );

  registration->StartRegistration();

  std::cout << "Optimizer Stopping Condition = "
            << optimizer->GetStopCondition() << std::endl;

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType finalParameters = registration->GetLastTransformParameters();
  unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  double bestValue = optimizer->GetValue();

  // Print out results
  std::cout << std::endl << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " versor X      = " << finalParameters[0] << std::endl;
  std::cout << " versor Y      = " << finalParameters[4] << std::endl;
  std::cout << " versor Z      = " << finalParameters[8]  << std::endl;
  std::cout << " Translation X = " << finalParameters[9]  << std::endl;
  std::cout << " Translation Y = " << finalParameters[10]  << std::endl;
  std::cout << " Translation Z = " << finalParameters[11]  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

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
