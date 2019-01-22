#include "itkAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCheckerBoardImageFilter.h"

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
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile outputImagefile ";
    std::cerr << " [hitchImageFile] [hitchOutputImageFile] "<<std::endl;
    return EXIT_FAILURE;
    }

  //ITK Paramters
  const   unsigned int    Dimension = 3;
  typedef int PixelType;
  typedef float InternalPixelType;
  typedef itk::Image< PixelType, Dimension > FixedImageType;
  typedef itk::Image< PixelType, Dimension > MovingImageType;
  typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
  typedef itk::AffineTransform< double, Dimension > TransformType;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::LinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
  typedef itk::MattesMutualInformationImageToImageMetric< InternalImageType, InternalImageType > MetricType;
  typedef OptimizerType::ScalesType OptimizerScalesType;
  typedef itk::MultiResolutionImageRegistrationMethod< InternalImageType, InternalImageType > RegistrationType;
  typedef itk::RecursiveMultiResolutionPyramidImageFilter< InternalImageType, InternalImageType > FixedImagePyramidType;
  typedef itk::RecursiveMultiResolutionPyramidImageFilter< InternalImageType, InternalImageType > MovingImagePyramidType;
  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  typedef itk::CastImageFilter< FixedImageType, InternalImageType > FixedCastFilterType;
  typedef itk::CastImageFilter< MovingImageType, InternalImageType > MovingCastFilterType;
  typedef itk::CenteredTransformInitializer< TransformType, FixedImageType, MovingImageType > TransformInitializerType;
  typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >    ResampleFilterType;
  typedef int OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  //ITK Filter Definitions
  OptimizerType::Pointer         optimizer         = OptimizerType::New();
  InterpolatorType::Pointer      interpolator      = InterpolatorType::New();
  RegistrationType::Pointer      registration      = RegistrationType::New();
  MetricType::Pointer            metric            = MetricType::New();
  TransformType::Pointer         transform         = TransformType::New();
  TransformInitializerType::Pointer initializer    = TransformInitializerType::New();
  FixedImageReaderType::Pointer     fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer    movingImageReader = MovingImageReaderType::New();
  MovingImageReaderType::Pointer  hitchImageReader  = FixedImageReaderType::New();
  FixedCastFilterType::Pointer   fixedCaster       = FixedCastFilterType::New();
  MovingCastFilterType::Pointer  movingCaster      = MovingCastFilterType::New();
  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();
  ResampleFilterType::Pointer resampleH = ResampleFilterType::New();
  WriterType::Pointer      writerH =  WriterType::New();
  CastFilterType::Pointer  casterH =  CastFilterType::New();
  //Parameter settings
  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
  if(argc==6)
    {
      hitchImageReader->SetFileName( argv[4] );
      writerH->SetFileName( argv[5] );   
    }

  //ITK Pipeline
  fixedCaster->SetInput( fixedImageReader->GetOutput() );
  movingCaster->SetInput( movingImageReader->GetOutput() );
  fixedCaster->Update();
  movingCaster->Update();
  registration->SetFixedImage( fixedCaster->GetOutput() );
  registration->SetMovingImage( movingCaster->GetOutput() );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  registration->SetMetric( metric  );
  registration->SetTransform( transform );
  registration->SetFixedImageRegion( fixedCaster->GetOutput()->GetBufferedRegion() );

  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  fixedImageReader->GetOutput() );
  initializer->SetMovingImage( movingImageReader->GetOutput() );
  initializer->MomentsOn();
  initializer->InitializeTransform();

  registration->SetInitialTransformParameters( transform->GetParameters() );
 
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  optimizerScales[0] = 1.0; // scale for M11
  optimizerScales[1] = 1.0; // scale for M12
  optimizerScales[2] = 1.0; // scale for M13
  optimizerScales[3] = 1.0; // scale for M21
  optimizerScales[4] = 1.0; // scale for M22
  optimizerScales[5] = 1.0; // scale for M23
  optimizerScales[6] = 1.0; // scale for M31
  optimizerScales[7] = 1.0; // scale for M32
  optimizerScales[8] = 1.0; // scale for M33
  optimizerScales[9] = 1.0 / 1e7;  // scale for translation on X
  optimizerScales[10] = 1.0 / 1e7; // scale for translation on Y
  optimizerScales[11] = 1.0 / 1e7; // scale for translation on Z
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

  try
    {
      registration->StartRegistration();
      std::cout << "Optimizer stop condition: "
		<< registration->GetOptimizer()->GetStopConditionDescription()
		<< std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout << "ExceptionObject caught !" << std::endl;
      std::cout << err << std::endl;
      return EXIT_FAILURE;
    }
  
  std::cout << "Optimizer Stopping Condition = "
            << optimizer->GetStopCondition() << std::endl;

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType finalParameters = registration->GetLastTransformParameters();
  unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  double bestValue = optimizer->GetValue();


  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  
  TransformType::Pointer finalTransform = TransformType::New();
  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );
    
  resample->SetTransform( finalTransform );
  resample->SetInput( movingImageReader->GetOutput() );
  resample->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  if(argc==6)
    {
      resampleH->SetTransform( finalTransform );
      resampleH->SetInput( hitchImageReader->GetOutput() );
      resampleH->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
      resampleH->SetOutputOrigin(  fixedImage->GetOrigin() );
      resampleH->SetOutputSpacing( fixedImage->GetSpacing() );
      resampleH->SetOutputDirection( fixedImage->GetDirection() );
      casterH->SetInput( resampleH->GetOutput() );
      writerH->SetInput( casterH->GetOutput()   );
      writerH->Update();
    }

  return EXIT_SUCCESS;
}
