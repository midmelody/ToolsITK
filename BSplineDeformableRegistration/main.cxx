#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSBOptimizerv4.h"
#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkBSplineTransformInitializer.h"
#include "itkTransformToDisplacementFieldFilter.h"
#include "itkCommand.h"

class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:
  typedef   itk::LBFGSBOptimizerv4   OptimizerType;
  typedef   const OptimizerType *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
  {
  Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
  {
  OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
  if( !(itk::IterationEvent().CheckEvent( &event )) )
    {
    return;
    }
  std::cout << optimizer->GetCurrentIteration() << "   ";
  std::cout << optimizer->GetCurrentMetricValue() << "   ";
  std::cout << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
  }
};


int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Missing Parameters " << std::endl;
      std::cerr << "Usage: " << argv[0];
      std::cerr << " fixedImageFile  movingImageFile outputImagefile ";
      std::cerr << " [hitchImageFile] [hitchOutputImageFile] ";
      return EXIT_FAILURE;
    }

  //ITK Paramters
  const   unsigned int ImageDimension = 3;
  const   unsigned int SpaceDimension = ImageDimension;
  const   unsigned int SplineOrder = 3;
  typedef float PixelType;
  typedef unsigned char OutputPixelType;
  typedef double CoordinateRepType;
  typedef itk::Image< PixelType, ImageDimension > ImageType;
  typedef itk::Image< OutputPixelType, ImageDimension > OutputImageType;
  typedef itk::BSplineTransform< CoordinateRepType, SpaceDimension, SplineOrder > TransformType;
  typedef itk::BSplineTransformInitializer< TransformType, ImageType> InitializerType;
  typedef itk::LBFGSBOptimizerv4 OptimizerType;
  typedef itk::MeanSquaresImageToImageMetricv4< ImageType, ImageType > MetricType;
  typedef itk::ImageRegistrationMethodv4< ImageType, ImageType > RegistrationType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
  typedef itk::CastImageFilter< ImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  //ITK Filter Definitions
  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  TransformType::Pointer      outputBSplineTransform = TransformType::New();
  InitializerType::Pointer    transformInitializer = InitializerType::New();
  ImageReaderType::Pointer    fixedImageReader  = ImageReaderType::New();
  ImageReaderType::Pointer    movingImageReader = ImageReaderType::New();
  ImageReaderType::Pointer    hitchImageReader  = ImageReaderType::New();
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  ResampleFilterType::Pointer resamplerH = ResampleFilterType::New();
  CastFilterType::Pointer     caster =  CastFilterType::New();
  CastFilterType::Pointer     casterH =  CastFilterType::New();
  WriterType::Pointer         writer =  WriterType::New();
  WriterType::Pointer         writerH =  WriterType::New();

  //Parameter settings
  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
  if(argc==6)
    {
      hitchImageReader->SetFileName( argv[4] );
      writerH->SetFileName( argv[5] );   
    }
  unsigned int numberOfGridNodesInOneDimension = 5;

  //ITK Pipeline
  //Read Inputs
  ImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
  fixedImageReader->Update();
  if(argc==6)
    {
      ImageType::ConstPointer hitchImage = hitchImageReader->GetOutput();
      hitchImageReader->Update();
    }

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetFixedImage(  fixedImage   );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );
  
  TransformType::MeshSizeType             meshSize;
  meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );

  transformInitializer->SetTransform( outputBSplineTransform );
  transformInitializer->SetImage( fixedImage );
  transformInitializer->SetTransformDomainMeshSize( meshSize );
  transformInitializer->InitializeTransform();

  typedef TransformType::ParametersType     ParametersType;
  const unsigned int numberOfParameters = outputBSplineTransform->GetNumberOfParameters();
  ParametersType parameters( numberOfParameters );
  parameters.Fill( 0.0 );
  outputBSplineTransform->SetParameters( parameters );
  
  registration->SetInitialTransform( outputBSplineTransform );
  registration->InPlaceOn();
  const unsigned int numberOfLevels = 1;

  RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize( numberOfLevels );
  shrinkFactorsPerLevel[0] = 1;

  RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize( numberOfLevels );
  smoothingSigmasPerLevel[0] = 0;

  registration->SetNumberOfLevels( numberOfLevels );
  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

  const unsigned int numParameters = outputBSplineTransform->GetNumberOfParameters();
  OptimizerType::BoundSelectionType boundSelect( numParameters );
  OptimizerType::BoundValueType upperBound( numParameters );
  OptimizerType::BoundValueType lowerBound( numParameters );
  
  boundSelect.Fill( OptimizerType::UNBOUNDED );
  upperBound.Fill( 0.0 );
  lowerBound.Fill( 0.0 );

  optimizer->SetBoundSelection( boundSelect );
  optimizer->SetUpperBound( upperBound );
  optimizer->SetLowerBound( lowerBound );

  optimizer->SetCostFunctionConvergenceFactor( 1e+12 );
  optimizer->SetGradientConvergenceTolerance( 1.0e-35 );
  optimizer->SetNumberOfIterations( 500 );
  optimizer->SetMaximumNumberOfFunctionEvaluations( 500 );
  optimizer->SetMaximumNumberOfCorrections( 5 );
  
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  
  //Start registration
  std::cout << "Starting Registration " << std::endl;
  try
    {
      registration->Update();
      std::cout << "Optimizer stop condition = " 
		<< registration->GetOptimizer()->GetStopConditionDescription()
                << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }

  OptimizerType::ParametersType finalParameters = outputBSplineTransform->GetParameters();

  //Use last transform to resample the input images
  resampler->SetTransform( outputBSplineTransform );
  resampler->SetInput( movingImageReader->GetOutput() );
  resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin( fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 100 );
  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput() );
  writer->Update();

  if(argc==6)
    {
      resamplerH->SetTransform( outputBSplineTransform );
      resamplerH->SetInput( hitchImageReader->GetOutput() );
      resamplerH->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
      resamplerH->SetOutputOrigin( fixedImage->GetOrigin() );
      resamplerH->SetOutputSpacing( fixedImage->GetSpacing() );
      resamplerH->SetOutputDirection( fixedImage->GetDirection() );
      resamplerH->SetDefaultPixelValue( 100 );
      casterH->SetInput( resamplerH->GetOutput() );
      writerH->SetInput( casterH->GetOutput() );
      writerH->Update();
    }

  return EXIT_SUCCESS;
}


  /*
  // Generate the explicit deformation field resulting from
  // the registration.
  typedef itk::Vector< float, ImageDimension >          VectorPixelType;
  typedef itk::Image< VectorPixelType, ImageDimension > DisplacementFieldImageType;

  typedef itk::TransformToDisplacementFieldFilter<
                        DisplacementFieldImageType,
                        CoordinateRepType >             DisplacementFieldGeneratorType;

  // Create an setup displacement field generator. 
  DisplacementFieldGeneratorType::Pointer dispfieldGenerator =
                                                 DisplacementFieldGeneratorType::New();
  dispfieldGenerator->UseReferenceImageOn();
  dispfieldGenerator->SetReferenceImage( fixedImage );
  dispfieldGenerator->SetTransform( outputBSplineTransform );
  try
    {
    dispfieldGenerator->Update();
    }
  catch ( itk::ExceptionObject & err )
    {
    std::cerr << "Exception detected while generating deformation field";
    std::cerr << " : "  << err << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter< DisplacementFieldImageType >  FieldWriterType;
  FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

  fieldWriter->SetInput( dispfieldGenerator->GetOutput() );

  if( argc >= 7 )
    {
    fieldWriter->SetFileName( argv[6] );
    try
      {
      fieldWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    }
*/
