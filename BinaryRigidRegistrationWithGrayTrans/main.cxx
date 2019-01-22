#include "itkImageRegistrationMethod.h"
#include "itkKappaStatisticImageToImageMetric.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <iostream>

// Command observer monitoring the evolution of the registration process.
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate   Self;
  typedef itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
  typedef const OptimizerType *              OptimizerPointer;
  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }
  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
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
  if( argc < 6 )
    {
      std::cerr << "Usage: " << argv[0];
      std::cerr << " fixedBinaryImageFile movingBinaryImageFile movingGrayScaleImageFile outputBinaryImageFile outputGrayScaleImageFile" << std::endl;
      return EXIT_FAILURE;
    }
  
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef short GrayPixelType;
  typedef itk::Image< GrayPixelType, Dimension >  GrayImageType;

  //Transform, optimizer, metric, interpolator: registration
  typedef itk::VersorRigid3DTransform< double > TransformType;
  typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
  typedef itk::KappaStatisticImageToImageMetric< ImageType, ImageType > MetricType;
  typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double > InterpolatorType;
  typedef itk::ImageRegistrationMethod< ImageType, ImageType > RegistrationType;
  MetricType::Pointer metric = MetricType::New();
  OptimizerType::Pointer optimizer = OptimizerType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  TransformType::Pointer  transform = TransformType::New();
  RegistrationType::Pointer registration = RegistrationType::New();
  metric->ComplementOn();
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetInterpolator( interpolator );
  registration->SetTransform( transform );

  //Input reader and rescaler
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer  fixedImageReader = ImageReaderType::New();
  ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
  fixedImageReader->SetFileName( argv[1] );
  movingImageReader->SetFileName( argv[2] );
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > InRescaleFilterType;
  InRescaleFilterType::Pointer fixedRescaleFilter = InRescaleFilterType::New();
  InRescaleFilterType::Pointer movingRescaleFilter = InRescaleFilterType::New();
  fixedRescaleFilter->SetInput( fixedImageReader->GetOutput() );
  movingRescaleFilter->SetInput( movingImageReader->GetOutput() );
  fixedRescaleFilter->SetOutputMinimum( 0 );
  fixedRescaleFilter->SetOutputMaximum( 255 );
  movingRescaleFilter->SetOutputMinimum( 0 );
  movingRescaleFilter->SetOutputMaximum( 255 );
  fixedRescaleFilter->Update();
  movingRescaleFilter->Update();

  //Initialize and start registration
  registration->SetFixedImage( fixedRescaleFilter->GetOutput() );
  registration->SetMovingImage( movingRescaleFilter->GetOutput() );
  registration->SetFixedImageRegion( fixedRescaleFilter->GetOutput()->GetBufferedRegion() );
  typedef itk::CenteredTransformInitializer< TransformType, ImageType, ImageType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform( transform );
  initializer->SetFixedImage( fixedRescaleFilter->GetOutput() );
  initializer->SetMovingImage( movingRescaleFilter->GetOutput() );
  //Use center of mass
  initializer->MomentsOn();
  initializer->InitializeTransform();
  //Set and get transform info
  typedef TransformType::VersorType VersorType;
  typedef VersorType::VectorType VectorType;
  VersorType rotation;
  VectorType axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  const double angle = 0;
  rotation.Set( axis, angle );
  transform->SetRotation( rotation );
  registration->SetInitialTransformParameters( transform->GetParameters() );
  //Optimizer parameters
  typedef OptimizerType::ScalesType OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  optimizer->SetScales( optimizerScales );
  optimizer->SetMaximumStepLength( 0.500  );
  optimizer->SetMinimumStepLength( 0.005 );
  optimizer->SetNumberOfIterations( 200 );
  //optimizer->SetParametersConvergenceTolerance( 0.1 );  // 1/10th pixel
  //optimizer->SetFunctionConvergenceTolerance(0.0001);    // 0.001 bits
  optimizer->SetGradientMagnitudeTolerance( 1e-5 );

  // Create the Command observer and register it with the optimizer.
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  try
    {
      registration->Update();
      std::cout << "Optimizer stop condition: "
		<< registration->GetOptimizer()->GetStopConditionDescription()
		<< std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
  const double versorX = finalParameters[0];
  const double versorY = finalParameters[1];
  const double versorZ = finalParameters[2];
  const double finalTranslationX = finalParameters[3];
  const double finalTranslationY = finalParameters[4];
  const double finalTranslationZ = finalParameters[5];
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const double bestValue = optimizer->GetValue();

  // Print out results
  std::cout << std::endl << std::endl;
  std::cout << " Result = " << std::endl;
  std::cout << " versor X      = " << versorX  << std::endl;
  std::cout << " versor Y      = " << versorY  << std::endl;
  std::cout << " versor Z      = " << versorZ  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

  //Transform moving image and print result
  transform->SetParameters( finalParameters );
  TransformType::MatrixType matrix = transform->GetMatrix();
  TransformType::OffsetType offset = transform->GetOffset();
  std::cout << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << std::endl << offset << std::endl;

  //Resample binary image according to final transform
  typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
  TransformType::Pointer finalTransform = TransformType::New();
  finalTransform->SetCenter( transform->GetCenter() );
  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetTransform( finalTransform );
  resampler->SetInput( movingRescaleFilter->GetOutput() );
  ImageType::Pointer fixedImage = fixedRescaleFilter->GetOutput();
  resampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 0 );
  resampler->SetInterpolator( interpolator );

  //Resample grayscale image according to final transform
  typedef itk::ImageFileReader< GrayImageType > GrayImageReaderType;
  GrayImageReaderType::Pointer grayImageReader = GrayImageReaderType::New();
  grayImageReader->SetFileName( argv[3] );
  grayImageReader->Update();
  typedef itk::ResampleImageFilter< GrayImageType, GrayImageType > GrayResampleFilterType;
  GrayResampleFilterType::Pointer grayResampler = GrayResampleFilterType::New();
  typedef itk::LinearInterpolateImageFunction< GrayImageType, double > GrayInterpolatorType;
  GrayInterpolatorType::Pointer grayInterpolator = GrayInterpolatorType::New();
  grayResampler->SetTransform( finalTransform );
  grayResampler->SetInput( grayImageReader->GetOutput() );
  grayResampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  grayResampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  grayResampler->SetOutputSpacing( fixedImage->GetSpacing() );
  grayResampler->SetOutputDirection( fixedImage->GetDirection() );
  grayResampler->SetDefaultPixelValue( 0 );
  grayResampler->SetInterpolator( grayInterpolator );

  //Write Mask
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > OutRescaleFilterType;
  WriterType::Pointer writer =  WriterType::New();
  OutRescaleFilterType::Pointer  outRescaleFilter =  OutRescaleFilterType::New();
  writer->SetFileName( argv[4] );
  outRescaleFilter->SetOutputMinimum( 0 );
  outRescaleFilter->SetOutputMaximum( 1 );
  outRescaleFilter->SetInput( resampler->GetOutput() );
  writer->SetInput( outRescaleFilter->GetOutput()   );
  writer->Update();

  //Write Grayscale Image
  typedef itk::ImageFileWriter< GrayImageType > GrayWriterType;
  GrayWriterType::Pointer grayWriter = GrayWriterType::New();
  grayWriter->SetFileName( argv[5] );
  grayWriter->SetInput( grayResampler->GetOutput() );
  grayWriter->Update();

  return EXIT_SUCCESS;
}
