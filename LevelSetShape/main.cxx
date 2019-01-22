// ARGUMENTS: 16 21 21 5 1.0 -0.5 3.0 0.05 1

#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkShapeDetectionLevelSetImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"


int main( int argc, char *argv[] )
{
  if( argc < 13 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage";
    std::cerr << " seedX seedY seedZ InitialDistance";
    std::cerr << " Sigma SigmoidAlpha SigmoidBeta ";
    std::cerr << " curvatureScaling propagationScaling fastMarchingImage" << std::endl;
    return 1;
    }

  typedef float           InternalPixelType;
  typedef unsigned char   OutputPixelType;
  const   unsigned int    Dimension = 3;

  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::ImageFileReader< InternalImageType > ReaderType;
  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< InternalImageType, InternalImageType >  GradientFilterType;
  typedef itk::SigmoidImageFilter< InternalImageType, InternalImageType >  SigmoidFilterType;
  typedef itk::FastMarchingImageFilter< InternalImageType, InternalImageType > FastMarchingFilterType;
  typedef FastMarchingFilterType::NodeContainer           NodeContainer;
  typedef FastMarchingFilterType::NodeType                NodeType;
  typedef itk::ShapeDetectionLevelSetImageFilter< InternalImageType, InternalImageType >    ShapeDetectionFilterType;
  typedef itk::RescaleIntensityImageFilter<InternalImageType, OutputImageType> CastFilterType; 
  typedef itk::BinaryThresholdImageFilter< InternalImageType, OutputImageType > ThresholdingFilterType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
  const double sigma = atof( argv[7] );
  gradientMagnitude->SetSigma(  sigma  );
  gradientMagnitude->SetInput( reader->GetOutput() );

  SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
  const double alpha =  atof( argv[8] );
  const double beta  =  atof( argv[9] );
  sigmoid->SetAlpha( alpha );
  sigmoid->SetBeta(  beta  );
  sigmoid->SetOutputMinimum(  0.0  );
  sigmoid->SetOutputMaximum(  1.0  );
  sigmoid->SetInput( gradientMagnitude->GetOutput() );

  FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
  NodeContainer::Pointer seeds = NodeContainer::New();
  InternalImageType::IndexType  seedPosition;
  seedPosition[0] = atoi( argv[3] );
  seedPosition[1] = atoi( argv[4] );
  seedPosition[2] = atoi( argv[5] );
  const double initialDistance = atof( argv[6] );
  NodeType node;
  const double seedValue = - initialDistance;
  node.SetValue( seedValue );
  node.SetIndex( seedPosition );
  seeds->Initialize();
  seeds->InsertElement( 0, node );
  fastMarching->SetInput( sigmoid->GetOutput() );
  fastMarching->SetTrialPoints(  seeds  );
  fastMarching->SetSpeedConstant( 1.0 );
  fastMarching->SetOutputSize( reader->GetOutput()->GetBufferedRegion().GetSize() );

  ShapeDetectionFilterType::Pointer shapeDetection = ShapeDetectionFilterType::New();
  shapeDetection->SetInput( fastMarching->GetOutput() );
  shapeDetection->SetFeatureImage( sigmoid->GetOutput() );
  const double curvatureScaling   = atof( argv[ 10 ] );
  const double propagationScaling = atof( argv[ 11 ] );
  shapeDetection->SetPropagationScaling(  propagationScaling );
  shapeDetection->SetCurvatureScaling( curvatureScaling );
  shapeDetection->SetMaximumRMSError( 0.02 );
  shapeDetection->SetNumberOfIterations( 800 );

  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  thresholder->SetLowerThreshold( -1000.0 );
  thresholder->SetUpperThreshold(     0.0 );
  thresholder->SetOutsideValue(  0  );
  thresholder->SetInsideValue(  255 );
  thresholder->SetInput( shapeDetection->GetOutput() );

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( thresholder->GetOutput() );
 
  CastFilterType::Pointer caster1 = CastFilterType::New();
  WriterType::Pointer writer1 = WriterType::New();
  caster1->SetInput( sigmoid->GetOutput() );
  writer1->SetInput( caster1->GetOutput() );
  writer1->SetFileName(argv[12]);
  caster1->SetOutputMinimum(   0 );
  caster1->SetOutputMaximum( 255 );

  // Print out some useful information
  std::cout << std::endl;
  std::cout << "Max. no. iterations: " << shapeDetection->GetNumberOfIterations() << std::endl;
  std::cout << "Max. RMS error: " << shapeDetection->GetMaximumRMSError() << std::endl;
  std::cout << std::endl;
  std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;
  std::cout << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;

  writer->Update();
  writer1->Update();

  return 0;
}
