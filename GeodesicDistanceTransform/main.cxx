#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << argv[0];
      std::cerr << " inputImage outputImage seedFile";
      std::cerr << std::endl;
      return 1;
    }

  typedef float           PixelType;
  const   unsigned int    Dimension = 3;
  typedef itk::Image< PixelType, Dimension > ImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
 
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }

  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  
  image->FillBuffer(255); 
 

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType, ImageType > GradientFilterType;
  typedef itk::SigmoidImageFilter< ImageType, ImageType > SigmoidFilterType;
  GradientFilterType::Pointer gradientMagnitude = GradientFilterType::New();
  SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
  sigmoid->SetOutputMinimum(  0.0  );
  sigmoid->SetOutputMaximum(  1.0  );

  typedef itk::FastMarchingImageFilter< ImageType, ImageType > FastMarchingFilterType;
  FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  RescaleFilterType::Pointer rescalePre = RescaleFilterType::New();
  rescalePre->SetOutputMinimum(   0 );
  rescalePre->SetOutputMaximum( 255 );

  rescale->SetOutputMinimum(   0 );
  rescale->SetOutputMaximum( 255 );

  rescalePre->SetInput( reader->GetOutput() );
  gradientMagnitude->SetInput( rescalePre->GetOutput() );
  sigmoid->SetInput( gradientMagnitude->GetOutput() );
  fastMarching->SetInput( sigmoid->GetOutput() );
  rescale->SetInput( fastMarching->GetOutput() );
 

  gradientMagnitude->SetSigma( 1.0 );
  sigmoid->SetAlpha( -0.5 );
  sigmoid->SetBeta( 3.0 );

  //Seed
  typedef FastMarchingFilterType::NodeContainer           NodeContainer;
  typedef FastMarchingFilterType::NodeType                NodeType;
  NodeContainer::Pointer seeds = NodeContainer::New();
  ImageType::IndexType  seedPosition;
  NodeType node;
  const double seedValue = 0.0;
  node.SetValue( seedValue );
  seeds->Initialize();
  //Read seed
  int seedCount;
  int x, y, z;
  std::cout<<"Reading seed file......."<<std::flush;
  FILE *seedFile = fopen(argv[3], "r");
  if(seedFile == NULL) 
    {
      std::cerr << argv[1] << " Seed file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  seedCount = 0;
  while(fscanf(seedFile, "%d %d %d", &x, &y, &z) == 3)  
    {
      seedPosition[0] = x;
      seedPosition[1] = y;
      seedPosition[2] = z;
      node.SetIndex( seedPosition );
      seeds->InsertElement( seedCount, node );
      seedCount++;
    }
  fclose(seedFile);
  std::cout<<"finished reading "<<seedCount<<" seeds"<<std::endl;
  std::cout<<std::endl;
  if(seedCount)
    writer->SetInput( rescale->GetOutput() );
  else
    writer->SetInput( image );

  fastMarching->SetTrialPoints(  seeds  );
  fastMarching->SetOutputSize( reader->GetOutput()->GetBufferedRegion().GetSize() );
  fastMarching->SetStoppingValue( 100 );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  return 0;
}
