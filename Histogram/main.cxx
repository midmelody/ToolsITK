#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageToHistogramFilter.h"
#include "itkImageRandomIteratorWithIndex.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile minValue maxValue numberOfBins" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::Statistics::ImageToHistogramFilter< ImageType > ImageToHistogramFilterType;

  const unsigned int MeasurementVectorSize = 1; //Grayscale
  const unsigned int minValue = static_cast< unsigned int >( atoi( argv[2] ) );
  const unsigned int maxValue = static_cast< unsigned int >( atoi( argv[3] ) );
  const unsigned int binsPerDimension = static_cast< unsigned int >( atoi( argv[4] ) );
    
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  ImageType::Pointer image = reader->GetOutput();

  ImageToHistogramFilterType::HistogramType::MeasurementVectorType lowerBound(binsPerDimension);
  lowerBound.Fill(minValue);
  ImageToHistogramFilterType::HistogramType::MeasurementVectorType upperBound(binsPerDimension);
  upperBound.Fill(maxValue);
  ImageToHistogramFilterType::HistogramType::SizeType size(MeasurementVectorSize);
  size.Fill(binsPerDimension);
 
  ImageToHistogramFilterType::Pointer imageToHistogramFilter = ImageToHistogramFilterType::New();
  imageToHistogramFilter->SetInput( image );
  imageToHistogramFilter->SetHistogramBinMinimum( lowerBound );
  imageToHistogramFilter->SetHistogramBinMaximum( upperBound );
  imageToHistogramFilter->SetHistogramSize( size );

  imageToHistogramFilter->Update();

  ImageToHistogramFilterType::HistogramType* histogram = imageToHistogramFilter->GetOutput();

  for(unsigned int i = 0; i < histogram->GetSize()[0]; ++i)
    std::cout << histogram->GetFrequency(i) << std::endl;

  return EXIT_SUCCESS;
}
