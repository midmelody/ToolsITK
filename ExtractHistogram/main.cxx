#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageToHistogramFilter.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile histogramFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::StatisticsImageFilter< ImageType > StatisticsImageFilterType;
  typedef itk::Statistics::ImageToHistogramFilter< ImageType > ImageToHistogramFilterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
  ImageToHistogramFilterType::Pointer imageToHistogramFilter = ImageToHistogramFilterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  reader->Update();
  statisticsImageFilter->SetInput( reader->GetOutput() );
  statisticsImageFilter->Update();
  imageToHistogramFilter->SetInput( reader->GetOutput() );
  int max, min;
  min = statisticsImageFilter->GetMinimum();
  max = statisticsImageFilter->GetMaximum();
  const unsigned int MeasurementVectorSize = 1; // Grayscale
  const unsigned int binsPerDimension = max-min+1;
  ImageToHistogramFilterType::HistogramType::MeasurementVectorType lowerBound(binsPerDimension);
  lowerBound.Fill(min);
  ImageToHistogramFilterType::HistogramType::MeasurementVectorType upperBound(binsPerDimension);
  upperBound.Fill(max);
  ImageToHistogramFilterType::HistogramType::SizeType size(MeasurementVectorSize);
  size.Fill(binsPerDimension);
 
  imageToHistogramFilter->SetHistogramBinMinimum( lowerBound );
  imageToHistogramFilter->SetHistogramBinMaximum( upperBound );
  imageToHistogramFilter->SetHistogramSize( size );
  imageToHistogramFilter->Update();
  
  ImageToHistogramFilterType::HistogramType* histogram = imageToHistogramFilter->GetOutput();
  FILE *histFile = fopen(argv[2], "w");
  int freq;
  std::cout<<"Range: "<<min<<" "<<max<<std::endl;
  for(unsigned int i = 0; i < histogram->GetSize()[0]; ++i)
    {
      freq = histogram->GetFrequency(i);
      fprintf(histFile, "%d %d \n", min+i, freq);
    }
  fclose(histFile);

  return EXIT_SUCCESS;
}
