#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageToHistogramFilter.h"
#include "itkImageRandomIteratorWithIndex.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileWriter.h"

int main(int argc, char* argv[])
{
  if( argc != 4 )
    {
    std::cerr << argv[0] << "InputImageName OutputImageName Mode(1: first peak at top; 0: first peak under top)" << std::endl;
    return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef short PixelType;
  typedef unsigned char OutPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< OutPixelType, Dimension > OutImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::Statistics::ImageToHistogramFilter< ImageType > ImageToHistogramFilterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, OutImageType > BinaryThresholdFilterType; 
  typedef itk::ImageFileWriter< OutImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ImageToHistogramFilterType::Pointer imageToHistogramFilter = ImageToHistogramFilterType::New();
  BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  const unsigned int MeasurementVectorSize = 1; 
  const unsigned int binsPerDimension = 100;
  ImageToHistogramFilterType::HistogramType::MeasurementVectorType lowerBound(binsPerDimension);
  lowerBound.Fill(0);
  ImageToHistogramFilterType::HistogramType::MeasurementVectorType upperBound(binsPerDimension);
  upperBound.Fill(10000);
  ImageToHistogramFilterType::HistogramType::SizeType size(MeasurementVectorSize);
  size.Fill(binsPerDimension);
  short thresholdLow;
  short thresholdHigh = 10000;
  bool mode = atoi( argv[3] );

  //Pipeline
  imageToHistogramFilter->SetInput( reader->GetOutput() );
  imageToHistogramFilter->SetHistogramBinMinimum( lowerBound );
  imageToHistogramFilter->SetHistogramBinMaximum( upperBound );
  imageToHistogramFilter->SetHistogramSize( size );
  imageToHistogramFilter->Update();

  ImageToHistogramFilterType::HistogramType* histogram = imageToHistogramFilter->GetOutput();
  float accu, start;
  if(mode)
    {
      accu = histogram->GetSize()[0]-5;
      start = histogram->GetSize()[0]-6;
    }
  else
    {
      accu = histogram->GetSize()[0]-1;
      start = histogram->GetSize()[0]-2;
    }
  float curr;
  int count=0;
  bool flag = true;
  for(unsigned int i = start; i > 0; i--)
    {
      curr = histogram->GetFrequency(i);
      std::cout<<curr<<" "<<accu<<std::endl;
      if((curr>accu)&&(curr>10000))
	{
	  if(mode)
	    thresholdLow = (i+5)*100;
	  else
	    thresholdLow = (i-5)*100;	    
	  std::cout<<"Threshold is "<<thresholdLow<<std::endl;
	  flag = false;
	  break;
	}
      accu = accu + histogram->GetFrequency(i);
      count++;
      if(count==5)
	{
	  count = 0;
	  accu = histogram->GetFrequency(i);
	}
    }

  if(flag)
    {
      thresholdLow = 9000;
      std::cout<<"No threshold found!!! Set arbitarily at 9000"<<std::endl;
    }
 
  thresholder->SetInput(reader->GetOutput());
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetLowerThreshold( thresholdLow );
  thresholder->SetUpperThreshold( thresholdHigh );
  writer->SetInput( thresholder->GetOutput() );
  try
    {
      writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
