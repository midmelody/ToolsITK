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
    std::cerr << argv[0] << "InputImageName OutputImageName percentage" << std::endl;
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
  short percentage = atoi( argv[3] );

  //Pipeline
  imageToHistogramFilter->SetInput( reader->GetOutput() );
  imageToHistogramFilter->SetHistogramBinMinimum( lowerBound );
  imageToHistogramFilter->SetHistogramBinMaximum( upperBound );
  imageToHistogramFilter->SetHistogramSize( size );
  imageToHistogramFilter->Update();

  ImageToHistogramFilterType::HistogramType* histogram = imageToHistogramFilter->GetOutput();
  unsigned int hist[histogram->GetSize()[0]];
  int i, j;
  for(i=0; i<histogram->GetSize()[0]; i++)
    {
      hist[i] = histogram->GetFrequency(histogram->GetSize()[0]-1-i);
      //std::cout<<histogram->GetSize()[0]-1-i<<" "<<hist[i]<<std::endl;
    }

  float pre, curr, post;
  int preIdx, postIdx;
  int peakCt = 0;
  bool flag = true;
  int peakIdx[2];
  for(i=0; i<histogram->GetSize()[0]; i++)
    {
      curr = hist[i];
      pre = 0;
      post = 0;
      for(j=1;j<=5;j++)
	{
	  preIdx = i-j;
	  postIdx = i+j;
	  if(preIdx>=0)
	    pre = pre+hist[preIdx];
	  if(postIdx<histogram->GetSize()[0])
	    post = post+hist[postIdx];
	}
      
      if((curr>pre)&&(curr>3000)&&(curr>post))
	{
	  peakIdx[peakCt] = histogram->GetSize()[0]-1-i;
	  peakCt++;
	}

      if(peakCt==2)
	{
	  /*
	  if(i>20)
	    thresholdLow = (histogram->GetSize()[0]-1-i+10)*100; 
	  else if(i>10)
 	    thresholdLow = (histogram->GetSize()[0]-1-i+5)*100;   
	  else
	  thresholdLow = (histogram->GetSize()[0]-1-i)*100;	  
	  */
	  thresholdLow = (peakIdx[1]+(peakIdx[0]-peakIdx[1])*percentage/100)*100;
	  std::cout<<"Threshold is "<<thresholdLow<<std::endl;
	  flag = false;
	  break;
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
