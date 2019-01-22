#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile binaryMaskImageFile outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int InputPixelType;
  typedef unsigned char OutputPixelType;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::LabelStatisticsImageFilter< InputImageType, InputImageType > LabelStatisticsFilterType;
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New();
  LabelStatisticsFilterType::Pointer labelStatisticsFilter = LabelStatisticsFilterType::New();
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  readerM->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );

  //Pipeline
  reader->Update();
  labelStatisticsFilter->SetInput(reader->GetOutput());
  labelStatisticsFilter->SetLabelInput( readerM->GetOutput() );  
  labelStatisticsFilter->Update();
  typedef LabelStatisticsFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  typedef LabelStatisticsFilterType::LabelPixelType                LabelPixelType;
  ValidLabelValuesType::const_iterator vIt=labelStatisticsFilter->GetValidLabelValues().begin()+1;
  LabelPixelType labelValue = *vIt;
  int minVal = labelStatisticsFilter->GetMinimum( labelValue );
  int maxVal = labelStatisticsFilter->GetMaximum( labelValue );
  std::cout<<"min: "<<minVal<<"; max: "<<maxVal<<std::endl;
  
  //Get image specs
  InputImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  InputImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  InputImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  InputImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  InputImageType::Pointer image = InputImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->Allocate();  
  //Perform intensity mapping
  InputImageType::IndexType pixelIndex;
  int i, j, k, value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  if(value < minVal)
	    value = minVal;
	  if(value > maxVal)
	    value = maxVal;
	  //Set image
	  image->SetPixel(pixelIndex,value);
	}
  //Rescale to 0~255
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum( 0 );
  rescaleFilter->SetOutputMaximum( 255 );
  //Write output  
  writer->SetInput( rescaleFilter->GetOutput() );
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
