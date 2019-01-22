#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkLabelStatisticsImageFilter.h"
 
#include <iostream>
#include "math.h"

int main(int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile inputLabelImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New();
  
  //Parameters
  reader->SetFileName( argv[1] );
  readerM->SetFileName( argv[2] );
  
  //Pipeline 
  typedef itk::LabelStatisticsImageFilter< ImageType, ImageType > LabelStatisticsImageFilterType;
  LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetInput(reader->GetOutput());
  labelStatisticsImageFilter->SetLabelInput( readerM->GetOutput() );  
  labelStatisticsImageFilter->Update();
    
  //std::cout << "Number of labels: " << labelStatisticsImageFilter->GetNumberOfLabels()-1 << std::endl;
  //std::cout << std::endl;
  
  typedef LabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  typedef LabelStatisticsImageFilterType::LabelPixelType                LabelPixelType;
 
  for(ValidLabelValuesType::const_iterator vIt=labelStatisticsImageFilter->GetValidLabelValues().begin()+1;
      vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
      ++vIt)
    {
    if ( labelStatisticsImageFilter->HasLabel(*vIt) )
      {
      LabelPixelType labelValue = *vIt;
      std::cout << "min: " << labelStatisticsImageFilter->GetMinimum( labelValue ) << "; ";
      std::cout << "max: " << labelStatisticsImageFilter->GetMaximum( labelValue ) << "; ";
      //std::cout << "median: " << labelStatisticsImageFilter->GetMedian( labelValue ) << std::endl;
      //std::cout << "mean: " << labelStatisticsImageFilter->GetMean( labelValue ) << std::endl;
      //std::cout << "sigma: " << labelStatisticsImageFilter->GetSigma( labelValue ) << std::endl;
      //std::cout << "variance: " << labelStatisticsImageFilter->GetVariance( labelValue ) << std::endl;
      //std::cout << "sum: " << labelStatisticsImageFilter->GetSum( labelValue ) << std::endl;
      std::cout << "count: " << labelStatisticsImageFilter->GetCount( labelValue ) << std::endl;
      //std::cout << "box: " << labelStatisticsImageFilter->GetBoundingBox( labelValue ) << std::endl; // can't output a box
      //std::cout << "region: " << labelStatisticsImageFilter->GetRegion( labelValue ) << std::endl;
      //std::cout << std::endl << std::endl;
 
      }
    }
 
  return EXIT_SUCCESS;
}
