#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkStatisticsImageFilter.h"
#include <iostream>
#include "math.h"
#include "string.h"
#include <iomanip>
#include <limits>

int main( int argc, char * argv[] )
{
  if( argc < 2 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
    
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  StatisticsImageFilterType::Pointer statistician = StatisticsImageFilterType::New ();
  reader->SetFileName( argv[1] );
  statistician->SetInput(reader->GetOutput());
  statistician->Update();
 
  std::cout << "Sum: " << std::setprecision (std::numeric_limits<double>::digits10 + 1)  << statistician->GetSum() << std::endl; 
  //std::cout << "Mean: " << statistician->GetMean() << std::endl;
  //std::cout << "Std.: " << statistician->GetSigma() << std::endl;
  //std::cout << "Min: " << statistician->GetMinimum() << std::endl;
  //std::cout << "Max: " << statistician->GetMaximum() << std::endl;


  return EXIT_SUCCESS;
}
