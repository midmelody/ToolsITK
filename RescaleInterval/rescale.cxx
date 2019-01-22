#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program threshold the original CT image with low and high threshold and map the voxel with HU within region to 0~255

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile thresholdLow thresholdHigh" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef double InputPixelType;
  typedef double OutputPixelType;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  int threLow = atoi(argv[3]);
  int threHigh = atoi(argv[4]);

  //Pipeline
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
	  if(value < threLow)
	    value = threLow;
	  if(value > threHigh)
	    value = threHigh;
	  //Set image
	  image->SetPixel(pixelIndex,value);
	}
  //Rescale to 0~65535
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(   0 );
  rescaleFilter->SetOutputMaximum( 65535 );
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
