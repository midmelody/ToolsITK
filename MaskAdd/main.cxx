#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile1 inputImageFile2 outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::StatisticsImageFilter< ImageType > StatFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  StatFilterType::Pointer stat = StatFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
 
  //Pipeline
  try
    {
      reader1->Update();
      reader2->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
 
  stat->SetInput(reader1->GetOutput());
  stat->Update();
  int maxlabel1 = stat->GetMaximum();

  //Get image specs 
  ImageType::SpacingType spacing = reader1->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader1->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader1->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader1->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  

  //Assign value
  ImageType::IndexType pixelIndex;
  unsigned int value, value1, value2;
  int i, j, k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value1 = reader1->GetOutput()->GetPixel(pixelIndex);
	  value2 = reader2->GetOutput()->GetPixel(pixelIndex);
	  if((value1==0)&&(value2==0))
	    value = 0;
	  else if(value1 > value2)
	    {
	      value = value1;
	      //std::cout<<value1<<" "<<value2<<" "<<pixelValue<<std::endl;
	    }
	  else
	    {
	      value = value2 + maxlabel1;
	      //std::cout<<value1<<" "<<value2<<" "<<pixelValue<<std::endl;
	    }
	  
	  image->SetPixel(pixelIndex, value);
	}


  StatFilterType::Pointer statO = StatFilterType::New();
  statO->SetInput(image);
  statO->Update();
  int test = statO->GetMaximum();
  std::cout<<test<<std::endl;

  //Write output  
  writer->SetInput( image );
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
