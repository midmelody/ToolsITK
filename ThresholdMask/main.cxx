#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 9 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile maskImageFile outputImageFile maskLabel1 threshold1 maskLabel2 threshold2" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef unsigned char OutPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< OutPixelType, Dimension > OutImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< OutImageType > WriterType;
    
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
   
  //Parameters
  reader->SetFileName( argv[1] );
  readerM->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
  int maskLabel1 = atoi( argv[4] );
  int threshold1 = atoi( argv[5] );
  int maskLabel2 = atoi( argv[6] );
  int threshold2 = atoi( argv[7] );
  int threshold3 = atoi( argv[8] );
  //Pipeline
  try
    {
      reader->Update();
      readerM->Update();     
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
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  OutImageType::Pointer image = OutImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  

  //Read image
  ImageType::IndexType pixelIndex;
  int i, j, k, value, label;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  label = readerM->GetOutput()->GetPixel(pixelIndex);
	  if(label == maskLabel1)
	    {
	      if(value>=threshold1)
		image->SetPixel(pixelIndex, 1);
	      else
		image->SetPixel(pixelIndex, 0);
	    }
	  else if(label == maskLabel2)
	    {
	      if(value>=threshold2)
		image->SetPixel(pixelIndex, 1);
	      else
		image->SetPixel(pixelIndex, 0);
	    }
	  else
	    {
	      if(value>=threshold3)
		image->SetPixel(pixelIndex, 1);
	      else
		image->SetPixel(pixelIndex, 0);		  
	    }
	}
  
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
