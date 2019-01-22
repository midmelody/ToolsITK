#include "itkImage.h"
#include "itkImageFileReader.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile inputMaskFile" << std::endl;
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
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  size[0] = size[0];
  size[1] = size[1];
  size[2] = size[2];
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
  float voxelSize = spacing[0]*spacing[1]*spacing[2];

  //Count;
  int ct = 0;
  double sum = 0;
  ImageType::IndexType index;
  ImageType::PixelType pixel;
  ImageType::PixelType pixelM;
  int i,j,k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  index[0] = i;
	  index[1] = j;
	  index[2] = k;
	  pixel = reader->GetOutput()->GetPixel( index );
	  pixelM = readerM->GetOutput()->GetPixel( index );
	  if(pixelM!=0)
	    {
	      ct++;
	      sum = sum + pixel; 
	    }
	}
  sum = sum/float(ct);
  
  std::cout<<"MeanHU "<<sum<<" "<<"VoxelCt "<<ct<<" voxels; volume "<<int(float(ct)*voxelSize/1000)<<"; "<<spacing[0]<<" "<<spacing[1]<<" "<<spacing[2]<<" "<< voxelSize <<std::endl;

  return EXIT_SUCCESS;
}
