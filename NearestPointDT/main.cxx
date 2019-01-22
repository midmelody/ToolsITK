#include <stdlib.h>
#include <stdio.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " distanceImage binaryImage ";
      std::cerr << " " << std::endl;
      return EXIT_FAILURE;
    }

  typedef float PixelType;
  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer readerDT = ReaderType::New();
  ReaderType::Pointer readerSeg = ReaderType::New();
  readerDT->SetFileName( argv[1] );
  readerSeg->SetFileName( argv[2] );

  try
    {
      readerDT->Update();
      readerSeg->Update();
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem reading image file : " << argv[1] << std::endl;
      std::cerr << excp << std::endl;
      return -1;
    }

  int i,j,k; //counters

  //Image specs
  int pRow,pCol,pSli; //Size of Image: row, column, slice
  int spaRow, spaCol, spaSli; //Spacing: row, column, slice
  ImageType::SpacingType spacing = readerDT->GetOutput()->GetSpacing();
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];
  ImageType::SizeType  size = readerDT->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::IndexType pixelIndex;

  //Find nearest in airway
  float min = 1000000;
  int label;
  float dist;
  int curI, curJ, curK;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  dist = readerDT->GetOutput()->GetPixel(pixelIndex);
	  label = readerSeg->GetOutput()->GetPixel(pixelIndex);
	  if(label==1)
	    if(dist<min)
	      {
		min = dist;
		curI = i;
		curJ = j;
		curK = k;
	      }
	}
  std::cout<<curI<<" "<<curJ<<" "<<curK<<": "<<min<<std::endl;

  return EXIT_SUCCESS;
}


