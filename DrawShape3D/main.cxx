#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#define PI 3.1415926
#include "DrawShape3D.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters
  float x, y;
  if( argc < 1 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] <<  std::endl;
      return EXIT_FAILURE;
    }
  ImageType::SpacingType spacing;
  spacing[0] = 0.5;
  spacing[1] = 0.5;
  spacing[2] = 0.5;
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];

  ImageType::SizeType  size;
  size[0] = 150;
  size[1] = 150;
  size[2] = 150;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  
  ImageType::RegionType region;
  ImageType::IndexType start = {0};
  region.SetSize( size );
  region.SetIndex( start );
  ImageType::IndexType pixelIndex;
  
  ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions( region );
  outImage->SetSpacing( spacing );
  outImage->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  outImage->SetPixel(pixelIndex,1);	  
	}

  for(k=10;k<140;k++)
    {
      std::cout<<"\r";
      std::cout<<(k-10+1)*100/130<<"%"<<std::flush;
      for(i=0;i<pRow;i++) 
	for(j=0;j<pCol;j++)
	  {
	    pixelIndex[0]=i;
	    pixelIndex[1]=j;
	    pixelIndex[2]=k;
	    x = float(i)-75;
	    y = float(j)-75;
	    if((x*x+y*y)<2500)
	      outImage->SetPixel(pixelIndex,3);
	  }
    }

  std::cout<<std::endl;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( outImage );
  writer->SetFileName("cylinder.img");
  writer->Update();
  return EXIT_SUCCESS;
}

