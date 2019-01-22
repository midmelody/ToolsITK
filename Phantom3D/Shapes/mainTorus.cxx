#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#define PI 3.1415926
#include "DrawShape3D.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters
  float x,y,z; 
  if( argc < 1 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] <<  std::endl;
      return EXIT_FAILURE;
    }
  ImageType::SpacingType spacing;
  spacing[0] = 1;
  spacing[1] = 1;
  spacing[2] = 1;
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];

  ImageType::SizeType  size;
  size[0] = 250;
  size[1] = 250;
  size[2] = 120;
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

  float torusR = 70;
  float tubeR = 50; 
  for(k=0;k<pSli;k++)
    {
      std::cout<<"\r";
      std::cout<<(k+1)*100/pSli<<"%"<<std::flush;
      for(i=0;i<pRow;i++) 
	for(j=0;j<pCol;j++)
	  {
	    pixelIndex[0]=i;
	    pixelIndex[1]=j;
	    pixelIndex[2]=k;
	    x = float(i)-size[0]/2;
	    y = float(j)-size[1]/2;
	    z = float(k)-size[2]/2;
	    if((torusR-sqrt(x*x+y*y))*(torusR-sqrt(x*x+y*y))+z*z<tubeR*tubeR)
	      outImage->SetPixel(pixelIndex,3);
	  }
    }

  std::cout<<std::endl;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( outImage );
  writer->SetFileName("torus.img");
  writer->Update();
  return EXIT_SUCCESS;
}

