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
  spacing[0] = 1;
  spacing[1] = 1;
  spacing[2] = 1;
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];

  ImageType::SizeType  size;
  size[0] = 250;
  size[1] = 250;
  size[2] = 250;
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

  for(k=10;k<240;k++)
    for(i=0;i<pRow;i++) 
      for(j=0;j<pCol;j++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  x = float(i)-50;
	  y = float(j)-50;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(i)-50;
	  y = float(j)-200;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(i)-200;
	  y = float(j)-50;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(i)-200;
	  y = float(j)-200;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	}

  for(i=10;i<240;i++)
    for(k=0;k<pSli;k++) 
      for(j=0;j<pCol;j++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  x = float(k)-50;
	  y = float(j)-50;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(k)-50;
	  y = float(j)-200;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(k)-200;
	  y = float(j)-50;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(k)-200;
	  y = float(j)-200;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	}

  for(j=10;j<240;j++)
    for(i=0;i<pRow;i++) 
      for(k=0;k<pCol;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  x = float(i)-50;
	  y = float(k)-50;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(i)-50;
	  y = float(k)-200;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(i)-200;
	  y = float(k)-50;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	  x = float(i)-200;
	  y = float(k)-200;
	  if((x*x+y*y)<25*25)
	    outImage->SetPixel(pixelIndex,3);
	}





  std::cout<<std::endl;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( outImage );
  writer->SetFileName("tubeBox.img");
  writer->Update();
  return EXIT_SUCCESS;
}

