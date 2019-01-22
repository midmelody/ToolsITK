#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#define PI 3.1415926
#include "DrawShape3D.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters
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


  int factor = 0.6;


  ImageType::SizeType  size;
  size[0] = 500;
  size[1] = 500;
  size[2] = 500;
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
  
  float x,y,z;
  int a,b,c;
  float sizeF = 500;
  float mag = 0.05*sizeF;
  i=0;
  j=0;
  k=0;
  while(k<pSli)
    {
      std::cout<<"\r";
      std::cout<<(k+1)*100/pSli<<"%"<<std::flush;
     
      for(i=0;i<pRow;i++)
	{
	  for(j=2;j<7;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*12*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	  for(j=15;j<20;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*12*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	  for(j=28;j<33;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*12*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	  for(j=89;j<99;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*10*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	  for(j=109;j<119;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*10*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	  for(j=129;j<139;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*10*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	  for(j=209;j<229;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*8*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	  for(j=249;j<269;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*8*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	  for(j=289;j<309;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*8*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,1);
	    }
	}
      int centerA0,centerA1,centerB0, centerB10, centerB11;    
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {	
	    x = i;
	    y = j;
	    centerA0 = round(mag*0.6*sin(float(k)/sizeF*8*PI)+70);
	    centerB0 = 425;
	    if(((x-centerA0)*(x-centerA0)+(y-centerB0)*(y-centerB0))<900)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    centerA1 = round(mag*0.6*sin(float(k)/sizeF*8*PI)+130);
	    centerB10 = 400;
	    centerB11 = 450;
	    if(((x-centerA1)*(x-centerA1)+(y-centerB10)*(y-centerB10))<400)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    if(((x-centerA1)*(x-centerA1)+(y-centerB11)*(y-centerB11))<400)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    centerA1 = round(mag*0.6*sin(float(k)/sizeF*8*PI)+180);
	    centerB10 = 380;
	    centerB0 = 425;
	    centerB11 = 470;
	    if(((x-centerA1)*(x-centerA1)+(y-centerB0)*(y-centerB0))<225)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    if(((x-centerA1)*(x-centerA1)+(y-centerB10)*(y-centerB10))<225)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    if(((x-centerA1)*(x-centerA1)+(y-centerB11)*(y-centerB11))<225)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	  }

      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {	
	    x = i;
	    y = j;
	    centerA0 = 350;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+380);
	    if((float((x-centerA0)*(x-centerA0))/1600+float((y-centerB0)*(y-centerB0))/225)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    centerA0 = 305;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+420);
	    if((float((x-centerA0)*(x-centerA0))/1600+float((y-centerB0)*(y-centerB0))/225)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    centerA0 = 395;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+420);
	    if((float((x-centerA0)*(x-centerA0))/1600+float((y-centerB0)*(y-centerB0))/225)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    centerA0 = 260;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+460);
	    if((float((x-centerA0)*(x-centerA0))/1600+float((y-centerB0)*(y-centerB0))/225)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    centerA0 = 350;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+460);
	    if((float((x-centerA0)*(x-centerA0))/1600+float((y-centerB0)*(y-centerB0))/225)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }
	    centerA0 = 440;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+460);
	    if((float((x-centerA0)*(x-centerA0))/1600+float((y-centerB0)*(y-centerB0))/225)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,1);
	      }


	  }
















      
      k=k+1;
    }
  std::cout<<std::endl;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( outImage );
  writer->SetFileName("phantom.img");
  writer->Update();
  return EXIT_SUCCESS;
}



/*
      for(i=0;i<pRow;k++) 
	for(j=0;j<pCol;j++)
	  {
	    pixelIndex[0]=i;
	    pixelIndex[1]=j;
	    pixelIndex[2]=k;
	  x = float(i)-49;
	  y = float(j)-49;
	  /*
	  z = float(k)-49;
	  if((x*x+y*y+z*z)<900)
	    outImage->SetPixel(pixelIndex,3);
	  else if((x*x+y*y+z*z)<1600)
	    outImage->SetPixel(pixelIndex,2);
	  else
	    outImage->SetPixel(pixelIndex,1);
	  */

/*
	  if((x*x+y*y)<900)
	    outImage->SetPixel(pixelIndex,3);
	  else if((x*x+y*y)<1600)
	    outImage->SetPixel(pixelIndex,2);
	  else
	    outImage->SetPixel(pixelIndex,1);

	  z = float(k)-49;
	  if((x*x+y*y+z*z)<400)
	    outImage->SetPixel(pixelIndex,4);
*/
