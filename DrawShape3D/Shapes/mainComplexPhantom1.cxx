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
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  outImage->SetPixel(pixelIndex,20);	  
	}

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
	  //Sin wall Group1
	  for(j=2;j<7;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*14*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	  for(j=15;j<20;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*14*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	  for(j=28;j<33;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*14*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	  //Sin wall Group2
	  for(j=89;j<97;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*12*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	  for(j=107;j<115;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*12*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	  for(j=125;j<134;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*12*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	  //Sin wall Group3
	  for(j=190;j<200;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*10*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	  for(j=210;j<220;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*10*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	  for(j=230;j<240;j++)
	    {
	      b = round(mag*sin(float(i)/sizeF*10*PI)+mag+j);
	      pixelIndex[0]=i;
	      pixelIndex[1]=b;
	      pixelIndex[2]=k;
	      outImage->SetPixel(pixelIndex,100);
	    }
	}

      //Circle sin tube
      int centerA0,centerA1,centerB0, centerB10, centerB11;    
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {	
	    x = i;
	    y = j;
	    centerA0 = round(mag*0.6*sin(float(k)/sizeF*8*PI)+415);
	    centerB0 = 455;
	    if(((x-centerA0)*(x-centerA0)+(y-centerB0)*(y-centerB0))<400)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerA1 = round(mag*0.6*sin(float(k)/sizeF*8*PI)+450);
	    centerB10 = 480;
	    centerB11 = 455;
	    if(((x-centerA1)*(x-centerA1)+(y-centerB10)*(y-centerB10))<100)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    if(((x-centerA1)*(x-centerA1)+(y-centerB11)*(y-centerB11))<100)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerB11 = 430;
	    if(((x-centerA1)*(x-centerA1)+(y-centerB11)*(y-centerB11))<100)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }


	    centerA1 = round(mag*0.6*sin(float(k)/sizeF*8*PI)+470);
	    centerB10 = 485;
	    centerB0 = 470;
	    centerB11 = 455;
	    if(((x-centerA1)*(x-centerA1)+(y-centerB0)*(y-centerB0))<25)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    if(((x-centerA1)*(x-centerA1)+(y-centerB10)*(y-centerB10))<25)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    if(((x-centerA1)*(x-centerA1)+(y-centerB11)*(y-centerB11))<25)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerB10 = 440;
	    centerB11 = 425;
	    if(((x-centerA1)*(x-centerA1)+(y-centerB10)*(y-centerB10))<25)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    if(((x-centerA1)*(x-centerA1)+(y-centerB11)*(y-centerB11))<25)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	  }

      //Ellipse sin tube
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {	
	    x = i;
	    y = j;
	    centerA0 = 70;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+310);
	    if((float((x-centerA0)*(x-centerA0))/225+float((y-centerB0)*(y-centerB0))/900)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerA0 = 25;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+360);
	    if((float((x-centerA0)*(x-centerA0))/400+float((y-centerB0)*(y-centerB0))/100)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerA0 = 70;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+360);
	    if((float((x-centerA0)*(x-centerA0))/400+float((y-centerB0)*(y-centerB0))/100)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerA0 = 115;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+360);
	    if((float((x-centerA0)*(x-centerA0))/400+float((y-centerB0)*(y-centerB0))/100)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }

	    centerA0 = 30;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+390);
	    if((float((x-centerA0)*(x-centerA0))/25+float((y-centerB0)*(y-centerB0))/100)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerA0 = 50;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+390);
	    if((float((x-centerA0)*(x-centerA0))/25+float((y-centerB0)*(y-centerB0))/100)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerA0 = 70;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+390);
	    if((float((x-centerA0)*(x-centerA0))/25+float((y-centerB0)*(y-centerB0))/100)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerA0 = 90;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+390);
	    if((float((x-centerA0)*(x-centerA0))/25+float((y-centerB0)*(y-centerB0))/100)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	    centerA0 = 110;
	    centerB0 = round(mag/2*sin(float(k)/sizeF*12*PI)+390);
	    if((float((x-centerA0)*(x-centerA0))/25+float((y-centerB0)*(y-centerB0))/100)<1)
	      {
		a = round(x);
		b = round(y);
		c = k;
		pixelIndex[0]=a;
		pixelIndex[1]=b;
		pixelIndex[2]=c;
		outImage->SetPixel(pixelIndex,100);
	      }
	  }

      //Circle walls
      centerA0 = 210;
      centerB0 = 430;
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {	
	    x = float(i);
	    y = float(j);
	    float length = (x-centerA0)*(x-centerA0)+(y-centerB0)*(y-centerB0);
	    if((length<15*15)||((length>25*25)&&(length<35*35))||((length>42*42)&&(length<50*50))||((length>55*55)&&(length<60*60)))
	      {
		pixelIndex[0]=i;
		pixelIndex[1]=j;
		pixelIndex[2]=k;
		outImage->SetPixel(pixelIndex,100);		
	      }
	  }

      //Block walls
      centerA0 = 360;
      centerB0 = 360;
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {	
	    x = float(i-centerA0)*sqrt(2)/2-float(j-centerB0)*sqrt(2)/2;
	    y = float(i-centerA0)*sqrt(2)/2+float(j-centerB0)*sqrt(2)/2;
	    float length;
	    if(fabs(x)>fabs(y))
	      length = fabs(x);
	    else
	      length = fabs(y);
	    if((length<15)||((length>25)&&(length<35))||((length>42)&&(length<50))||((length>55)&&(length<60)))
	      {
		pixelIndex[0]=i;
		pixelIndex[1]=j;
		pixelIndex[2]=k;
		outImage->SetPixel(pixelIndex,100);		
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
