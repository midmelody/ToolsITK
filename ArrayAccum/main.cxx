#include "itkImage.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h"

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputFile1 inputFile2 binNum1 binNum2 outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  FILE * pFileX;
  long lSizeX;
  double * bufferX;
  size_t resultX;
  pFileX = fopen ( argv[1] , "rb" );
  if (pFileX==NULL) {fputs ("File error",stderr); exit (1);}
  // obtain file size:
  fseek (pFileX , 0 , SEEK_END);
  lSizeX = ftell (pFileX);
  lSizeX = lSizeX/8;
  rewind (pFileX);
  // allocate memory to contain the whole file:
  bufferX = (double*) malloc (sizeof(double)*lSizeX);
  if (bufferX == NULL) {fputs ("Memory error",stderr); exit (2);}
  // copy the file doubleo the bufferX:
  resultX = fread (bufferX,sizeof(double),lSizeX,pFileX);
  std::cout<<lSizeX<<" "<<resultX<<std::endl;
  if (resultX != lSizeX) {fputs ("Reading error",stderr); exit (3);}
  // the whole file is now loaded in the memory bufferX
  fclose (pFileX);

  FILE * pFileY;
  long lSizeY;
  double * bufferY;
  size_t resultY;
  pFileY = fopen ( argv[2] , "rb" );
  if (pFileY==NULL) {fputs ("File error",stderr); exit (1);}
  // obtain file size:
  fseek (pFileY , 0 , SEEK_END);
  lSizeY = ftell (pFileY);
  lSizeY = lSizeY/8; 
  rewind (pFileY);
  // allocate memory to contain the whole file:
  bufferY = (double*) malloc (sizeof(double)*lSizeY);
  if (bufferY == NULL) {fputs ("Memory error",stderr); exit (2);}
  // copy the file doubleo the bufferY:
  resultY = fread (bufferY,sizeof(double),lSizeY,pFileY);
  std::cout<<lSizeY<<" "<<resultY<<std::endl;
  if (resultY != lSizeY) {fputs ("Reading error",stderr); exit (3);}
  // the whole file is now loaded in the memory bufferY
  fclose (pFileY);

  //Accumulator
  int x,y;
  int binNum1 = atoi(argv[3]);
  int binNum2 = atoi(argv[4]);
  int max1 = 0;
  int max2 = 0;
  int min1 = 100000;
  int min2 = 100000;

  int i;
  for(i=0;i<lSizeY;i++)
    {
      x = int(bufferX[i]);
      y = int(bufferY[i]);
      if(max1<x) 
	max1 = x;
      if(max2<y)
	max2 = y;
      if(min1>x) 
	min1 = x;
      if(min2>y)
	min2 = y;
    }

  std::cout<<"bin1 "<<binNum1<<" max1 "<<max1<<" min1 "<<min1<<std::endl;
  std::cout<<"bin2 "<<binNum2<<" max2 "<<max2<<" min2 "<<min2<<std::endl;

  double accum[256][256];
  for(x=0;x<255;x++)
    for(y=0;y<255;y++)
      {
	accum[x][y] = 0;
      }

  for(i=0;i<lSizeY;i++)
    {
      x = int(bufferX[i]);
      y = int(bufferY[i]);
      if(x<0) x=0;
      if(y<0) y=0;
      x = round(float(x)/float(max1)*float(binNum1-1));
      y = round(float(y)/float(max2)*float(binNum2-1));
      if((x>=0)&&(x<256)&&(y>=0)&&(y<256))
	accum[x][y] = accum[x][y]+1;
      else
	{
	  std::cout<<bufferX[i]<<" "<<bufferY[i]<<std::endl;	  
	  std::cout<<x<<" "<<y<<std::endl;
	}
    }

  FILE *outFile = fopen(argv[5], "w");
  long add = 0;
  for(x=0;x<21;x++)
    for(y=0;y<256;y++)
      {
	fprintf(outFile, "%d %d %f\n", x, y, accum[x][y]);
	add = add+accum[x][y];
      }
  fclose(outFile);


  std::cout<<add<<std::endl;
  std::cout<<"--------------------------------"<<std::endl;
  free (bufferX);
  free (bufferY);

  return EXIT_SUCCESS;
}
