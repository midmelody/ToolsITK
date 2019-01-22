#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "Compare.h"

int main(int argc, char *argv[])
{
  int i,j; //counters

  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputFileName outputImage" << std::endl;
      return EXIT_FAILURE;
    }
 
  //Read file
  FILE *file;
  file = fopen(argv[1],"r");  
  fscanf(file,"%d\n%d\n", &pRow, &pCol);
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );

  //Allocate Image
  ImageType::SizeType size;
  size[0] = pRow;
  size[1] = pCol; 
  ImageType::RegionType region;
  region.SetSize( size );
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->Allocate();  

  //Output
  float *value;
  value = (float *)malloc(sizeof(float) * pRow * pCol);
  float temp;
  float max = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	fscanf(file,"%f\n",&temp);
	if(max < temp)
	  max = temp;
	value[(j)*pRow+(i)] = temp;
      }
  fclose(file);

  ImageType::IndexType pixelIndex;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex[0] = i;
	pixelIndex[1] = j;
	image->SetPixel(pixelIndex, int(value[(j)*pRow+(i)]/max*255));
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


