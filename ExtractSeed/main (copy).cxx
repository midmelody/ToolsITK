#include "itkImage.h"
#include "itkImageFileReader.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " seedImageFile seedTxtFile seedLabelThreshold(use 0 for all seeds)" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  int seedLabel = atoi(argv[3]);  

  //Pipeline
  try
    {
      reader->Update();
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
  int pRow = size[0];
  int pCol = size[1];
  int pSli = size[2]; 
  
  int i,j,k,value;
  ImageType::IndexType pixelIndex;
  FILE *seedFile = fopen(argv[2], "w");
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  if(value)
	    {
	      if(seedLabel==0)
		fprintf(seedFile, "%d %d %d 1\n", i, j, k) ;
	      else
		{
		  if(value <= seedLabel)
		    fprintf(seedFile, "%d %d %d %d\n", i, j, k, value) ;		    
		}
	    }
	}
  fclose(seedFile);
  return EXIT_SUCCESS;
}
