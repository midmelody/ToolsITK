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
      std::cerr << argv[0] << " inputImageFile maskImageFile seedFile " << std::endl;
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
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = size[0];
  int pCol = size[1];
  int pSli = size[2]; 
  
  int i,j,k,value,valueM;
  ImageType::IndexType pixelIndex;
  int max = 0;
  int iM, jM, kM;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  valueM = readerM->GetOutput()->GetPixel(pixelIndex);
	  if(valueM)
	    {
	      value = reader->GetOutput()->GetPixel(pixelIndex);
	      if(max<value)
		{
		  max = value;
		  iM = i;
		  jM = j;
		  kM = k;
		}
	    }
	}
  
  FILE *seedFile = fopen(argv[3], "w");
  fprintf(seedFile, "%d %d %d \n", iM, jM, kM) ;
  fclose(seedFile);
  
  return EXIT_SUCCESS;
}
