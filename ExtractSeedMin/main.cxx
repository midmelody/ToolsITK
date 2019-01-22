#include "itkImage.h"
#include "itkImageFileReader.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile maskImageFile seedFile sliceRangeMin sliceRangeMax" << std::endl;
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
  ImageType::SizeType  sizeM = reader->GetOutput()->GetRequestedRegion().GetSize();
  if((pRow != sizeM[0]) || (pCol != sizeM[1]) || (pSli != sizeM[2]))
    {
      std::cout<<"Mask dimension doesn't match!"<<std::endl;
      return EXIT_FAILURE;
    }
  int i,j,k,value,valueM;
  ImageType::IndexType pixelIndex;
  int min = 1000000;
  int iM, jM, kM;
  int sliceRangeMin = atoi(argv[4]);
  int sliceRangeMax = atoi(argv[5]);
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=sliceRangeMin-1;k<sliceRangeMax;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  valueM = readerM->GetOutput()->GetPixel(pixelIndex);
	  if(valueM>0)
	    {
	      value = reader->GetOutput()->GetPixel(pixelIndex);
	      if(min>value)
		{
		  min = value;
		  iM = i;
		  jM = j;
		  kM = k;
		}
	    }
	}
  std::cout<<min<<std::endl;
  FILE *seedFile = fopen(argv[3], "w");
  fprintf(seedFile, "%d %d %d \n", iM, jM, kM) ;
  fclose(seedFile);
  
  return EXIT_SUCCESS;
}
