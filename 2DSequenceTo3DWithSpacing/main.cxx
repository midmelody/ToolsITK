#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

int main( int argc, char * argv[] )
{
  if( argc < 9 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFolder outputImageFile totalCount spacingX spacingY spacingZ originZ InverseFlag" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int InputDimension = 2;
  const unsigned int OutputDimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, InputDimension > InputImageType;
  typedef itk::Image< PixelType, OutputDimension > OutputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  
  //Parameters
  int inverseFlag = atoi(argv[8]);
  char fileNamePre[50];
  strcpy(fileNamePre, argv[1]);
  writer->SetFileName( argv[2] );

  char fileCt[5];
  int n;
  int ct = 1;
  n=sprintf (fileCt, "%02d.tif", ct);
  char fileName[60];
  strcpy(fileName, fileNamePre);
  strcat(fileName, fileCt);
  
  reader->SetFileName( fileName );
  reader->Update();
  
  //Set image specs
  OutputImageType::SpacingType spacing;
  spacing[0] = atof(argv[4]);
  spacing[1] = atof(argv[5]);
  spacing[2] = atof(argv[6]);
  OutputImageType::PointType origin; 
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = atof(argv[7]);
  OutputImageType::DirectionType direction;
  direction.SetIdentity();
  InputImageType::SizeType  sizeIn = reader->GetOutput()->GetRequestedRegion().GetSize();
  OutputImageType::SizeType  size;
  size[0] = sizeIn[0];
  size[1] = sizeIn[1];
  size[2] = atoi(argv[3]);

  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  OutputImageType::RegionType region;
  region.SetSize( size ); 

  //Allocate new image
  OutputImageType::Pointer image = OutputImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  

  //Read image 1 by 1
  InputImageType::Pointer imageTemp = InputImageType::New();
  InputImageType::IndexType pixelInIndex;
  OutputImageType::IndexType pixelOutIndex;
  int i, j, k, value;
  inverseFlag = atoi(argv[8]);
  std::cout<<inverseFlag<<std::endl;
  bool flag = inverseFlag;
  for(k=0;k<pSli;k++)
    {
      n=sprintf (fileCt, "%02d.tif", k+1);
      char fileNameTemp[50];
      strcpy(fileNameTemp, fileNamePre);
      strcat(fileNameTemp, fileCt);
      std::cout<<fileNameTemp<<std::endl;
      reader->SetFileName( fileNameTemp );
      reader->Update();
      imageTemp = reader->GetOutput();
      imageTemp->DisconnectPipeline();

      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {
	    pixelInIndex[0] = i;
	    pixelInIndex[1] = j;	    
	    pixelOutIndex[0]=i;
	    pixelOutIndex[1]=j;
	    pixelOutIndex[2]=k;

	    //value = 255- imageTemp->GetPixel(pixelInIndex);
	    //image->SetPixel(pixelOutIndex, value);
	    if(flag)
	      {
		value =  255-imageTemp->GetPixel(pixelInIndex);
		image->SetPixel(pixelOutIndex, value);
	      }
	    else
	      {
		value  = imageTemp->GetPixel(pixelInIndex);
		image->SetPixel(pixelOutIndex, value);
	      }	    
	  }
    }

  writer->SetInput( image );
  writer->Update();

  return EXIT_SUCCESS;
}
