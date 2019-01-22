#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//------------------Coordinate convert-------------------------------------
#define localSize(i,j,k) localSize[(k)*pRow*pCol+(j)*pRow+(i)]

int main(int argc, char *argv[])
{
  int i,j,k; //Global counters
  int a,b,c; //Local counters
  int value; //DT value
  int mask; //Local Maxima mask
  if( argc < 4 ) 
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " DTImage LocalMaxImage outputImage" << std::endl;
      return EXIT_FAILURE;
    }

  //------------------ITK Part-----------------------------------------------
  const unsigned int Dimension = 3;
  typedef float PixelType; 
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;

  //Read image
  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  ImageReaderType::Pointer readerM = ImageReaderType::New();
  readerM->SetFileName( argv[2] );
  readerM->Update();
  
  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );

  //Allocate new image
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate(); 

  //Get info and 
  ImageType::IndexType pixelIndex;
  float *localSize;
  localSize = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	localSize(i,j,k) = 0;
 
  std::cout<<"Find local structure size: "<<std::endl;
  for(i=0;i<pRow;i++)
    {   
      std::cout<<"\r";
      std::cout<<i+1<<"/"<<pRow<<std::flush;
      for(j=0;j<pCol;j++)
     	for(k=0;k<pSli;k++) 
	  {
	    pixelIndex[0]=i;
	    pixelIndex[1]=j;
	    pixelIndex[2]=k;
	    value = reader->GetOutput()->GetPixel(pixelIndex);
	    mask = readerM->GetOutput()->GetPixel(pixelIndex);
	    if(mask)
	      {
		int search = ceil(value);
		float dist;
		for(a=i-search;a<i+search;a++)
		  for(b=j-search;b<j+search;b++)
		    for(c=k-search;c<k+search;c++)
		      if((a>=0)&&(a<pRow)&&(b>=0)&&(b<pCol)&&(c>=0)&&(c<pSli))
			{ 
			  dist = (a-i)*(a-i)+(b-j)*(b-j)+(c-k)*(c-k);		    
			  if((dist<search*search)&&(localSize(a,b,c)<value))
			    localSize(a,b,c) = value;
			}		  
	      }
	  }
      }

  std::cout<<std::endl;

  //Output structure size
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  image->SetPixel(pixelIndex,localSize(i,j,k));
      }
  
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( image );
  writer->SetFileName( argv[3] );
  writer->Update();
  
  free(localSize);
  return EXIT_SUCCESS;
}
