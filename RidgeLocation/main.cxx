#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>

#define distImage(i,j,k) distImage[(k)*pRow*pCol+(j)*pRow+(i)]

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage outputImage " << std::endl;
      return EXIT_FAILURE;
    }
 
  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
    
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New(); 

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

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
 
  float *distImage;
  distImage = (float *)malloc(sizeof(float) * pRow * pCol * pSli);

  //Read image and set boundaries to 1
  ImageType::IndexType pixelIndex;
  int i, j, k;
  int a, b, c;
  float valueP, valueQ, distPQ;
  int accu;
  bool flag;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  distImage(i,j,k) = reader->GetOutput()->GetPixel(pixelIndex);
	}
  /* 
  int obsX=282;
  int obsY=353;
  int obsZ=216;
  for(i=obsX;i<obsX+1;i++)
    for(j=obsY;j<obsY+1;j++)
      for(k=obsZ;k<obsZ+1;k++)
  */
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  valueP = distImage(i,j,k);
	  flag = true;
	  accu = 0;
	  if(valueP<5)
	    flag = false;
	  else
	    {
	      for(a=i-1; a<=i+1; a++)
		for(b=j-1; b<=j+1; b++)
		  for(c=k-1; c<=k+1; c++)
		    if((a>=0)&&(a<pRow)&&(b>=0)&&(b<pCol)&&(c>=0)&&(c<pSli))
		      {
			//std::cout<<a<<" "<<b<<" "<<c<<": ";
			//distPQ = (a-i)*(a-i)+(b-j)*(b-j)+(c-k)*(c-k);
			//distPQ = sqrt(distPQ);
			valueQ = distImage(a,b,c);
			if(valueP<valueQ)
			  accu++;
			//std::cout<<valueP<<" "<<valueQ<<" "<<distPQ<<std::endl;
		      }
	      if(accu>3)
		flag = false;
	    }
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;	    
	  if(flag)
	    image->SetPixel(pixelIndex, 255);
	  else
	    image->SetPixel(pixelIndex, 0);
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
  
  free(distImage);
  return EXIT_SUCCESS;
}
