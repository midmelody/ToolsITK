#include <iostream>
#include <math.h>
#include <string>
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//------------------ITK Part-----------------------------------------------
typedef unsigned char PixelType; 
typedef itk::Image< PixelType, 2 > ImageType2D;
typedef itk::Image< PixelType, 3 > ImageType3D;
typedef itk::ImageFileReader< ImageType2D > ReaderType;
typedef itk::ImageFileWriter< ImageType3D > WriterType;
//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image: row, column, slice
//------------------Coordinate convert-------------------------------------
#define data(i,j,k) data[(k)*pRow*pCol+(j)*pRow+(i)]

int main(int argc, char *argv[])
{
  int i, j, k; //counters

  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImagePrefix outputImage sliceCount" << std::endl;
      return EXIT_FAILURE;
    }
 
  //Read image
  char fileNameRef[50];    
  sprintf(fileNameRef, "%s_0.png", argv[1]);
  std::cout<<fileNameRef<<std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( fileNameRef );
  reader->Update();
  
  //Get image specs
  ImageType2D::SizeType size = reader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = atoi(argv[3]);
  std::cout<<pRow<<" "<<pCol<<" "<<pSli<<std::endl;
  
  //Get info
  ImageType2D::IndexType index2D;
  std::cout<<"here"<<std::endl;
  unsigned char *data;
  data = new unsigned char[pRow*pCol*pSli];
  std::cout<<"here"<<std::endl;
  
  for(k=0;k<pSli;k++)
    {
      std::cout<<k<<std::endl;
      char fileName[50];    
      sprintf(fileName, "%s_%d.png", argv[1], k);
      std::cout<<fileName<<std::endl;
      reader->SetFileName( fileName );
      reader->Update();
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {
	    index2D[0] = i;
	    index2D[1] = j;
	    data(i,j,k) = reader->GetOutput()->GetPixel(index2D);
	  }  
    }

  //Set value;
  ImageType3D::SizeType size3D;
  size3D[0] = pRow;
  size3D[1] = pCol;
  size3D[2] = pSli;
  ImageType3D::RegionType region;
  region.SetSize( size3D );
  ImageType3D::Pointer outImage = ImageType3D::New();
  outImage->SetRegions( region );
  outImage->Allocate(); 
  ImageType3D::IndexType index3D;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  index3D[0]=i;
	  index3D[1]=j;
	  index3D[2]=k;
	  outImage->SetPixel(index3D, data(i,j,k));
	}

  //Write
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( outImage );
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

  delete [] data;
  return EXIT_SUCCESS;
}


