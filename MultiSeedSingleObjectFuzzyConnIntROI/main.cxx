//The program compute multi-seed Iterative Rlative Fuzzy Connectedness
#include "FuzzyCon.h"

int main( int argc, char * argv[] )
{
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile ROIImageFile ObjectMaskImageFile outputImageFile sigma  mean" << std::endl;
      return EXIT_FAILURE;
    }
 
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerROI = ReaderType::New();
  ReaderType::Pointer readerMSK = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  readerROI->SetFileName( argv[2] );
  readerMSK->SetFileName( argv[3] );
  writer->SetFileName( argv[4] );
  float sigma = atof(argv[5]);
  float mean = atoi(argv[6]);
 
  //Pipeline
  std::cout<<std::endl<<"Reading images......."<<std::flush;
  try
    {
      reader->Update();
      readerROI->Update();
      readerMSK->Update();
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
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  OutImageType::Pointer image = OutImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->Allocate();  

  //Read image
  float *inputImage;
  inputImage = ( float *)malloc(sizeof( float ) * pRow * pCol * pSli); 
  bool *ROIImage, *MSKImage;
  ROIImage = (bool *)malloc(sizeof(bool) * pRow * pCol * pSli); 
  MSKImage = (bool *)malloc(sizeof(bool) * pRow * pCol * pSli);

  ImageType::IndexType pixelIndex;
  int i, j, k, seedCount = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  inputImage(i,j,k) = reader->GetOutput()->GetPixel(pixelIndex);
	  ROIImage(i,j,k) = bool(readerROI->GetOutput()->GetPixel(pixelIndex));
	  MSKImage(i,j,k) = bool(readerMSK->GetOutput()->GetPixel(pixelIndex));
	  if((ROIImage(i,j,k))&&(MSKImage(i,j,k))&&(inputImage(i,j,k)>-900))
	    seedCount++;
	}
  std::cout<<"finished"<<std::endl;
  std::cout<<std::endl;
  
  //Read seed
  int index;
  int seed[seedCount];
  int ct = 0;
  std::cout<<"Generating "<<seedCount<<" seeds......."<<std::flush;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  if((ROIImage(i,j,k))&&(MSKImage(i,j,k))&&(inputImage(i,j,k)>-900))
	    {
	      index = k*pRow*pCol+j*pRow+i;
	      seed[ct] = index;
	      ct++;
	    }
	}
  std::cout<<ct<<" finished "<<std::endl;
  std::cout<<std::endl;
  
  //multi-seed single-object kFOEMS
  //Compute fuzzy connectedness
  float *fuzzyConImage;
  fuzzyConImage = ( float *)malloc(sizeof( float) * pRow * pCol * pSli); 
  memset(fuzzyConImage,0,pRow*pCol*pSli*sizeof( float));
  //Compute use kFOEMS
  std::cout <<"Computing kFOEMS......."<<std::flush;
  kFOEMS( inputImage, spacing, sigma, mean, seed, seedCount, fuzzyConImage, ROIImage);
  std::cout <<"finished"<<std::endl;
  std::cout<<std::endl;
  //Set output
  std::cout <<"Writing output......."<<std::flush;
  float value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = float(fuzzyConImage(i,j,k));
	  image->SetPixel(pixelIndex, value);
	}
  //Write output  
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
  std::cout <<"finished"<<std::endl;
  std::cout<<std::endl;
  return EXIT_SUCCESS;
}
