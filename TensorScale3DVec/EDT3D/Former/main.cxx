#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "DT3D.h"

int main(int argc, char *argv[])
{
  int i,j,k,x,y,z,u,v,w,t; //counters
  int numberOfZeroCrossing; //number of zero crossings
  VOXEL *zeroCrossing; //zero crossing locations;

  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage sigmaLoGSmall gradientThreHigh sigmaLoGLarge gradientThreLow neighborhoodSize" << std::endl;
      return EXIT_FAILURE;
    }

  if(atoi(argv[6])%2!=1)
    {
      std::cerr << "neighborhood size needs to be odd" << std::endl;
      return EXIT_FAILURE;  
    }

  //Read image
  ImageReaderType::Pointer ImageReader = ImageReaderType::New();
  ImageReader->SetFileName( argv[1] );
  ImageReader->Update();
  ImageType::SpacingType spacing = ImageReader->GetOutput()->GetSpacing();
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];
  ImageType::SizeType  size = ImageReader->GetOutput()->GetRequestedRegion().GetSize();
  size[0] = size[0] + 1;
  size[1] = size[1] + 1;
  size[2] = size[2] + 1;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::RegionType region;
  ImageType::IndexType start = ImageReader->GetOutput()->GetRequestedRegion().GetIndex();
  region.SetSize( size );
  region.SetIndex( start );
  //extend input data
  float *inputData;
  inputData = (float *)malloc(sizeof(float) * pRow * pCol * pSli);  
  ImageType::IndexType pixelIndex;
  float minimum=10000;
  //Core
  for(i=0;i<pRow-1;i++)
    for(j=0;j<pCol-1;j++)
      for(k=0;k<pSli-1;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  inputData(i,j,k)=ImageReader->GetOutput()->GetPixel(pixelIndex);
	  if(minimum>inputData(i,j,k)) minimum=inputData(i,j,k);
	}
  //3 Slabs
  for(i=0;i<pRow-1;i++)
    for(j=0;j<pCol-1;j++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=j;
	pixelIndex[2]=pSli-3;
	inputData(i,j,pSli-1) = ImageReader->GetOutput()->GetPixel(pixelIndex);
      }
  for(j=0;j<pCol-1;j++)
    for(k=0;k<pSli-1;k++)
      {
	pixelIndex[0]=pRow-3;
	pixelIndex[1]=j;
	pixelIndex[2]=k;
	inputData(pRow-1,j,k) = ImageReader->GetOutput()->GetPixel(pixelIndex);
      }
  for(i=0;i<pRow-1;i++)
    for(k=0;k<pSli-1;k++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=pCol-3;
	pixelIndex[2]=k;
	inputData(i,pCol-1,k) = ImageReader->GetOutput()->GetPixel(pixelIndex);
      }
  //3 Lines
  for(i=0;i<pRow-1;i++)
    {
      pixelIndex[0]=i;
      pixelIndex[1]=pCol-3;
      pixelIndex[2]=pSli-3;  
      inputData(i,pCol-1,pSli-1) = ImageReader->GetOutput()->GetPixel(pixelIndex);
    }
  for(j=0;j<pCol-1;j++)
    {
      pixelIndex[0]=pRow-3;
      pixelIndex[1]=j;
      pixelIndex[2]=pSli-3;
      inputData(pRow-1,j,pSli-1) = ImageReader->GetOutput()->GetPixel(pixelIndex);
    }
  for(k=0;k<pSli-1;k++)
    {
      pixelIndex[0]=pRow-3;
      pixelIndex[1]=pCol-3;
      pixelIndex[2]=k;
      inputData(pRow-1,pCol-1,k) = ImageReader->GetOutput()->GetPixel(pixelIndex);
    }
  //1 Point
  pixelIndex[0]=pRow-3;
  pixelIndex[1]=pCol-3;
  pixelIndex[2]=pSli-3;
  inputData(pRow-1,pCol-1,pSli-1) = ImageReader->GetOutput()->GetPixel(pixelIndex);

  //Assign image
  ImageType::Pointer InputImage = ImageType::New();
  InputImage->SetRegions( region );
  InputImage->SetSpacing( spacing );
  InputImage->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  InputImage->SetPixel(pixelIndex,inputData(i,j,k)+absl(minimum));
	}
  free(inputData);
  /*
  ImageWriterType::Pointer inputWrt = ImageWriterType::New();
  inputWrt->SetInput( InputImage );
  inputWrt->SetFileName("input.img");
  inputWrt->Update();
  */
  /////////////////////////////////////////////////////////////////////////
  //--------------------Edge-detection Part--------------------------------
  //Assign space (largest possible) for zero crossing
  zeroCrossing = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli); 
  //Apply gradient magnitude filter
  GraFilterType::Pointer GraFilter = GraFilterType::New();
  GraFilter->SetInput( InputImage );
  GraFilter->Update();
  //Apply LoG filter (High threshold - small sigma; Low threshold -large sigma)
  LoGFilterType::Pointer LoGFilterHigh = LoGFilterType::New();
  LoGFilterHigh->SetNormalizeAcrossScale( false );
  LoGFilterHigh->SetInput( InputImage );
  const double sigmaHigh = atof( argv[2] );
  //std::cout<<sigmaHigh<<" "<<spaRow<<std::endl;
  LoGFilterHigh->SetSigma( sigmaHigh*spaRow );
  LoGFilterHigh->Update();
  LoGFilterType::Pointer LoGFilterLow = LoGFilterType::New();
  LoGFilterLow->SetNormalizeAcrossScale( false );
  LoGFilterLow->SetInput( InputImage );
  const double sigmaLow = atof( argv[4] );
  //std::cout<<sigmaLow<<" "<<spaRow<<std::endl;
  LoGFilterLow->SetSigma( sigmaLow*spaRow );
  LoGFilterLow->Update();
  //LoG response and gradient magnitude
  float *LoGImageHigh;
  LoGImageHigh = (float *)malloc(sizeof(float) * (pRow+1) * (pCol+1) * (pSli+1));  
  float *LoGImageLow;
  LoGImageLow  = (float *)malloc(sizeof(float) * (pRow+1) * (pCol+1) * (pSli+1)); 
  float *gradientImage;
  gradientImage =(float *)malloc(sizeof(float) * (pRow+1) * (pCol+1) * (pSli+1)); 
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  //Get LoG result and gradient magnitude result
	  LoGImageHigh(i,j,k) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
	  LoGImageLow(i,j,k) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
	  gradientImage(i,j,k)= GraFilter->GetOutput()->GetPixel(pixelIndex);
	} 
   //3 Slabs
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=j;
	pixelIndex[2]=pSli-2;
	LoGImageHigh(i,j,pSli) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
	LoGImageLow(i,j,pSli)= LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
	gradientImage(i,j,pSli)= GraFilter->GetOutput()->GetPixel(pixelIndex);
      }
  for(j=0;j<pCol;j++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=pRow-2;
	pixelIndex[1]=j;
	pixelIndex[2]=k;
	LoGImageHigh(pRow,j,k) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
	LoGImageLow(pRow,j,k)= LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
	gradientImage(pRow,j,k)= GraFilter->GetOutput()->GetPixel(pixelIndex);
      }
  for(i=0;i<pRow;i++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=pCol-2;
	pixelIndex[2]=k;
	LoGImageHigh(i,pCol,k) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
	LoGImageLow(i,pCol,k)= LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
	gradientImage(i,pCol,k)= GraFilter->GetOutput()->GetPixel(pixelIndex);
      }
  //3 Lines
  for(i=0;i<pRow;i++)
    {
      pixelIndex[0]=i;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=pSli-2;  
      LoGImageHigh(i,pCol,pSli) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
      LoGImageLow(i,pCol,pSli)= LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
      gradientImage(i,pCol,pSli)= GraFilter->GetOutput()->GetPixel(pixelIndex);
    }
  for(j=0;j<pCol;j++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=j;
      pixelIndex[2]=pSli-2;
      LoGImageHigh(pRow,j,pSli) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
      LoGImageLow(pRow,j,pSli)= LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
      gradientImage(pRow,j,pSli)= GraFilter->GetOutput()->GetPixel(pixelIndex);
    }
  for(k=0;k<pSli;k++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=k;
      LoGImageHigh(pRow,pCol,k) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
      LoGImageLow(pRow,pCol,k)= LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
      gradientImage(pRow,pCol,k)= GraFilter->GetOutput()->GetPixel(pixelIndex);
    }
  //1 Point
  pixelIndex[0]=pRow-2;
  pixelIndex[1]=pCol-2;
  pixelIndex[2]=pSli-2;
  LoGImageHigh(pRow,pCol,pSli) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
  LoGImageLow(pRow,pCol,pSli) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
  gradientImage(pRow,pCol,pSli) = GraFilter->GetOutput()->GetPixel(pixelIndex);
 
  //write output for LoG and gradient
  /*
  ImageType::SizeType  size1;
  size1[0] = size[0] + 1;
  size1[1] = size[1] + 1;
  size1[2] = size[2] + 1;
  ImageType::RegionType region1;
  region1.SetSize( size1 );
  ImageType::Pointer LoGHImage = ImageType::New();
  LoGHImage->SetRegions( region1 );
  LoGHImage->SetSpacing( spacing );
  LoGHImage->Allocate();
  ImageType::Pointer LoGLImage = ImageType::New();
  LoGLImage->SetRegions( region1 );
  LoGLImage->SetSpacing( spacing );
  LoGLImage->Allocate();
  ImageType::Pointer GraImage = ImageType::New();
  GraImage->SetRegions( region1 );
  GraImage->SetSpacing( spacing );
  GraImage->Allocate();
  for(i=0;i<=pRow;i++)
    for(j=0;j<=pCol;j++)
      for(k=0;k<=pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  LoGHImage->SetPixel(pixelIndex,LoGImageHigh(i,j,k));
	  LoGLImage->SetPixel(pixelIndex,LoGImageLow(i,j,k));
	  GraImage->SetPixel(pixelIndex,gradientImage(i,j,k));  
	}
  //Writes result for LoGLow, LoGHigh, and gradient
  ImageWriterType::Pointer writerLoGH = ImageWriterType::New();
  writerLoGH->SetInput( LoGHImage );
  writerLoGH->SetFileName("LoGH.img");
  writerLoGH->Update();
  ImageWriterType::Pointer writerLoGL = ImageWriterType::New();
  writerLoGL->SetInput( LoGLImage );  
  writerLoGL->SetFileName("LoGL.img");
  writerLoGL->Update();
  ImageWriterType::Pointer writerGra = ImageWriterType::New();
  writerGra->SetInput( GraImage );  
  writerGra->SetFileName("Gra.img");
  writerGra->Update();
  */

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Compute zero crossing
  const float graThreHigh = atof( argv[3] );
  const float graThreLow  = atof( argv[5] );
  std::cout<<"Start computing zero crossings: "<<std::endl;
  numberOfZeroCrossing = ComputeZeroCrossing(LoGImageHigh,LoGImageLow,gradientImage,graThreHigh,graThreLow,zeroCrossing);
  std::cout<<"Total number of zero crossings: "<<numberOfZeroCrossing<<std::endl;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //--------------------Distance Transform Part----------------------------
  //Compute distance transform 
  float *DT;
  DT = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  VOXEL *nearestEdgePoint;
  nearestEdgePoint =(VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);    
  ImageType::Pointer DTImage = ImageType::New();
  DTImage->SetRegions( region );
  DTImage->SetSpacing( spacing );
  DTImage->Allocate();
  //Initialize
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  DT(i,j,k) = 20000;
	  nearestEdgePoint(i,j,k).x = -1;
	  nearestEdgePoint(i,j,k).y = -1;
	  nearestEdgePoint(i,j,k).z = -1;
	}

  //Compute raster scan

  std::cout<<"Start raster scan EDT: "<<std::endl;
  const int neighborCount =  atoi( argv[6] );
  RasterScanEDT(numberOfZeroCrossing, zeroCrossing, DT, nearestEdgePoint, neighborCount); 
  //write output for DT
  //compute point-NEP vector orientation
  float maxDT = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  if(maxDT<DT(i,j,k)) maxDT=DT(i,j,k);
	  DTImage->SetPixel(pixelIndex,sqrt(sqrt(DT(i,j,k))));
	}
  //output DT
  ImageWriterType::Pointer writerDT = ImageWriterType::New();
  writerDT->SetInput( DTImage );
  writerDT->SetFileName("temp/DT.img");
  writerDT->Update();
 
  free(LoGImageHigh);
  free(LoGImageLow);
  free(gradientImage);
  free(DT);
  free(nearestEdgePoint);
  return EXIT_SUCCESS;
}
