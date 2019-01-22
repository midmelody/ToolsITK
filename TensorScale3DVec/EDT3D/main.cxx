#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "DT3D.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters
  int numberOfZeroCrossing; //number of zero crossings
  VOXEL *zeroCrossing; //zero crossing locations;

  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage sigmaLoGSmall gradientThreHigh sigmaLoGLarge gradientThreLow RasterNeighborhoodSize" << std::endl;
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
  //Core
  for(i=0;i<pRow-1;i++)
    for(j=0;j<pCol-1;j++)
      for(k=0;k<pSli-1;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  inputData(i,j,k)=ImageReader->GetOutput()->GetPixel(pixelIndex);
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
  //ImageReader->Delete();
  
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
	  InputImage->SetPixel(pixelIndex,inputData(i,j,k));
	}
  free(inputData);

  /////////////////////////////////////////////////////////////////////////
  //--------------------Edge-detection Part--------------------------------
  const double sigmaHigh = atof( argv[2] )*spaRow;
  const double sigmaLow = atof( argv[4] )*spaRow;
  float *DoGImageHigh;
  DoGImageHigh =(float *)malloc(sizeof(float) * (pRow+1) * (pCol+1) * (pSli+1)); 
  float *LoGImageHigh;
  LoGImageHigh = (float *)malloc(sizeof(float) * (pRow+1) * (pCol+1) * (pSli+1)); 
  float *DoGImageLow;
  DoGImageLow =(float *)malloc(sizeof(float) * (pRow+1) * (pCol+1) * (pSli+1)); 
  float *LoGImageLow;
  LoGImageLow = (float *)malloc(sizeof(float) * (pRow+1) * (pCol+1) * (pSli+1)); 

  //---------------High------------------
  //------Apply DoG filter--------
  DevFilterType::Pointer DoGFilterHighX = DevFilterType::New();
  DevFilterType::Pointer DoGFilterHighY = DevFilterType::New();
  DevFilterType::Pointer DoGFilterHighZ = DevFilterType::New();
  DoGFilterHighX->SetOrder( DevFilterType::ZeroOrder );
  DoGFilterHighY->SetOrder( DevFilterType::ZeroOrder );
  DoGFilterHighZ->SetOrder( DevFilterType::ZeroOrder );  
  DoGFilterHighX->SetDirection( 0 );
  DoGFilterHighY->SetDirection( 1 );
  DoGFilterHighZ->SetDirection( 2 );
  DoGFilterHighX->SetSigma( sigmaHigh );
  DoGFilterHighY->SetSigma( sigmaHigh );
  DoGFilterHighZ->SetSigma( sigmaHigh );
  DoGFilterHighX->SetInput( InputImage );  
  DoGFilterHighY->SetInput( DoGFilterHighX->GetOutput() );
  DoGFilterHighZ->SetInput( DoGFilterHighY->GetOutput() );
  GraFilterType::Pointer GraFilterHigh = GraFilterType::New();
  GraFilterHigh->SetInput( DoGFilterHighZ->GetOutput() );
  GraFilterHigh->Update();

  //Core
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i; 
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  DoGImageHigh(i,j,k) = GraFilterHigh->GetOutput()->GetPixel(pixelIndex);
	} 

  //3 Slabs
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=j;
	pixelIndex[2]=pSli-2;
	DoGImageHigh(i,j,pSli) = GraFilterHigh->GetOutput()->GetPixel(pixelIndex);
      }
  for(j=0;j<pCol;j++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=pRow-2;
	pixelIndex[1]=j;
	pixelIndex[2]=k;
	DoGImageHigh(pRow,j,k) = GraFilterHigh->GetOutput()->GetPixel(pixelIndex);
      }
  for(i=0;i<pRow;i++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=pCol-2;
	pixelIndex[2]=k;
	DoGImageHigh(i,pCol,k) = GraFilterHigh->GetOutput()->GetPixel(pixelIndex);
      }

  //3 Lines
  for(i=0;i<pRow;i++)
    {
      pixelIndex[0]=i;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=pSli-2;  
      DoGImageHigh(i,pCol,pSli) = GraFilterHigh->GetOutput()->GetPixel(pixelIndex);
    }
  for(j=0;j<pCol;j++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=j;
      pixelIndex[2]=pSli-2;
      DoGImageHigh(pRow,j,pSli) = GraFilterHigh->GetOutput()->GetPixel(pixelIndex);
    }
  for(k=0;k<pSli;k++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=k;
      DoGImageHigh(pRow,pCol,k) = GraFilterHigh->GetOutput()->GetPixel(pixelIndex);
    }

  //1 Point
  pixelIndex[0]=pRow-2;
  pixelIndex[1]=pCol-2;
  pixelIndex[2]=pSli-2;
  DoGImageHigh(pRow,pCol,pSli) = GraFilterHigh->GetOutput()->GetPixel(pixelIndex);

  //free space
  //DoGFilterHighX->Delete(); 
  //DoGFilterHighY->Delete();
  //DoGFilterHighZ->Delete();
  //GraFilterHigh->Delete();
 
  //------Apply LoG filter--------
  LoGFilterType::Pointer LoGFilterHigh = LoGFilterType::New();
  LoGFilterHigh->SetNormalizeAcrossScale( false );
  LoGFilterHigh->SetInput( InputImage );
  LoGFilterHigh->SetSigma( sigmaHigh );
  LoGFilterHigh->Update();
  //Core
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  LoGImageHigh(i,j,k) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
	} 
  //3 Slabs
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=j;
	pixelIndex[2]=pSli-2;
	LoGImageHigh(i,j,pSli) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
      }
  for(j=0;j<pCol;j++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=pRow-2;
	pixelIndex[1]=j;
	pixelIndex[2]=k;
	LoGImageHigh(pRow,j,k) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
      }
  for(i=0;i<pRow;i++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=pCol-2;
	pixelIndex[2]=k;
	LoGImageHigh(i,pCol,k) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
      }
  //3 Lines
  for(i=0;i<pRow;i++)
    {
      pixelIndex[0]=i;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=pSli-2;  
      LoGImageHigh(i,pCol,pSli) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
    }
  for(j=0;j<pCol;j++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=j;
      pixelIndex[2]=pSli-2;
      LoGImageHigh(pRow,j,pSli) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
    }
  for(k=0;k<pSli;k++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=k;
      LoGImageHigh(pRow,pCol,k) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
    }
  //1 Point
  pixelIndex[0]=pRow-2;
  pixelIndex[1]=pCol-2;
  pixelIndex[2]=pSli-2;
  LoGImageHigh(pRow,pCol,pSli) = LoGFilterHigh->GetOutput()->GetPixel(pixelIndex);
  //free space
  //LoGFilterHigh->Delete(); 

 
  //---------------Low------------------
  //------Apply DoG filter--------
  DevFilterType::Pointer DoGFilterLowX = DevFilterType::New();
  DevFilterType::Pointer DoGFilterLowY = DevFilterType::New();
  DevFilterType::Pointer DoGFilterLowZ = DevFilterType::New();
  DoGFilterLowX->SetOrder( DevFilterType::ZeroOrder );
  DoGFilterLowY->SetOrder( DevFilterType::ZeroOrder );
  DoGFilterLowZ->SetOrder( DevFilterType::ZeroOrder );  
  DoGFilterLowX->SetDirection( 0 );
  DoGFilterLowY->SetDirection( 1 );
  DoGFilterLowZ->SetDirection( 2 );
  DoGFilterLowX->SetSigma( sigmaLow );
  DoGFilterLowY->SetSigma( sigmaLow );
  DoGFilterLowZ->SetSigma( sigmaLow );
  DoGFilterLowX->SetInput( InputImage );  
  DoGFilterLowY->SetInput( DoGFilterLowX->GetOutput() );
  DoGFilterLowZ->SetInput( DoGFilterLowY->GetOutput() );
  GraFilterType::Pointer GraFilterLow = GraFilterType::New();
  GraFilterLow->SetInput( DoGFilterLowZ->GetOutput() );
  GraFilterLow->Update();
  //Core
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  DoGImageLow(i,j,k) = GraFilterLow->GetOutput()->GetPixel(pixelIndex);
	} 
  //3 Slabs
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=j;
	pixelIndex[2]=pSli-2;
	DoGImageLow(i,j,pSli) = GraFilterLow->GetOutput()->GetPixel(pixelIndex);
      }
  for(j=0;j<pCol;j++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=pRow-2;
	pixelIndex[1]=j;
	pixelIndex[2]=k;
	DoGImageLow(pRow,j,k) = GraFilterLow->GetOutput()->GetPixel(pixelIndex);
      }
  for(i=0;i<pRow;i++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=pCol-2;
	pixelIndex[2]=k;
	DoGImageLow(i,pCol,k) = GraFilterLow->GetOutput()->GetPixel(pixelIndex);
      }
  //3 Lines
  for(i=0;i<pRow;i++)
    {
      pixelIndex[0]=i;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=pSli-2;  
      DoGImageLow(i,pCol,pSli) = GraFilterLow->GetOutput()->GetPixel(pixelIndex);
    }
  for(j=0;j<pCol;j++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=j;
      pixelIndex[2]=pSli-2;
      DoGImageLow(pRow,j,pSli) = GraFilterLow->GetOutput()->GetPixel(pixelIndex);
    }
  for(k=0;k<pSli;k++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=k;
      DoGImageLow(pRow,pCol,k) = GraFilterLow->GetOutput()->GetPixel(pixelIndex);
    }
  //1 Point
  pixelIndex[0]=pRow-2;
  pixelIndex[1]=pCol-2;
  pixelIndex[2]=pSli-2;
  DoGImageLow(pRow,pCol,pSli) = GraFilterLow->GetOutput()->GetPixel(pixelIndex);
  //free space
  //DoGFilterLowX->Delete(); 
  //DoGFilterLowY->Delete();
  //DoGFilterLowZ->Delete();
  //GraFilterLow->Delete();

  //------Apply LoG filter--------
  LoGFilterType::Pointer LoGFilterLow = LoGFilterType::New();
  LoGFilterLow->SetNormalizeAcrossScale( false );
  LoGFilterLow->SetInput( InputImage );
  LoGFilterLow->SetSigma( sigmaLow );
  LoGFilterLow->Update();
  //Core
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  LoGImageLow(i,j,k) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
	} 
  //3 Slabs
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=j;
	pixelIndex[2]=pSli-2;
	LoGImageLow(i,j,pSli) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
      }
  for(j=0;j<pCol;j++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=pRow-2;
	pixelIndex[1]=j;
	pixelIndex[2]=k;
	LoGImageLow(pRow,j,k) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
      }
  for(i=0;i<pRow;i++)
    for(k=0;k<pSli;k++)
      {
	pixelIndex[0]=i;
	pixelIndex[1]=pCol-2;
	pixelIndex[2]=k;
	LoGImageLow(i,pCol,k) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
      }
  //3 Lines
  for(i=0;i<pRow;i++)
    {
      pixelIndex[0]=i;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=pSli-2;  
      LoGImageLow(i,pCol,pSli) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
    }
  for(j=0;j<pCol;j++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=j;
      pixelIndex[2]=pSli-2;
      LoGImageLow(pRow,j,pSli) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
    }
  for(k=0;k<pSli;k++)
    {
      pixelIndex[0]=pRow-2;
      pixelIndex[1]=pCol-2;
      pixelIndex[2]=k;
      LoGImageLow(pRow,pCol,k) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
    }
  //1 Point
  pixelIndex[0]=pRow-2;
  pixelIndex[1]=pCol-2;
  pixelIndex[2]=pSli-2;
  LoGImageLow(pRow,pCol,pSli) = LoGFilterLow->GetOutput()->GetPixel(pixelIndex);
  //free space
  //LoGFilterLow->Delete();
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Compute zero crossing
  //Assign space (largest possible) for zero crossing
  zeroCrossing = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli); 
  const float graThreHigh = atof( argv[3] );
  const float graThreLow  = atof( argv[5] );
  std::cout<<"Start computing zero crossings: "<<std::endl;
  numberOfZeroCrossing = ComputeZeroCrossing(LoGImageHigh,LoGImageLow,DoGImageHigh,DoGImageLow,graThreHigh,graThreLow,zeroCrossing);
  std::cout<<"Total number of zero crossings: "<<numberOfZeroCrossing<<std::endl;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  free(LoGImageHigh);
  free(LoGImageLow);
  free(DoGImageHigh);
  free(DoGImageLow);
 
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
	  DTImage->SetPixel(pixelIndex,DT(i,j,k));
	}

  //output local scale
  LocMaxFilterType::Pointer LocMaxFilter = LocMaxFilterType::New();
  LocMaxFilter->SetInput(DTImage);
  LocMaxFilter->SetHeight(0.1);
  LocMaxFilter->Update();
  ImageWriterType::Pointer writerLocMax = ImageWriterType::New();
  writerLocMax->SetInput(LocMaxFilter->GetOutput() );
  writerLocMax->SetFileName("temp/LocMax.img");
  writerLocMax->Update();


  //output DT
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
  ImageWriterType::Pointer writerDT = ImageWriterType::New();
  writerDT->SetInput( DTImage );
  writerDT->SetFileName("temp/DT.img");
  writerDT->Update();




  int a,b,c;
  float r; //update radius length
  int currentX, currentY, currentZ; //float index of pixel 
  float *localSize;
  localSize = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  localSize(i,j,k)=0;
	}

  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  if( LocMaxFilter->GetOutput()->GetPixel(pixelIndex))
	    {
	      //std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	      r = sqrt(DT(i,j,k));
	      //std::cout<<r<<std::endl;
	      for(a=floor(-r);a<=ceil(r);a++)
		for(b=floor(-r);b<=ceil(r);b++)
		  for(c=floor(-r);c<=ceil(r);c++)	      
		    {
		      currentX = i+a;
		      currentY = j+b;
		      currentZ = k+c;		  
		      if((currentX>=0)&&(currentX<pRow)&&(currentY>=0)&&(currentY<pCol)&&(currentZ>=0)&&(currentZ<pSli))
			if((a*a+b*b+c*c)<=(r+1)*(r+1))
			  if(r>localSize(currentX,currentY,currentZ))
			    localSize(currentX,currentY,currentZ) = r;
		    }
	    }
	}
  
  ImageType::Pointer LocSizeImage = ImageType::New();
  LocSizeImage->SetRegions( region );
  LocSizeImage->SetSpacing( spacing );
  LocSizeImage->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  LocSizeImage->SetPixel(pixelIndex,localSize(i,j,k));
	}
  //output local size
  ImageWriterType::Pointer writerLocSize = ImageWriterType::New();
  writerLocSize->SetInput( LocSizeImage );
  writerLocSize->SetFileName("temp/LocSize.img");
  writerLocSize->Update();


  free(localSize);
  free(DT);
  free(nearestEdgePoint);
  return EXIT_SUCCESS;
}
