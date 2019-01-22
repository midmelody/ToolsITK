#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "GraEDT3D.h"

int main(int argc, char *argv[])
{
  int i,j,k,a,b,t; //counters
 
  if( argc < 1 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " neighborhoodSize " << std::endl;
      return EXIT_FAILURE;
    }

  //Read image
  ImageReaderType::Pointer ImageReader = ImageReaderType::New();
  ImageReader->SetFileName( "temp/DT.img" );
  ImageReader->Update();
  ImageType::SpacingType spacing = ImageReader->GetOutput()->GetSpacing();
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];
  ImageType::SizeType  size = ImageReader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::RegionType region;
  ImageType::IndexType start = ImageReader->GetOutput()->GetRequestedRegion().GetIndex();
  region.SetSize( size );
  region.SetIndex( start );
  ImageType::IndexType pixelIndex;
  //Smoothed DT
  FilterType::Pointer filterX = FilterType::New();
  FilterType::Pointer filterY = FilterType::New();
  FilterType::Pointer filterZ = FilterType::New();
  filterX->SetDirection( 0 );   // 0 --> X direction
  filterY->SetDirection( 1 );   // 1 --> Y direction
  filterZ->SetDirection( 2 );   // 2 --> Z direction
  filterX->SetOrder( FilterType::ZeroOrder );
  filterY->SetOrder( FilterType::ZeroOrder );
  filterZ->SetOrder( FilterType::ZeroOrder );
  filterX->SetNormalizeAcrossScale( false );
  filterY->SetNormalizeAcrossScale( false ); 
  filterZ->SetNormalizeAcrossScale( false ); 
  filterX->SetInput( ImageReader->GetOutput() );
  filterY->SetInput( filterX->GetOutput() );
  filterZ->SetInput( filterY->GetOutput() );
  const double sigma = 1;
  filterX->SetSigma( sigma );
  filterY->SetSigma( sigma );
  filterZ->SetSigma( sigma );
  filterZ->Update();
  ImageWriterType::Pointer writerSmoothDT = ImageWriterType::New();
  writerSmoothDT->SetInput( filterZ->GetOutput() );
  writerSmoothDT->SetFileName("temp/SmoothDT.img");
  writerSmoothDT->Update();

  //Get distance transform 
  float *DT;
  DT = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  //DT(i,j,k) = ImageReader->GetOutput()->GetPixel(pixelIndex);
	  DT(i,j,k) = filterZ->GetOutput()->GetPixel(pixelIndex);
	  DT(i,j,k) = DT(i,j,k)*DT(i,j,k);
	}

  //***********************DT Analyze Part********************************* 
  //Apply gradient - use 3D Sobel operator
  float *GradientDTX;
  GradientDTX = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  float *GradientDTY;
  GradientDTY = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  float *GradientDTZ;
  GradientDTZ = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  float *vecOriGraZenith;
  vecOriGraZenith = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  float *vecOriGraAzimuth;
  vecOriGraAzimuth = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  float posX,negX,posY,negY,posZ,negZ; // positive and negative part of gradient computation
  float weight;
  float norm; //used for normalize the orientation vector
  //------------------------Sobel Operator----------------------------------
  //core
  for(i=1;i<pRow-1;i++)
    for(j=1;j<pCol-1;j++)
      for(k=1;k<pSli-1;k++) 
	{ 
	  posX = 0;
	  negX = 0;
	  posY = 0;
	  negY = 0;
	  posZ = 0;
	  negZ = 0;
	  for(a=-1;a<=1;a++)
	    for(b=-1;b<=1;b++)
	      {
		if((abs(a)+abs(b))==0) weight = 4.0;
		if((abs(a)+abs(b))==1) weight = 2.0;
		if((abs(a)+abs(b))==2) weight = 1.0;
		posX = posX + DT(i+1,j+a,k+b)*weight;
		negX = negX + DT(i-1,j+a,k+b)*weight;
		posY = posY + DT(i+a,j+1,k+b)*weight;
		negY = negY + DT(i+a,j-1,k+b)*weight;
		posZ = posZ + DT(i+a,j+b,k+1)*weight;
		negZ = negZ + DT(i+a,j+b,k-1)*weight;
	      }
	  GradientDTX(i,j,k) = (posX-negX)/32.0;
	  GradientDTY(i,j,k) = (posY-negY)/32.0;
	  GradientDTZ(i,j,k) = (posZ-negZ)/32.0;
	}
  //6 slabs 
  for(i=1;i<pRow-1;i++)
    for(j=1;j<pCol-1;j++)
      {
	posX = DT(i+1,j-1,0)+2.0*DT(i+1,j,0)+DT(i+1,j+1,0);
	negX = DT(i-1,j-1,0)+2.0*DT(i-1,j,0)+DT(i-1,j+1,0);
	GradientDTX(i,j,0) = (posX-negX)/8.0;
	posX = DT(i+1,j-1,pSli-1)+2.0*DT(i+1,j,pSli-1)+DT(i+1,j+1,pSli-1);
	negX = DT(i-1,j-1,pSli-1)+2.0*DT(i-1,j,pSli-1)+DT(i-1,j+1,pSli-1);
	GradientDTX(i,j,pSli-1) = (posX-negX)/8.0;

	posY = DT(i-1,j+1,0)+2.0*DT(i,j+1,0)+DT(i+1,j+1,0);
	negY = DT(i-1,j-1,0)+2.0*DT(i,j-1,0)+DT(i+1,j-1,0);
	GradientDTY(i,j,0) = (posY-negY)/8.0;
	posY = DT(i-1,j+1,pSli-1)+2.0*DT(i,j+1,pSli-1)+DT(i+1,j+1,pSli-1);
	negY = DT(i-1,j-1,pSli-1)+2.0*DT(i,j-1,pSli-1)+DT(i+1,j-1,pSli-1);
	GradientDTY(i,j,pSli-1) = (posY-negY)/8.0;

	GradientDTZ(i,j,0) = DT(i,j,1)-DT(i,j,0);
	GradientDTZ(i,j,pSli-1) = DT(i,j,pSli-1)-DT(i,j,pSli-2);
      }
  for(j=1;j<pCol-1;j++)
    for(k=1;k<pSli-1;k++) 
      {
	GradientDTX(0,j,k) = DT(1,j,k)-DT(0,j,k);
	GradientDTX(pRow-1,j,k) = DT(pRow-1,j,k)-DT(pRow-2,j,k);

	posY = DT(0,j+1,k-1)+2.0*DT(0,j+1,k)+DT(0,j+1,k+1);
	negY = DT(0,j-1,k-1)+2.0*DT(0,j-1,k)+DT(0,j-1,k+1);
	GradientDTY(0,j,k) = (posY-negY)/8.0;
	posY = DT(pRow-1,j+1,k-1)+2.0*DT(pRow-1,j+1,k)+DT(pRow-1,j+1,k+1);
	negY = DT(pRow-1,j-1,k-1)+2.0*DT(pRow-1,j-1,k)+DT(pRow-1,j-1,k+1);
	GradientDTY(pRow-1,j,k) = (posY-negY)/8.0;

	posZ = DT(0,j-1,k+1)+2.0*DT(0,j,k+1)+DT(0,j+1,k+1);
	negZ = DT(0,j-1,k-1)+2.0*DT(0,j,k-1)+DT(0,j+1,k-1);
	GradientDTZ(0,j,k) = (posZ-negZ)/8.0;
	posZ = DT(pRow-1,j-1,k+1)+2.0*DT(pRow-1,j,k+1)+DT(pRow-1,j+1,k+1);
	negZ = DT(pRow-1,j-1,k-1)+2.0*DT(pRow-1,j,k-1)+DT(pRow-1,j+1,k-1);
	GradientDTZ(pRow-1,j,k) = (posZ-negZ)/8.0;
      }
  for(i=1;i<pRow-1;i++)
    for(k=1;k<pSli-1;k++)
      {
	posX = DT(i+1,0,k-1)+2.0*DT(i+1,0,k)+DT(i+1,0,k+1);
	negX = DT(i-1,0,k-1)+2.0*DT(i-1,0,k)+DT(i-1,0,k+1);
	GradientDTX(i,0,k) = (posX-negX)/8.0;
	posX = DT(i+1,pCol-1,k-1)+2.0*DT(i+1,pCol-1,k)+DT(i+1,pCol-1,k+1);
	negX = DT(i-1,pCol-1,k-1)+2.0*DT(i-1,pCol-1,k)+DT(i-1,pCol-1,k+1);
	GradientDTX(i,pCol-1,k) = (posX-negX)/8.0;

	GradientDTY(i,0,k) = DT(i,1,k)-DT(i,0,k);
	GradientDTY(i,pCol-1,k) = DT(i,pCol-1,k)-DT(i,pCol-2,k);

	posZ = DT(i-1,0,k+1)+2.0*DT(i,0,k+1)+DT(i+1,0,k+1);
	negZ = DT(i-1,0,k-1)+2.0*DT(i,0,k-1)+DT(i+1,0,k-1);
	GradientDTZ(i,0,k) = (posZ-negZ)/8.0;
	posZ = DT(i-1,pCol-1,k+1)+2.0*DT(i,pCol-1,k+1)+DT(i+1,pCol-1,k+1);
	negZ = DT(i-1,pCol-1,k-1)+2.0*DT(i,pCol-1,k-1)+DT(i+1,pCol-1,k-1);
	GradientDTZ(i,pCol-1,k) = (posZ-negZ)/8.0;	
      }
  //12 lines
  for(i=1;i<pRow-1;i++)
    {
      GradientDTX(i,0,0) = (DT(i+1,0,0)-DT(i-1,0,0))/2.0;
      GradientDTY(i,0,0) = DT(i,1,0)-DT(i,0,0);
      GradientDTZ(i,0,0) = DT(i,0,1)-DT(i,0,0);
      GradientDTX(i,pCol-1,0) = (DT(i+1,pCol-1,0)-DT(i-1,pCol-1,0))/2.0;
      GradientDTY(i,pCol-1,0) = DT(i,pCol-1,0)-DT(i,pCol-2,0);
      GradientDTZ(i,pCol-1,0) = DT(i,pCol-1,1)-DT(i,pCol-1,0);
      GradientDTX(i,0,pSli-1) = (DT(i+1,0,pSli-1)-DT(i-1,0,pSli-1))/2.0;
      GradientDTY(i,0,pSli-1) = DT(i,1,pSli-1)-DT(i,0,pSli-1);
      GradientDTZ(i,0,pSli-1) = DT(i,0,pSli-1)-DT(i,0,pSli-2);
      GradientDTX(i,pCol-1,pSli-1) = (DT(i+1,pCol-1,pSli-1)-DT(i-1,pCol-1,pSli-1))/2.0;
      GradientDTY(i,pCol-1,pSli-1) = DT(i,pCol-1,pSli-1)-DT(i,pCol-2,pSli-1);
      GradientDTZ(i,pCol-1,pSli-1) = DT(i,pCol-1,pSli-1)-DT(i,pCol-1,pSli-2);     
    }
  for(j=1;j<pCol-1;j++)
    {
      GradientDTX(0,j,0) = DT(1,j,0)-DT(0,j,0);
      GradientDTY(0,j,0) = (DT(0,j+1,0)-DT(0,j-1,0))/2.0;
      GradientDTZ(0,j,0) = DT(0,j,1)-DT(0,j,0);
      GradientDTX(pRow-1,j,0) = DT(pRow-1,j,0)-DT(pRow-2,j,0);
      GradientDTY(pRow-1,j,0) = (DT(pRow-1,j+1,0)-DT(pRow-1,j-1,0))/2.0;
      GradientDTZ(pRow-1,j,0) = DT(pRow-1,j,1)-DT(pRow-1,j,0);
      GradientDTX(0,j,pSli-1) = DT(1,j,pSli-1)-DT(0,j,pSli-1);
      GradientDTY(0,j,pSli-1) = (DT(0,j+1,pSli-1)-DT(0,j-1,pSli-1))/2.0;
      GradientDTZ(0,j,pSli-1) = DT(0,j,pSli-1)-DT(0,j,pSli-2);
      GradientDTX(pRow-1,j,pSli-1) = DT(pRow-1,j,pSli-1)-DT(pRow-2,j,pSli-1);
      GradientDTY(pRow-1,j,pSli-1) = (DT(pRow-1,j+1,pSli-1)-DT(pRow-1,j-1,pSli-1))/2.0;
      GradientDTZ(pRow-1,j,pSli-1) = DT(pRow-1,j,pSli-1)-DT(pRow-1,j,pSli-2);      
    }
  for(k=1;k<pSli-1;k++) 
    {
      GradientDTX(0,0,k) = DT(1,0,k)-DT(0,0,k);
      GradientDTY(0,0,k) = DT(0,1,k)-DT(0,0,k);
      GradientDTZ(0,0,k) = (DT(0,0,k+1)-DT(0,0,k-1))/2.0;
      GradientDTX(pRow-1,0,k) = DT(pRow-1,0,k)-DT(pRow-2,0,k);
      GradientDTY(pRow-1,0,k) = DT(pRow-1,1,k)-DT(pRow-1,0,k);
      GradientDTZ(pRow-1,0,k) = (DT(pRow-1,0,k+1)-DT(pRow-1,0,k-1))/2.0;
      GradientDTX(0,pCol-1,k) = DT(1,pCol-1,k)-DT(0,pCol-1,k);
      GradientDTY(0,pCol-1,k) = DT(0,pCol-1,k)-DT(0,pCol-2,k);     
      GradientDTZ(0,pCol-1,k) = (DT(0,pCol-1,k+1)-DT(0,pCol-1,k-1))/2.0;
      GradientDTX(pRow-1,pCol-1,k) = DT(pRow-1,pCol-1,k)-DT(pRow-2,pCol-1,k);
      GradientDTY(pRow-1,pCol-1,k) = DT(pRow-1,pCol-1,k)-DT(pRow-1,pCol-2,k);      
      GradientDTZ(pRow-1,pCol-1,k) = (DT(pRow-1,pCol-1,k+1)-DT(pRow-1,pCol-1,k-1))/2.0;  
    }
  //8 points
  GradientDTX(0,0,0) = DT(1,0,0)-DT(0,0,0);
  GradientDTY(0,0,0) = DT(0,1,0)-DT(0,0,0);
  GradientDTZ(0,0,0) = DT(0,0,1)-DT(0,0,0);
  GradientDTX(0,0,pSli-1) = DT(1,0,pSli-1)-DT(0,0,pSli-1);
  GradientDTY(0,0,pSli-1) = DT(0,1,pSli-1)-DT(0,0,pSli-1);
  GradientDTZ(0,0,pSli-1) = DT(0,0,pSli-1)-DT(0,0,pSli-2);
  GradientDTX(0,pCol-1,0) = DT(1,pCol-1,0)-DT(0,pCol-1,0);
  GradientDTY(0,pCol-1,0) = DT(0,pCol-1,0)-DT(0,pCol-2,0);
  GradientDTZ(0,pCol-1,0) = DT(0,pCol-1,1)-DT(0,pCol-1,0);
  GradientDTX(pRow-1,0,0) = DT(pRow-1,0,0)-DT(pRow-2,0,0);
  GradientDTY(pRow-1,0,0) = DT(pRow-1,1,0)-DT(pRow-1,0,0);
  GradientDTZ(pRow-1,0,0) = DT(pRow-1,0,1)-DT(pRow-1,0,0);
  GradientDTX(0,pCol-1,pSli-1) = DT(1,pCol-1,pSli-1)-DT(0,pCol-1,pSli-1);
  GradientDTY(0,pCol-1,pSli-1) = DT(0,pCol-1,pSli-1)-DT(0,pCol-2,pSli-1);
  GradientDTZ(0,pCol-1,pSli-1) = DT(0,pCol-1,pSli-1)-DT(0,pCol-1,pSli-2);
  GradientDTX(pRow-1,0,pSli-1) = DT(pRow-1,0,pSli-1)-DT(pRow-2,0,pSli-1);
  GradientDTY(pRow-1,0,pSli-1) = DT(pRow-1,1,pSli-1)-DT(pRow-1,0,pSli-1);
  GradientDTZ(pRow-1,0,pSli-1) = DT(pRow-1,0,pSli-1)-DT(pRow-1,0,pSli-2);
  GradientDTX(pRow-1,pCol-1,0) = DT(pRow-1,pCol-1,0)-DT(pRow-2,pCol-1,0);
  GradientDTY(pRow-1,pCol-1,0) = DT(pRow-1,pCol-1,0)-DT(pRow-1,pCol-2,0);
  GradientDTZ(pRow-1,pCol-1,0) = DT(pRow-1,pCol-1,1)-DT(pRow-1,pCol-1,0);
  GradientDTX(pRow-1,pCol-1,pSli-1) = DT(pRow-1,pCol-1,pSli-1)-DT(pRow-2,pCol-1,pSli-1);
  GradientDTY(pRow-1,pCol-1,pSli-1) = DT(pRow-1,pCol-1,pSli-1)-DT(pRow-1,pCol-2,pSli-1);
  GradientDTZ(pRow-1,pCol-1,pSli-1) = DT(pRow-1,pCol-1,pSli-1)-DT(pRow-1,pCol-1,pSli-2);
 

  //------------------Write Sobel result---------------------------------------
  ImageType::Pointer GradientDTXImage = ImageType::New();
  GradientDTXImage->SetRegions( region );
  GradientDTXImage->SetSpacing( spacing );
  GradientDTXImage->Allocate();
  ImageType::Pointer GradientDTYImage = ImageType::New();
  GradientDTYImage->SetRegions( region );
  GradientDTYImage->SetSpacing( spacing );
  GradientDTYImage->Allocate();
  ImageType::Pointer GradientDTZImage = ImageType::New();
  GradientDTZImage->SetRegions( region );
  GradientDTZImage->SetSpacing( spacing );
  GradientDTZImage->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k; 
	  if(absl(GradientDTX(i,j,k))<0.00001) GradientDTX(i,j,k)=0;
 	  if(absl(GradientDTY(i,j,k))<0.00001) GradientDTY(i,j,k)=0;
	  if(absl(GradientDTZ(i,j,k))<0.00001) GradientDTZ(i,j,k)=0;
	  norm = GradientDTX(i,j,k)*GradientDTX(i,j,k);
	  norm = norm + GradientDTY(i,j,k)*GradientDTY(i,j,k);
	  norm = norm + GradientDTZ(i,j,k)*GradientDTZ(i,j,k);
	  norm = sqrt(norm);
	  if(norm == 0) norm = 1;
	  GradientDTXImage->SetPixel(pixelIndex,GradientDTX(i,j,k)/norm);
	  GradientDTYImage->SetPixel(pixelIndex,GradientDTY(i,j,k)/norm);
	  GradientDTZImage->SetPixel(pixelIndex,GradientDTZ(i,j,k)/norm);
	}
  ImageWriterType::Pointer writerGradientDTX = ImageWriterType::New();
  writerGradientDTX->SetInput( GradientDTXImage );
  writerGradientDTX->SetFileName("temp/GradientDTX.img");
  writerGradientDTX->Update();
  ImageWriterType::Pointer writerGradientDTY = ImageWriterType::New();
  writerGradientDTY->SetInput( GradientDTYImage );
  writerGradientDTY->SetFileName("temp/GradientDTY.img");
  writerGradientDTY->Update();
  ImageWriterType::Pointer writerGradientDTZ = ImageWriterType::New();
  writerGradientDTZ->SetInput( GradientDTZImage );
  writerGradientDTZ->SetFileName("temp/GradientDTZ.img");
  writerGradientDTZ->Update();
  
  //----------------------Smooth the gradient field ----------------------------------------------
  //Gradient magnitude weighted (center line attenuation)
  int x,y,z;
  int smoothSize=atoi(argv[1]); //neighborhood size: 2n+1 * 2n+1 * 2n+1
  float *SmoothGraDTX;
  SmoothGraDTX = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  float *SmoothGraDTY;
  SmoothGraDTY = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  float *SmoothGraDTZ;
  SmoothGraDTZ = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  //First,  compute gradient magnitude
  //Second, record gradient pairs weighted by gradient magnitude in matrix
  //Third,  add additional 0s of the neighborhood size to the matrix
  //Fourth, compute SVD for orientation result
  vnl_matrix<float> graPairs((smoothSize*2+1)*(smoothSize*2+1)*(smoothSize*2+1)*2,3,0.0);  
  float tempLength;
  for(i=0;i<pRow;i++)
    {
      for(j=0;j<pCol;j++)
	for(k=0;k<pSli;k++)
	  {
	    for(t=0;t<(smoothSize*2+1)*(smoothSize*2+1)*(smoothSize*2+1)*2;t++)
	      {
		graPairs(t,0)=0.0;
		graPairs(t,1)=0.0;
		graPairs(t,2)=0.0;
	      }
	    t=0;
	    //compute gradient magnitude and record in pairs
	    for(x=-smoothSize;x<=smoothSize;x++)
	      for(y=-smoothSize;y<=smoothSize;y++)	  
		for(z=-smoothSize;z<=smoothSize;z++)	  
		  if(((i+x)>=0)&&((i+x)<pRow)&&((j+y)>=0)&&((j+y)<pCol)&&((k+z)>=0)&&((k+z)<pSli))
		    {
		      tempLength = GradientDTX(i+x,j+y,k+z)*GradientDTX(i+x,j+y,k+z);
		      tempLength = tempLength + GradientDTY(i+x,j+y,k+z)*GradientDTY(i+x,j+y,k+z);
		      tempLength = tempLength + GradientDTZ(i+x,j+y,k+z)*GradientDTZ(i+x,j+y,k+z);
		      graPairs(t,0) = GradientDTX(i+x,j+y,k+z)*tempLength;
		      graPairs(t,1) = GradientDTY(i+x,j+y,k+z)*tempLength;
		      graPairs(t,2) = GradientDTZ(i+x,j+y,k+z)*tempLength;
		      t = t+1;
		    }
	    //SVD
	    vnl_matrix<float> SVDV = vnl_svd<float> (graPairs).V(); 
	    SmoothGraDTX(i,j,k) = SVDV(0,0);
	    SmoothGraDTY(i,j,k) = SVDV(1,0);
	    SmoothGraDTZ(i,j,k) = SVDV(2,0);
	    if(signf(SmoothGraDTX(i,j,k))!=signf(GradientDTX(i,j,k))) 
	      SmoothGraDTX(i,j,k) = - SmoothGraDTX(i,j,k);
	    if(signf(SmoothGraDTY(i,j,k))!=signf(GradientDTY(i,j,k))) 
	      SmoothGraDTY(i,j,k) = - SmoothGraDTY(i,j,k);
	    if(signf(SmoothGraDTZ(i,j,k))!=signf(GradientDTZ(i,j,k))) 
	      SmoothGraDTZ(i,j,k) = - SmoothGraDTZ(i,j,k);
	  }
      std::cout<<i<<std::endl;
    }
  
  //------------------Write Smoothed result---------------------------------------
  ImageType::Pointer SmoothGraDTXImage = ImageType::New();
  SmoothGraDTXImage->SetRegions( region );
  SmoothGraDTXImage->SetSpacing( spacing );
  SmoothGraDTXImage->Allocate();
  ImageType::Pointer SmoothGraDTYImage = ImageType::New();
  SmoothGraDTYImage->SetRegions( region );
  SmoothGraDTYImage->SetSpacing( spacing );
  SmoothGraDTYImage->Allocate();
  ImageType::Pointer SmoothGraDTZImage = ImageType::New();
  SmoothGraDTZImage->SetRegions( region );
  SmoothGraDTZImage->SetSpacing( spacing );
  SmoothGraDTZImage->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k; 
	  if(absl(SmoothGraDTX(i,j,k))<0.00001) SmoothGraDTX(i,j,k)=0;
 	  if(absl(SmoothGraDTY(i,j,k))<0.00001) SmoothGraDTY(i,j,k)=0;
 	  if(absl(SmoothGraDTZ(i,j,k))<0.00001) SmoothGraDTZ(i,j,k)=0;
 	  if(absl(absl(SmoothGraDTX(i,j,k)-1))<0.00001) 
	    SmoothGraDTX(i,j,k)=1.0*signf(SmoothGraDTX(i,j,k));
 	  if(absl(absl(SmoothGraDTY(i,j,k)-1))<0.00001) 
	    SmoothGraDTY(i,j,k)=1.0*signf(SmoothGraDTY(i,j,k));
 	  if(absl(absl(SmoothGraDTZ(i,j,k)-1))<0.00001) 
	    SmoothGraDTZ(i,j,k)=1.0*signf(SmoothGraDTZ(i,j,k));

	  norm = SmoothGraDTX(i,j,k)*SmoothGraDTX(i,j,k);
	  norm = norm + SmoothGraDTY(i,j,k)*SmoothGraDTY(i,j,k);
	  norm = norm + SmoothGraDTZ(i,j,k)*SmoothGraDTZ(i,j,k);
	  norm = sqrt(norm);
	  if(norm == 0) norm = 1;

	  SmoothGraDTXImage->SetPixel(pixelIndex,SmoothGraDTX(i,j,k)/norm);
	  SmoothGraDTYImage->SetPixel(pixelIndex,SmoothGraDTY(i,j,k)/norm);
	  SmoothGraDTZImage->SetPixel(pixelIndex,SmoothGraDTZ(i,j,k)/norm);
	}
  ImageWriterType::Pointer writerSmoothGraDTX = ImageWriterType::New();
  writerSmoothGraDTX->SetInput( SmoothGraDTXImage );
  writerSmoothGraDTX->SetFileName("temp/GradientDTXSmooth.img");
  writerSmoothGraDTX->Update();
  ImageWriterType::Pointer writerSmoothGraDTY = ImageWriterType::New();
  writerSmoothGraDTY->SetInput( SmoothGraDTYImage );
  writerSmoothGraDTY->SetFileName("temp/GradientDTYSmooth.img");
  writerSmoothGraDTY->Update();
  ImageWriterType::Pointer writerSmoothGraDTZ = ImageWriterType::New();
  writerSmoothGraDTZ->SetInput( SmoothGraDTZImage );
  writerSmoothGraDTZ->SetFileName("temp/GradientDTZSmooth.img");
  writerSmoothGraDTZ->Update();

  free(DT);
  free(GradientDTX);
  free(GradientDTY);
  free(GradientDTZ);
  free(SmoothGraDTX);
  free(SmoothGraDTY);
  free(SmoothGraDTZ);
  
  return EXIT_SUCCESS;
}
