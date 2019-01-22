#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#define PI 3.1415926
#include "Wkt.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage outputImage iterationTime" << std::endl;
      return EXIT_FAILURE;
    }
  
  //Read image
  ImageReaderType::Pointer ImageReader = ImageReaderType::New();
  ImageReader->SetFileName( argv[1] );
  ImageReader->Update();
  ImageType::SpacingType spacing = ImageReader->GetOutput()->GetSpacing();  
  ImageType::SizeType  size = ImageReader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];

  ImageType::RegionType region;
  ImageType::IndexType start = {0};
  region.SetSize( size );
  region.SetIndex( start );

  int iterationTime = atoi(argv[3]);
  ImageType::IndexType pixelIndex;
  float *oriImage;
  oriImage =(float *)malloc(sizeof(float *) * pRow * pCol * pSli);  

  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)      
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  oriImage(i,j,k) = ImageReader->GetOutput()->GetPixel(pixelIndex);
	}

  //Filtering Parameters
  float sigma = 0.5;
  float rho = 4;
  float alpha = 0.5;
  float Ci = 1;
  float tau = 0.1;
  
  int iter;
  std::cout<<"0%"<<std::flush;
  for(iter=0;iter<iterationTime;iter++)
    {
      //Define filters
      ImageType::Pointer imageOri = ImageType::New();
      imageOri->SetRegions( region );
      imageOri->Allocate();
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  for(k=0;k<pSli;k++)
	    {
	      pixelIndex[0]=i;
	      pixelIndex[1]=j;
	      pixelIndex[2]=k;	
	      imageOri->SetPixel(pixelIndex,oriImage(i,j,k));
	    }
      
      //Regularize image by Gaussian filtering with sigma 
      GaussFilterType::Pointer RegFilterX = GaussFilterType::New();
      RegFilterX->SetInput( imageOri );
      RegFilterX->SetSigma( sigma );
      RegFilterX->SetDirection( 0 );
      RegFilterX->SetOrder( GaussFilterType::ZeroOrder );
      RegFilterX->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer RegFilterY = GaussFilterType::New();
      RegFilterY->SetInput( RegFilterX->GetOutput() );
      RegFilterY->SetSigma( sigma );
      RegFilterY->SetDirection( 1 );
      RegFilterY->SetOrder( GaussFilterType::ZeroOrder );
      RegFilterY->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer RegFilterZ = GaussFilterType::New();
      RegFilterZ->SetInput( RegFilterY->GetOutput() );
      RegFilterZ->SetSigma( sigma );
      RegFilterZ->SetDirection( 2 );
      RegFilterZ->SetOrder( GaussFilterType::ZeroOrder );
      RegFilterZ->SetNormalizeAcrossScale( false );

      //Gradient filter
      DevFilterType::Pointer DevFilterX = DevFilterType::New();
      DevFilterX->SetInput(RegFilterZ->GetOutput());
      DevFilterX->SetOrder(1);
      DevFilterX->SetDirection(0);
      DevFilterX->Update();
      DevFilterType::Pointer DevFilterY = DevFilterType::New();
      DevFilterY->SetInput(RegFilterZ->GetOutput());
      DevFilterY->SetOrder(1);
      DevFilterY->SetDirection(1);
      DevFilterY->Update();
      DevFilterType::Pointer DevFilterZ = DevFilterType::New();
      DevFilterZ->SetInput(RegFilterZ->GetOutput());
      DevFilterZ->SetOrder(1);
      DevFilterZ->SetDirection(2);
      DevFilterZ->Update();

      //Structure Tensor Images
      ImageType::Pointer XXImage = ImageType::New();
      XXImage->SetRegions( region );
      XXImage->Allocate();
      ImageType::Pointer XYImage = ImageType::New();
      XYImage->SetRegions( region );
      XYImage->Allocate();
      ImageType::Pointer XZImage = ImageType::New();
      XZImage->SetRegions( region );
      XZImage->Allocate();
      ImageType::Pointer YYImage = ImageType::New();
      YYImage->SetRegions( region );
      YYImage->Allocate();
      ImageType::Pointer YZImage = ImageType::New();
      YZImage->SetRegions( region );
      YZImage->Allocate();
      ImageType::Pointer ZZImage = ImageType::New();
      ZZImage->SetRegions( region );
      ZZImage->Allocate();
      float tempX, tempY, tempZ;
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  for(k=0;k<pSli;k++)	  
	    {
	      pixelIndex[0]=i;
	      pixelIndex[1]=j;
	      pixelIndex[2]=k;	      
	      tempX = DevFilterX->GetOutput()->GetPixel(pixelIndex);
	      tempY = DevFilterY->GetOutput()->GetPixel(pixelIndex);
	      tempZ = DevFilterZ->GetOutput()->GetPixel(pixelIndex);
	      XXImage->SetPixel(pixelIndex,tempX*tempX);
	      XYImage->SetPixel(pixelIndex,tempX*tempY);
	      XZImage->SetPixel(pixelIndex,tempX*tempZ);
	      YYImage->SetPixel(pixelIndex,tempY*tempY);
	      YZImage->SetPixel(pixelIndex,tempY*tempZ);
	      ZZImage->SetPixel(pixelIndex,tempZ*tempZ);
	    }
      
      //Componentwise Gaussian average
      //XX
      GaussFilterType::Pointer AvgFilterXXX = GaussFilterType::New();
      AvgFilterXXX->SetInput( XXImage );
      AvgFilterXXX->SetSigma( rho );
      AvgFilterXXX->SetDirection( 0 );
      AvgFilterXXX->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXXX->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterXXY = GaussFilterType::New();
      AvgFilterXXY->SetInput( AvgFilterXXX->GetOutput() );
      AvgFilterXXY->SetSigma( rho );
      AvgFilterXXY->SetDirection( 1 );
      AvgFilterXXY->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXXY->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterXXZ = GaussFilterType::New();
      AvgFilterXXZ->SetInput( AvgFilterXXY->GetOutput() );
      AvgFilterXXZ->SetSigma( rho );
      AvgFilterXXZ->SetDirection( 2 );
      AvgFilterXXZ->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXXZ->SetNormalizeAcrossScale( false );
      AvgFilterXXZ->Update();
      //XY
      GaussFilterType::Pointer AvgFilterXYX = GaussFilterType::New();
      AvgFilterXYX->SetInput( XYImage );
      AvgFilterXYX->SetSigma( rho );
      AvgFilterXYX->SetDirection( 0 );
      AvgFilterXYX->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXYX->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterXYY = GaussFilterType::New();
      AvgFilterXYY->SetInput( AvgFilterXYX->GetOutput() );
      AvgFilterXYY->SetSigma( rho );
      AvgFilterXYY->SetDirection( 1 );
      AvgFilterXYY->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXYY->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterXYZ = GaussFilterType::New();
      AvgFilterXYZ->SetInput( AvgFilterXYY->GetOutput() );
      AvgFilterXYZ->SetSigma( rho );
      AvgFilterXYZ->SetDirection( 2 );
      AvgFilterXYZ->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXYZ->SetNormalizeAcrossScale( false );
      AvgFilterXYZ->Update();
      //XZ
      GaussFilterType::Pointer AvgFilterXZX = GaussFilterType::New();
      AvgFilterXZX->SetInput( XZImage );
      AvgFilterXZX->SetSigma( rho );
      AvgFilterXZX->SetDirection( 0 );
      AvgFilterXZX->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXZX->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterXZY = GaussFilterType::New();
      AvgFilterXZY->SetInput( AvgFilterXZX->GetOutput() );
      AvgFilterXZY->SetSigma( rho );
      AvgFilterXZY->SetDirection( 1 );
      AvgFilterXZY->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXZY->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterXZZ = GaussFilterType::New();
      AvgFilterXZZ->SetInput( AvgFilterXZY->GetOutput() );
      AvgFilterXZZ->SetSigma( rho );
      AvgFilterXZZ->SetDirection( 2 );
      AvgFilterXZZ->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterXZZ->SetNormalizeAcrossScale( false );
      AvgFilterXZZ->Update();
      //YY
      GaussFilterType::Pointer AvgFilterYYX = GaussFilterType::New();
      AvgFilterYYX->SetInput( YYImage );
      AvgFilterYYX->SetSigma( rho );
      AvgFilterYYX->SetDirection( 0 );
      AvgFilterYYX->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterYYX->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterYYY = GaussFilterType::New();
      AvgFilterYYY->SetInput( AvgFilterYYX->GetOutput() );
      AvgFilterYYY->SetSigma( rho );
      AvgFilterYYY->SetDirection( 1 );
      AvgFilterYYY->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterYYY->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterYYZ = GaussFilterType::New();
      AvgFilterYYZ->SetInput( AvgFilterYYY->GetOutput() );
      AvgFilterYYZ->SetSigma( rho );
      AvgFilterYYZ->SetDirection( 2 );
      AvgFilterYYZ->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterYYZ->SetNormalizeAcrossScale( false );
      AvgFilterYYZ->Update();
      //YZ
      GaussFilterType::Pointer AvgFilterYZX = GaussFilterType::New();
      AvgFilterYZX->SetInput( YZImage );
      AvgFilterYZX->SetSigma( rho );
      AvgFilterYZX->SetDirection( 0 );
      AvgFilterYZX->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterYZX->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterYZY = GaussFilterType::New();
      AvgFilterYZY->SetInput( AvgFilterYZX->GetOutput() );
      AvgFilterYZY->SetSigma( rho );
      AvgFilterYZY->SetDirection( 1 );
      AvgFilterYZY->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterYZY->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterYZZ = GaussFilterType::New();
      AvgFilterYZZ->SetInput( AvgFilterYZY->GetOutput() );
      AvgFilterYZZ->SetSigma( rho );
      AvgFilterYZZ->SetDirection( 2 );
      AvgFilterYZZ->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterYZZ->SetNormalizeAcrossScale( false );
      AvgFilterYZZ->Update();
      //ZZ
      GaussFilterType::Pointer AvgFilterZZX = GaussFilterType::New();
      AvgFilterZZX->SetInput( ZZImage );
      AvgFilterZZX->SetSigma( rho );
      AvgFilterZZX->SetDirection( 0 );
      AvgFilterZZX->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterZZX->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterZZY = GaussFilterType::New();
      AvgFilterZZY->SetInput( AvgFilterZZX->GetOutput() );
      AvgFilterZZY->SetSigma( rho );
      AvgFilterZZY->SetDirection( 1 );
      AvgFilterZZY->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterZZY->SetNormalizeAcrossScale( false );
      GaussFilterType::Pointer AvgFilterZZZ = GaussFilterType::New();
      AvgFilterZZZ->SetInput( AvgFilterZZY->GetOutput() );
      AvgFilterZZZ->SetSigma( rho );
      AvgFilterZZZ->SetDirection( 1 );
      AvgFilterZZZ->SetOrder( GaussFilterType::ZeroOrder );
      AvgFilterZZZ->SetNormalizeAcrossScale( false );
      AvgFilterZZZ->Update();

      //Diffusion Tensor and Update image
      ImageType::Pointer diffXImage = ImageType::New();
      diffXImage->SetRegions( region );
      diffXImage->Allocate();
      ImageType::Pointer diffYImage = ImageType::New();
      diffYImage->SetRegions( region );
      diffYImage->Allocate();
      ImageType::Pointer diffZImage = ImageType::New();
      diffZImage->SetRegions( region );
      diffZImage->Allocate();
      ImageType::Pointer updateXImage = ImageType::New();
      updateXImage->SetRegions( region );
      updateXImage->Allocate();
      ImageType::Pointer updateYImage = ImageType::New();
      updateYImage->SetRegions( region );
      updateYImage->Allocate();
      ImageType::Pointer updateZImage = ImageType::New();
      updateZImage->SetRegions( region );
      updateZImage->Allocate();

      //Gradient of original image filter
      DevFilterType::Pointer DevOriFilterX = DevFilterType::New();
      DevOriFilterX->SetInput(imageOri);
      DevOriFilterX->SetOrder(1);
      DevOriFilterX->SetDirection(0);
      DevOriFilterX->Update();
      DevFilterType::Pointer DevOriFilterY = DevFilterType::New();
      DevOriFilterY->SetInput(imageOri);
      DevOriFilterY->SetOrder(1);
      DevOriFilterY->SetDirection(1);
      DevOriFilterY->Update();
      DevFilterType::Pointer DevOriFilterZ = DevFilterType::New();
      DevOriFilterZ->SetInput(imageOri);
      DevOriFilterZ->SetOrder(1);
      DevOriFilterZ->SetDirection(2);
      DevOriFilterZ->Update();

      //Weikert's Diffusion
      float j11,j12,j13,j22,j23,j33;
      float PCACheck; //check for singular
      float a1,a2,a3,b1,b2,b3,c1,c2,c3;
      float eig1, eig2, eig3;
      float lambda1,lambda2,lambda3;
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  for(k=0;k<pSli;k++)
	    {
	      pixelIndex[0]=i;
	      pixelIndex[1]=j;
	      pixelIndex[2]=k;
	      
	      j11 = AvgFilterXXZ->GetOutput()->GetPixel(pixelIndex);
	      j12 = AvgFilterXYZ->GetOutput()->GetPixel(pixelIndex);
	      j13 = AvgFilterXZZ->GetOutput()->GetPixel(pixelIndex);
	      j22 = AvgFilterYYZ->GetOutput()->GetPixel(pixelIndex);
	      j23 = AvgFilterYZZ->GetOutput()->GetPixel(pixelIndex);
	      j33 = AvgFilterZZZ->GetOutput()->GetPixel(pixelIndex);
	    
	      if((j11!=j11)||(j12!=j12)||(j13!=j13)||(j22!=j22)||(j23!=j23)||(j33!=j33))
		{
		  std::cout<<j11<<" "<<j12<<" "<<j13<<" "<<j22<<" "<<j23<<" "<<j33<<std::endl;
		  break;
		}
	      PCACheck = j11*(j22*j33+j23*j23)-j12*(j12*j33-j13*j23)+j13*(j12*j23-j13*j22);
	      if(fabs(PCACheck)<0.001)
		{
		  a1 = 1;
		  a2 = 0;
		  a3 = 0;
		  b1 = 0;
		  b2 = 1;
		  b3 = 0;
		  c1 = 0;
		  c2 = 0;
		  c3 = 1;
		  eig1 = 1;
		  eig2 = 1;
		  eig3 = 1;
		}
	      else
		{
		  //----------------------------------------------------------------
		  //3x3 eigensystem analysis of the structure tensor matrix
		  vnl_matrix<float> samplePCA(3,3,0.0);
		  samplePCA(0,0) = j11;
		  samplePCA(0,1) = j12;
		  samplePCA(0,2) = j13;
		  samplePCA(1,0) = j12;
		  samplePCA(1,1) = j22;
		  samplePCA(1,2) = j23;
		  samplePCA(2,0) = j13;
		  samplePCA(2,1) = j23;
		  samplePCA(2,2) = j33;
	      
		  vnl_matrix<float> SVDV = vnl_svd<float> (samplePCA).V(); 
		  vnl_matrix<float> SVDW = vnl_svd<float> (samplePCA).W();
		  a1 = SVDV(0,0);
		  a2 = SVDV(1,0);
		  a3 = SVDV(2,0);
		  b1 = SVDV(0,1);
		  b2 = SVDV(1,1);
		  b3 = SVDV(2,1);
		  c1 = SVDV(0,2);
		  c2 = SVDV(1,2);
		  c3 = SVDV(2,2);
		  eig1 = SVDW(0,0);
		  eig2 = SVDW(1,1);
		  eig3 = SVDW(2,2);
		}
	      
	      //------------------------------------------------------------------
	      //Convert to diffusion matrix
	      lambda1 = alpha;
	      lambda2 = alpha;
	      if((eig2==0)||(eig3==0))
		lambda3 = 1;
	      else
		lambda3 = alpha+(1-alpha)*exp(-log(2.0)*Ci/pow(eig2/(alpha+eig3),4));
	      float D11,D12,D13,D22,D23,D33;
	      D11 = lambda1*a1*a1+lambda2*b1*b1+lambda3*c1*c1;
	      D12 = lambda1*a1*a2+lambda2*b1*b2+lambda3*c1*c2;	      
	      D13 = lambda1*a1*a3+lambda2*b1*b3+lambda3*c1*c3;
	      D22 = lambda1*a2*a2+lambda2*b2*b2+lambda3*c2*c2;
	      D23 = lambda1*a2*a3+lambda2*b2*b3+lambda3*c2*c3;
	      D33 = lambda1*a3*a3+lambda2*b3*b3+lambda3*c3*c3;

	      //Compute update image
	      tempX = DevOriFilterX->GetOutput()->GetPixel(pixelIndex);
	      tempY = DevOriFilterY->GetOutput()->GetPixel(pixelIndex);
	      tempZ = DevOriFilterZ->GetOutput()->GetPixel(pixelIndex);	      
	      updateXImage->SetPixel(pixelIndex,D11*tempX+D12*tempY+D13*tempZ);
	      updateYImage->SetPixel(pixelIndex,D12*tempX+D22*tempY+D23*tempZ);
	      updateZImage->SetPixel(pixelIndex,D13*tempX+D23*tempY+D33*tempZ);	      
	    }
      
      //Gradient of update image filter
      DevFilterType::Pointer DevUpdFilterX = DevFilterType::New();
      DevUpdFilterX->SetInput(updateXImage);
      DevUpdFilterX->SetOrder(1);
      DevUpdFilterX->SetDirection(0);
      DevUpdFilterX->Update();
      DevFilterType::Pointer DevUpdFilterY = DevFilterType::New();
      DevUpdFilterY->SetInput(updateYImage);
      DevUpdFilterY->SetOrder(1);
      DevUpdFilterY->SetDirection(1);
      DevUpdFilterY->Update();
      DevFilterType::Pointer DevUpdFilterZ = DevFilterType::New();
      DevUpdFilterZ->SetInput(updateZImage);
      DevUpdFilterZ->SetOrder(1);
      DevUpdFilterZ->SetDirection(2);
      DevUpdFilterZ->Update();
      
      //Update
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  for(k=0;k<pSli;k++)
	    {
	      pixelIndex[0]=i;
	      pixelIndex[1]=j;
	      pixelIndex[2]=k;
	      tempX = DevUpdFilterX->GetOutput()->GetPixel(pixelIndex);
	      tempY = DevUpdFilterY->GetOutput()->GetPixel(pixelIndex);
	      tempZ = DevUpdFilterZ->GetOutput()->GetPixel(pixelIndex);
	      //if((i==12)&&(j==35)&&(k==62))
	      //std::cout<< oriImage(i,j,k)<<" ";
	      oriImage(i,j,k) = oriImage(i,j,k)+tau*(tempX+tempY+tempZ);
	      //if((i==12)&&(j==35)&&(k==62))
	      //std::cout<< oriImage(i,j,k)<<std::endl;
	    }
      imageOri->DisconnectPipeline(); 
      std::cout<<"\r";
      std::cout<<int((iter+1)*100/iterationTime)<<"%"<<std::flush;
    }
  std::cout<<std::endl;

  ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions( region );
  outImage->SetSpacing( spacing);
  outImage->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  outImage->SetPixel(pixelIndex,oriImage(i,j,k));
	} 
  
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( outImage );
  writer->SetFileName( argv[2] );
  writer->Update();
}
