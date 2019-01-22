#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "ShapeTensor3D.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters
 
  if( argc < 2 ) 
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " smoothFlag" << std::endl;
      return EXIT_FAILURE;
    }

  //Read image
  ImageReaderType::Pointer DTImageReader = ImageReaderType::New();
  DTImageReader->SetFileName( "temp/DT.img" );
  DTImageReader->Update();
  ImageReaderType::Pointer GraXImageReader = ImageReaderType::New();
  ImageReaderType::Pointer GraYImageReader = ImageReaderType::New();
  ImageReaderType::Pointer GraZImageReader = ImageReaderType::New();
  GraXImageReader->SetFileName( "temp/GradientDTX.img" );
  GraYImageReader->SetFileName( "temp/GradientDTY.img" );
  GraZImageReader->SetFileName( "temp/GradientDTZ.img" );
  GraXImageReader->Update();
  GraYImageReader->Update();
  GraZImageReader->Update();
  ImageType::SpacingType spacing = DTImageReader->GetOutput()->GetSpacing();
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];
  ImageType::SizeType  size = DTImageReader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::RegionType region;
  ImageType::IndexType start = DTImageReader->GetOutput()->GetRequestedRegion().GetIndex();
  region.SetSize( size );
  region.SetIndex( start );
  ImageType::IndexType pixelIndex;

  float *DT;
  DT = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  VOXEL *vecOri;
  vecOri = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);
  VOXEL *tangentOri1;
  tangentOri1 = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);
  VOXEL *tangentOri2;
  tangentOri2 = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  DT(i,j,k) = DTImageReader->GetOutput()->GetPixel(pixelIndex);
	  DT(i,j,k) = DT(i,j,k)*DT(i,j,k);
	  vecOri(i,j,k).x = GraXImageReader->GetOutput()->GetPixel(pixelIndex);
	  vecOri(i,j,k).y = GraYImageReader->GetOutput()->GetPixel(pixelIndex);
	  vecOri(i,j,k).z = GraZImageReader->GetOutput()->GetPixel(pixelIndex);
	}

  //------------------------Trace in orthogonal plane----------------------------------------------- 
  //1. Compute 2 unit tracing orientation vector in tangent plane
  //2. Trace structure size along the 2 orientation
  float norm;
  //Get tangent vector
  bool smoothFlag = atoi(argv[1]);
  FindTangentOri(DT, vecOri, tangentOri1, tangentOri2, smoothFlag);
  
  ImageType::Pointer tangOri1X = ImageType::New();
  tangOri1X->SetRegions( region );
  tangOri1X->Allocate();
  ImageType::Pointer tangOri1Y = ImageType::New();
  tangOri1Y->SetRegions( region );
  tangOri1Y->Allocate();
  ImageType::Pointer tangOri1Z = ImageType::New();
  tangOri1Z->SetRegions( region );
  tangOri1Z->Allocate();
  ImageType::Pointer tangOri2X = ImageType::New();
  tangOri2X->SetRegions( region );
  tangOri2X->Allocate();
  ImageType::Pointer tangOri2Y = ImageType::New();
  tangOri2Y->SetRegions( region );
  tangOri2Y->Allocate();
  ImageType::Pointer tangOri2Z = ImageType::New();
  tangOri2Z->SetRegions( region );
  tangOri2Z->Allocate();


  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;	  
	  norm = tangentOri1(i,j,k).x*tangentOri1(i,j,k).x;
	  norm = norm + tangentOri1(i,j,k).y*tangentOri1(i,j,k).y;
	  norm = norm + tangentOri1(i,j,k).z*tangentOri1(i,j,k).z;
	  norm = sqrt(norm);
	  if(norm == 0) norm = 1;
	  tangentOri1(i,j,k).x = tangentOri1(i,j,k).x/norm;
	  tangentOri1(i,j,k).y = tangentOri1(i,j,k).y/norm;
	  tangentOri1(i,j,k).z = tangentOri1(i,j,k).z/norm;
  	  tangOri1X->SetPixel(pixelIndex,tangentOri1(i,j,k).x);
	  tangOri1Y->SetPixel(pixelIndex,tangentOri1(i,j,k).y);
	  tangOri1Z->SetPixel(pixelIndex,tangentOri1(i,j,k).z);

	  norm = tangentOri2(i,j,k).x*tangentOri2(i,j,k).x;
	  norm = norm + tangentOri2(i,j,k).y*tangentOri2(i,j,k).y;
	  norm = norm + tangentOri2(i,j,k).z*tangentOri2(i,j,k).z;
	  norm = sqrt(norm);
	  if(norm == 0) norm = 1;
	  tangentOri2(i,j,k).x = tangentOri2(i,j,k).x/norm;
	  tangentOri2(i,j,k).y = tangentOri2(i,j,k).y/norm;
	  tangentOri2(i,j,k).z = tangentOri2(i,j,k).z/norm;	  
	  tangOri2X->SetPixel(pixelIndex,tangentOri2(i,j,k).x);
	  tangOri2Y->SetPixel(pixelIndex,tangentOri2(i,j,k).y);
	  tangOri2Z->SetPixel(pixelIndex,tangentOri2(i,j,k).z);
	}
  
  ImageWriterType::Pointer writertangOri1X = ImageWriterType::New();
  writertangOri1X->SetInput( tangOri1X );
  writertangOri1X->SetFileName("temp/tangOri1X.img");
  writertangOri1X->Update(); 
  ImageWriterType::Pointer writertangOri1Y = ImageWriterType::New();
  writertangOri1Y->SetInput( tangOri1Y );
  writertangOri1Y->SetFileName("temp/tangOri1Y.img");
  writertangOri1Y->Update();
  ImageWriterType::Pointer writertangOri1Z = ImageWriterType::New();
  writertangOri1Z->SetInput( tangOri1Z );
  writertangOri1Z->SetFileName("temp/tangOri1Z.img");
  writertangOri1Z->Update();
  ImageWriterType::Pointer writertangOri2X = ImageWriterType::New();
  writertangOri2X->SetInput( tangOri2X );
  writertangOri2X->SetFileName("temp/tangOri2X.img");
  writertangOri2X->Update(); 
  ImageWriterType::Pointer writertangOri2Y = ImageWriterType::New();
  writertangOri2Y->SetInput( tangOri2Y );
  writertangOri2Y->SetFileName("temp/tangOri2Y.img");
  writertangOri2Y->Update();
  ImageWriterType::Pointer writertangOri2Z = ImageWriterType::New();
  writertangOri2Z->SetInput(tangOri2Z  );
  writertangOri2Z->SetFileName("temp/tangOri2Z.img");
  writertangOri2Z->Update(); 
  
  
  free(DT);
  free(tangentOri1);
  free(tangentOri2);
  
  return EXIT_SUCCESS;
}
