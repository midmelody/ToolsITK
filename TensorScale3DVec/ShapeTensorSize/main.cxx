#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "ShapeTensor3D.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters
 
  if( argc < 1 ) 
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] <<  std::endl;
      return EXIT_FAILURE;
    }

  //Read image
  ImageReaderType::Pointer DTImageReader = ImageReaderType::New();
  DTImageReader->SetFileName( "temp/DT.img" );
  DTImageReader->Update();
  ImageReaderType::Pointer readertangOri1X = ImageReaderType::New();
  readertangOri1X->SetFileName("temp/tangOri1X.img");
  readertangOri1X->Update(); 
  ImageReaderType::Pointer readertangOri1Y = ImageReaderType::New();
  readertangOri1Y->SetFileName("temp/tangOri1Y.img");
  readertangOri1Y->Update();
  ImageReaderType::Pointer readertangOri1Z = ImageReaderType::New();
  readertangOri1Z->SetFileName("temp/tangOri1Z.img");
  readertangOri1Z->Update();
  ImageReaderType::Pointer readertangOri2X = ImageReaderType::New();
  readertangOri2X->SetFileName("temp/tangOri2X.img");
  readertangOri2X->Update(); 
  ImageReaderType::Pointer readertangOri2Y = ImageReaderType::New();
  readertangOri2Y->SetFileName("temp/tangOri2Y.img");
  readertangOri2Y->Update();
  ImageReaderType::Pointer readertangOri2Z = ImageReaderType::New();
  readertangOri2Z->SetFileName("temp/tangOri2Z.img");
  readertangOri2Z->Update(); 

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
  VOXEL *tangentOri1;
  tangentOri1 = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);
  VOXEL *tangentOri2;
  tangentOri2 = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);
  float *localSize1;
  localSize1 = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  float *localSize2;
  localSize2 = (float *)malloc(sizeof(float) * pRow * pCol * pSli);
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  DT(i,j,k) = DTImageReader->GetOutput()->GetPixel(pixelIndex);
	  DT(i,j,k) = DT(i,j,k)*DT(i,j,k);
	  tangentOri1(i,j,k).x = readertangOri1X->GetOutput()->GetPixel(pixelIndex);
	  tangentOri1(i,j,k).y = readertangOri1Y->GetOutput()->GetPixel(pixelIndex);
	  tangentOri1(i,j,k).z = readertangOri1Z->GetOutput()->GetPixel(pixelIndex);
	  tangentOri2(i,j,k).x = readertangOri2X->GetOutput()->GetPixel(pixelIndex);
	  tangentOri2(i,j,k).y = readertangOri2Y->GetOutput()->GetPixel(pixelIndex);
	  tangentOri2(i,j,k).z = readertangOri2Z->GetOutput()->GetPixel(pixelIndex);
	}
  
  //Trace structure size
  int maxTrace = 200;
  std::cout<<"Trace structure size in tangent plane: "<<std::endl;
  TraceStructureSize(DT, tangentOri1, localSize1, maxTrace);
  TraceStructureSize(DT, tangentOri2, localSize2, maxTrace);

  //Output structure size
  ImageType::Pointer structSizeImage1 = ImageType::New();
  structSizeImage1->SetRegions( region );
  structSizeImage1->Allocate();
  ImageType::Pointer structSizeImage2 = ImageType::New();  
  structSizeImage2->SetRegions( region );
  structSizeImage2->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  structSizeImage1->SetPixel(pixelIndex,localSize1(i,j,k));
	  structSizeImage2->SetPixel(pixelIndex,localSize2(i,j,k));
      }
  
  ImageWriterType::Pointer writerSize1 = ImageWriterType::New();
  writerSize1->SetInput( structSizeImage1 );
  writerSize1->SetFileName("temp/StructSize1.img");
  writerSize1->Update();
  ImageWriterType::Pointer writerSize2 = ImageWriterType::New();
  writerSize2->SetInput( structSizeImage2 );
  writerSize2->SetFileName("temp/StructSize2.img");
  writerSize2->Update();
  
  free(DT);
  free(tangentOri1);
  free(tangentOri2);
  free(localSize1);
  free(localSize2);
  
  return EXIT_SUCCESS;
}
