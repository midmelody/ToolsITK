#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "Write.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters


  //Read image
  ImageReaderType::Pointer Size1Reader = ImageReaderType::New();
  Size1Reader->SetFileName( "temp/DT.img" );
  Size1Reader->Update();
  ImageReaderType::Pointer Ori1XReader = ImageReaderType::New();
  Ori1XReader->SetFileName( "temp/GradientDTX.img" );
  Ori1XReader->Update();
  ImageReaderType::Pointer Ori1YReader = ImageReaderType::New();
  Ori1YReader->SetFileName( "temp/GradientDTY.img" );
  Ori1YReader->Update();
  ImageReaderType::Pointer Ori1ZReader = ImageReaderType::New();
  Ori1ZReader->SetFileName( "temp/GradientDTZ.img" );
  Ori1ZReader->Update();
  ImageReaderType::Pointer Size2Reader = ImageReaderType::New();
  Size2Reader->SetFileName( "temp/StructSize1.img" );
  Size2Reader->Update();
  ImageReaderType::Pointer Ori2XReader = ImageReaderType::New();
  Ori2XReader->SetFileName( "temp/tangOri1X.img" );
  Ori2XReader->Update();
  ImageReaderType::Pointer Ori2YReader = ImageReaderType::New();
  Ori2YReader->SetFileName( "temp/tangOri1Y.img" );
  Ori2YReader->Update();
  ImageReaderType::Pointer Ori2ZReader = ImageReaderType::New();
  Ori2ZReader->SetFileName( "temp/tangOri1Z.img" );
  Ori2ZReader->Update();
  ImageReaderType::Pointer Size3Reader = ImageReaderType::New();
  Size3Reader->SetFileName( "temp/StructSize2.img" );
  Size3Reader->Update();
  ImageReaderType::Pointer Ori3XReader = ImageReaderType::New();
  Ori3XReader->SetFileName( "temp/tangOri2X.img" );
  Ori3XReader->Update();
  ImageReaderType::Pointer Ori3YReader = ImageReaderType::New();
  Ori3YReader->SetFileName( "temp/tangOri2Y.img" );
  Ori3YReader->Update();
  ImageReaderType::Pointer Ori3ZReader = ImageReaderType::New();
  Ori3ZReader->SetFileName( "temp/tangOri2Z.img" );
  Ori3ZReader->Update();
  ImageReaderType::Pointer LocSizeReader = ImageReaderType::New();
  LocSizeReader->SetFileName( "temp/LocSize.img" );
  LocSizeReader->Update();

  //Get image specs
  ImageType::SpacingType spacing = Size1Reader->GetOutput()->GetSpacing();
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];
  ImageType::SizeType  size = Size1Reader->GetOutput()->GetRequestedRegion().GetSize();
  size[0] = size[0] - 1;
  size[1] = size[1] - 1;
  size[2] = size[2] - 1;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::RegionType region;
  ImageType::IndexType start = Size1Reader->GetOutput()->GetRequestedRegion().GetIndex();
  region.SetSize( size );
  region.SetIndex( start );
 
  //Assign image
  ImageType::Pointer Size1Image = ImageType::New();
  Size1Image->SetRegions( region );
  Size1Image->SetSpacing( spacing );
  Size1Image->Allocate();
  ImageType::Pointer Ori1XImage = ImageType::New();
  Ori1XImage->SetRegions( region );
  Ori1XImage->SetSpacing( spacing );
  Ori1XImage->Allocate();
  ImageType::Pointer Ori1YImage = ImageType::New();
  Ori1YImage->SetRegions( region );
  Ori1YImage->SetSpacing( spacing );
  Ori1YImage->Allocate();
  ImageType::Pointer Ori1ZImage = ImageType::New();
  Ori1ZImage->SetRegions( region );
  Ori1ZImage->SetSpacing( spacing );
  Ori1ZImage->Allocate();
  ImageType::Pointer Size2Image = ImageType::New();
  Size2Image->SetRegions( region );
  Size2Image->SetSpacing( spacing );
  Size2Image->Allocate();
  ImageType::Pointer Ori2XImage = ImageType::New();
  Ori2XImage->SetRegions( region );
  Ori2XImage->SetSpacing( spacing );
  Ori2XImage->Allocate();
  ImageType::Pointer Ori2YImage = ImageType::New();
  Ori2YImage->SetRegions( region );
  Ori2YImage->SetSpacing( spacing );
  Ori2YImage->Allocate();
  ImageType::Pointer Ori2ZImage = ImageType::New();
  Ori2ZImage->SetRegions( region );
  Ori2ZImage->SetSpacing( spacing );
  Ori2ZImage->Allocate();
  ImageType::Pointer Size3Image = ImageType::New();
  Size3Image->SetRegions( region );
  Size3Image->SetSpacing( spacing );
  Size3Image->Allocate();
  ImageType::Pointer Ori3XImage = ImageType::New();
  Ori3XImage->SetRegions( region );
  Ori3XImage->SetSpacing( spacing );
  Ori3XImage->Allocate();
  ImageType::Pointer Ori3YImage = ImageType::New();
  Ori3YImage->SetRegions( region );
  Ori3YImage->SetSpacing( spacing );
  Ori3YImage->Allocate();
  ImageType::Pointer Ori3ZImage = ImageType::New();
  Ori3ZImage->SetRegions( region );
  Ori3ZImage->SetSpacing( spacing );
  Ori3ZImage->Allocate();
  ImageType::Pointer LocSizeImage = ImageType::New();
  LocSizeImage->SetRegions( region );
  LocSizeImage->SetSpacing( spacing );
  LocSizeImage->Allocate();
  ImageType::IndexType pixelIndex;
  float size1, size2, size3;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  Size1Image->SetPixel(pixelIndex,Size1Reader->GetOutput()->GetPixel(pixelIndex)*Size1Reader->GetOutput()->GetPixel(pixelIndex));
	  Ori1XImage->SetPixel(pixelIndex,Ori1XReader->GetOutput()->GetPixel(pixelIndex));
	  Ori1YImage->SetPixel(pixelIndex,Ori1YReader->GetOutput()->GetPixel(pixelIndex));
	  Ori1ZImage->SetPixel(pixelIndex,Ori1ZReader->GetOutput()->GetPixel(pixelIndex));
	  size1 = Size1Reader->GetOutput()->GetPixel(pixelIndex)*Size1Reader->GetOutput()->GetPixel(pixelIndex);
	  size2 = Size2Reader->GetOutput()->GetPixel(pixelIndex);
	  size3 = Size3Reader->GetOutput()->GetPixel(pixelIndex);
	  if(size2<size1) size2 = size1;
	  if(size3<size1) size3 = size1;
	  if(size2<=size3)
	    {
	      Size2Image->SetPixel(pixelIndex,Size2Reader->GetOutput()->GetPixel(pixelIndex));
	      Ori2XImage->SetPixel(pixelIndex,Ori2XReader->GetOutput()->GetPixel(pixelIndex));
	      Ori2YImage->SetPixel(pixelIndex,Ori2YReader->GetOutput()->GetPixel(pixelIndex));
	      Ori2ZImage->SetPixel(pixelIndex,Ori2ZReader->GetOutput()->GetPixel(pixelIndex));
	      Size3Image->SetPixel(pixelIndex,Size3Reader->GetOutput()->GetPixel(pixelIndex));
	      Ori3XImage->SetPixel(pixelIndex,Ori3XReader->GetOutput()->GetPixel(pixelIndex));
	      Ori3YImage->SetPixel(pixelIndex,Ori3YReader->GetOutput()->GetPixel(pixelIndex));
	      Ori3ZImage->SetPixel(pixelIndex,Ori3ZReader->GetOutput()->GetPixel(pixelIndex));
	    }
	  else
	    {
	      Size2Image->SetPixel(pixelIndex,Size3Reader->GetOutput()->GetPixel(pixelIndex));
	      Ori2XImage->SetPixel(pixelIndex,Ori3XReader->GetOutput()->GetPixel(pixelIndex));
	      Ori2YImage->SetPixel(pixelIndex,Ori3YReader->GetOutput()->GetPixel(pixelIndex));
	      Ori2ZImage->SetPixel(pixelIndex,Ori3ZReader->GetOutput()->GetPixel(pixelIndex));
	      Size3Image->SetPixel(pixelIndex,Size2Reader->GetOutput()->GetPixel(pixelIndex));
	      Ori3XImage->SetPixel(pixelIndex,Ori2XReader->GetOutput()->GetPixel(pixelIndex));
	      Ori3YImage->SetPixel(pixelIndex,Ori2YReader->GetOutput()->GetPixel(pixelIndex));
	      Ori3ZImage->SetPixel(pixelIndex,Ori2ZReader->GetOutput()->GetPixel(pixelIndex));	      
	    }
	  LocSizeImage->SetPixel(pixelIndex,LocSizeReader->GetOutput()->GetPixel(pixelIndex));
	}

  //output 
  ImageWriterType::Pointer Size1Writer = ImageWriterType::New();
  Size1Writer->SetInput( Size1Image );
  Size1Writer->SetFileName("TensorScale/Size1.img");
  Size1Writer->Update();
  ImageWriterType::Pointer Ori1XWriter = ImageWriterType::New();
  Ori1XWriter->SetInput( Ori1XImage );
  Ori1XWriter->SetFileName("TensorScale/Ori1X.img");
  Ori1XWriter->Update(); 
  ImageWriterType::Pointer Ori1YWriter = ImageWriterType::New();
  Ori1YWriter->SetInput( Ori1YImage );
  Ori1YWriter->SetFileName("TensorScale/Ori1Y.img");
  Ori1YWriter->Update();
  ImageWriterType::Pointer Ori1ZWriter = ImageWriterType::New();
  Ori1ZWriter->SetInput( Ori1ZImage );
  Ori1ZWriter->SetFileName("TensorScale/Ori1Z.img");
  Ori1ZWriter->Update();
  ImageWriterType::Pointer Size2Writer = ImageWriterType::New();
  Size2Writer->SetInput( Size2Image );
  Size2Writer->SetFileName("TensorScale/Size2.img");
  Size2Writer->Update();
  ImageWriterType::Pointer Ori2XWriter = ImageWriterType::New();
  Ori2XWriter->SetInput( Ori2XImage );
  Ori2XWriter->SetFileName("TensorScale/Ori2X.img");
  Ori2XWriter->Update(); 
  ImageWriterType::Pointer Ori2YWriter = ImageWriterType::New();
  Ori2YWriter->SetInput( Ori2YImage );
  Ori2YWriter->SetFileName("TensorScale/Ori2Y.img");
  Ori2YWriter->Update();
  ImageWriterType::Pointer Ori2ZWriter = ImageWriterType::New();
  Ori2ZWriter->SetInput( Ori2ZImage );
  Ori2ZWriter->SetFileName("TensorScale/Ori2Z.img");
  Ori2ZWriter->Update();
  ImageWriterType::Pointer Size3Writer = ImageWriterType::New();
  Size3Writer->SetInput( Size3Image );
  Size3Writer->SetFileName("TensorScale/Size3.img");
  Size3Writer->Update();
  ImageWriterType::Pointer Ori3XWriter = ImageWriterType::New();
  Ori3XWriter->SetInput( Ori3XImage );
  Ori3XWriter->SetFileName("TensorScale/Ori3X.img");
  Ori3XWriter->Update(); 
  ImageWriterType::Pointer Ori3YWriter = ImageWriterType::New();
  Ori3YWriter->SetInput( Ori3YImage );
  Ori3YWriter->SetFileName("TensorScale/Ori3Y.img");
  Ori3YWriter->Update();
  ImageWriterType::Pointer Ori3ZWriter = ImageWriterType::New();
  Ori3ZWriter->SetInput( Ori3ZImage );
  Ori3ZWriter->SetFileName("TensorScale/Ori3Z.img");
  Ori3ZWriter->Update();
  ImageWriterType::Pointer LocSizeWriter = ImageWriterType::New();
  LocSizeWriter->SetInput( LocSizeImage );
  LocSizeWriter->SetFileName("TensorScale/LocSize.img");
  LocSizeWriter->Update();
  return EXIT_SUCCESS;
}
