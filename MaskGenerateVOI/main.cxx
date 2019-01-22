#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " maskImage outputImage ";
      std::cerr << " " << std::endl;
      return EXIT_FAILURE;
    }

  typedef short PixelType;

  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::Pointer inputImage = ImageType::New();
  ImageType::Pointer maskImage = ImageType::New();
  ImageType::Pointer outputImage = ImageType::New();

  inputImage = reader->GetOutput();

  ImageType::SizeType sizeOri = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = sizeOri[0];
  int pCol = sizeOri[1];
  int pSli = sizeOri[2]; 

  int Smin = 10000000;
  int Smax = 0;

  int i, j, k, value;
  ImageType::IndexType pixelIndex;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = inputImage->GetPixel(pixelIndex); 
	  if(value>0)
	    {
	      if(Smin > k)
		Smin = k;
	      if(Smax < k)
		Smax = k;
	    }
	}

  if(Smin > 0) Smin = Smin - 1;
  if(Smax < pSli-1) Smax = Smax + 1;

  ImageType::SizeType size;
  size[0] = pRow;
  size[1] = pCol;
  size[2] = pSli;

  outputImage->SetRegions( size );
  outputImage->SetSpacing( inputImage->GetSpacing() );
  outputImage->SetOrigin( inputImage->GetOrigin() );
  outputImage->SetDirection( inputImage->GetDirection() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0.0 );
  
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=Smin;k<=Smax;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  outputImage->SetPixel( pixelIndex, 1 );
	}

  writer->SetFileName( argv[2] );
  writer->SetInput( outputImage );
  writer->Update();

  return EXIT_SUCCESS;
}

