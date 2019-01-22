#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " referenceImage maskImage candidateImage outputImage ";
      std::cerr << " " << std::endl;
      return EXIT_FAILURE;
    }

  typedef short PixelType;

  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ReaderType::Pointer readerR = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New(); 
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  readerR->SetFileName( argv[1] );
  readerM->SetFileName( argv[2] );
  reader->SetFileName( argv[3] );
  readerR->Update();
  readerM->Update();
  reader->Update();

  ImageType::Pointer referenceImage = ImageType::New();
  ImageType::Pointer maskImage = ImageType::New();
  ImageType::Pointer candidateImage = ImageType::New();
  ImageType::Pointer outputImage = ImageType::New();

  referenceImage = readerR->GetOutput();
  maskImage = readerM->GetOutput();
  candidateImage = reader->GetOutput();

  //Get image specs
  ImageType::SpacingType spacing = referenceImage->GetSpacing(); 
  ImageType::PointType origin = referenceImage->GetOrigin(); 
  ImageType::DirectionType direction = referenceImage->GetDirection();
  ImageType::SizeType  sizeRef = referenceImage->GetRequestedRegion().GetSize();
  int pRow = sizeRef[0];
  int pCol = sizeRef[1];
  int pSli = sizeRef[2]; 
  ImageType::SizeType  sizeM = maskImage->GetRequestedRegion().GetSize();
  if((pRow != sizeM[0]) || (pCol != sizeM[1]) || (pSli != sizeM[2]))
    {
      std::cout<<"Mask dimension doesn't match!"<<std::endl;
      return EXIT_FAILURE;
    }

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
	  value = maskImage->GetPixel(pixelIndex); 
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
  size[2] = Smax - Smin + 1;

  ImageType::SizeType sizeC = candidateImage->GetRequestedRegion().GetSize();
  if((sizeC[0]!= size[0]) || (sizeC[1]!= size[1]) || (sizeC[2]!= size[2]))
    {
      std::cout<<"Candidate dimension doesn't match!"<<std::endl;
      return EXIT_FAILURE;
    }

  outputImage->SetRegions( sizeRef );
  outputImage->SetSpacing( spacing );
  outputImage->SetOrigin( origin );
  outputImage->SetDirection( direction );
  outputImage->Allocate();
  outputImage->FillBuffer( 0.0 );

  ImageType::PixelType pixel;
  ImageType::IndexType inputIndex;
  ImageType::IndexType outputIndex;
  
  for( inputIndex[2] = 0, outputIndex[2] = Smin; inputIndex[2] < size[2], outputIndex[2] <= Smax; inputIndex[2]++, outputIndex[2]++ )
    {
      for( inputIndex[1] = 0, outputIndex[1] = 0; inputIndex[1] < size[1], outputIndex[1] <= size[1]; inputIndex[1]++, outputIndex[1]++ )
	{
	  for( inputIndex[0] = 0, outputIndex[0] = 0; inputIndex[0] < size[0], outputIndex[0] <= size[0]; inputIndex[0]++, outputIndex[0]++ )
	    {
	      pixel = candidateImage->GetPixel( inputIndex );
	      outputImage->SetPixel( outputIndex, pixel );
	    }
	}
    }

  writer->SetFileName( argv[4] );
  writer->SetInput( outputImage );
  writer->Update();

  return EXIT_SUCCESS;
}

