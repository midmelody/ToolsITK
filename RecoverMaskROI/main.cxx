#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main( int argc, char * argv[] )
{
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " referenceImage maskImage candidateImage outputImage mode(1:mask; 2:above; 3:below; 4:between) extendPercent(mode1)";
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

  int Rmin = 10000000;
  int Rmax = 0;
  int Cmin = 10000000;
  int Cmax = 0;
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
	      if(Rmin > i)
		Rmin = i;
	      if(Rmax < i)
		Rmax = i;
	      if(Cmin > j)
		Cmin = j;
	      if(Cmax < j)
		Cmax = j;
	      if(Smin > k)
		Smin = k;
	      if(Smax < k)
		Smax = k;
	    }
	}

  //if(Rmin > 0) Rmin = Rmin - 1;
  //if(Rmax < pRow-1) Rmax = Rmax + 1;     
  //if(Cmin > 0) Cmin = Cmin - 1;
  //if(Cmax < pCol-1) Cmax = Cmax + 1;
  //if(Smin > 0) Smin = Smin - 1;
  //if(Smax < pSli-1) Smax = Smax + 1;

  int mode = atoi(argv[5]);
  if(mode == 1)
    {
      float ratio = atof(argv[6]);
      int sizeAdj;
      sizeAdj = Rmax - Rmin;
      //std::cout<<sizeAdj<<std::endl;
      sizeAdj = sizeAdj*ratio/2;
      //std::cout<<sizeAdj<<std::endl;
      Rmin = Rmin - sizeAdj;
      Rmax = Rmax + sizeAdj;
      sizeAdj = Cmax - Cmin;
      sizeAdj = sizeAdj*ratio/2;
      Cmin = Cmin - sizeAdj;
      Cmax = Cmax + sizeAdj;
      sizeAdj = Smax - Smin;
      sizeAdj = sizeAdj*ratio/2;
      Smin = Smin - sizeAdj;
      Smax = Smax + sizeAdj;
      if(Rmin > 0) Rmin = Rmin - 1;
      if(Rmax < pRow-1) Rmax = Rmax + 1;     
      if(Cmin > 0) Cmin = Cmin - 1;
      if(Cmax < pCol-1) Cmax = Cmax + 1;
      if(Smin > 0) Smin = Smin - 1;
      if(Smax < pSli-1) Smax = Smax + 1;
    }
  else if(mode == 2)
    {
      Rmin = 0;
      Rmax = pRow-1;
      Cmin = 0;
      Cmax = pCol-1;
      Smin = Smax;
      Smax = pSli-1;
    }
  else if(mode == 3)
    {
      Rmin = 0;
      Rmax = pRow-1;
      Cmin = 0;
      Cmax = pCol-1;
      Smax = Smin;
      Smin = 0;
    }
  else if(mode == 4)
    {
      Rmin = 0;
      Rmax = pRow-1;
      Cmin = 0;
      Cmax = pCol-1;
      Smax = Smax;
      Smin = Smin;
    }
  else
    {
      std::cerr << "Invalid mode" << std::endl;
      return EXIT_FAILURE;
    }

  ImageType::SizeType size;
  size[0] = Rmax - Rmin + 1;
  size[1] = Cmax - Cmin + 1;
  size[2] = Smax - Smin + 1;

  ImageType::SizeType sizeC = candidateImage->GetRequestedRegion().GetSize();

  if((sizeC[0]!= size[0]) || (sizeC[1]!= size[1]) || (sizeC[2]!= size[2]))
    {
      std::cout<<"Candidate dimension doesn't match!"<<sizeC<<";"<<size<<std::endl;
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
      for( inputIndex[1] = 0, outputIndex[1] = Cmin; inputIndex[1] < size[1], outputIndex[1] <= Cmax; inputIndex[1]++, outputIndex[1]++ )
	{
	  for( inputIndex[0] = 0, outputIndex[0] = Rmin; inputIndex[0] < size[0], outputIndex[0] <= Rmax; inputIndex[0]++, outputIndex[0]++ )
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

