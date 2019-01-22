#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main( int argc, char * argv[] )
{
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage maskImage outputImage mode(1:above; 2:below; 3:between; 4:mask/Ratio; 5:mask/Fixed) extendRatio/length(mode4/5) SquareFlag(mode4)";
      std::cerr << " " << std::endl;
      return EXIT_FAILURE;
    }

  typedef int PixelType;

  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New(); 
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  readerM->SetFileName( argv[2] );
  reader->Update();
  readerM->Update();

  ImageType::Pointer inputImage = ImageType::New();
  ImageType::Pointer maskImage = ImageType::New();
  ImageType::Pointer outputImage = ImageType::New();

  inputImage = reader->GetOutput();
  maskImage = readerM->GetOutput();

  ImageType::SizeType sizeOri = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = sizeOri[0];
  int pCol = sizeOri[1];
  int pSli = sizeOri[2];
  //std::cout<<pRow<<" "<<pCol<<" "<<pSli<<std::endl;
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
  float spa = spacing[0];
  
  int Rmin = 10000000;
  int Rmax = 0;
  int Cmin = 10000000;
  int Cmax = 0;
  int Smin = 10000000;
  int Smax = 0;

  int i, j, k, value;
  ImageType::IndexType pixelIndex;
  int ct = 0;
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
	      ct++;
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

  //std::cout<<Rmin<<" "<<Rmax<<" "<<Cmin<<" "<<Cmax<<" "<<Smin<<" "<<Smax<<std::endl;
  //Check if mask exist
  if(ct==0)
    {
      std::cerr << "Blank mask " << std::endl;
      return EXIT_FAILURE;
    }    
  
  //if(Rmin > 0) Rmin = Rmin - 1;
  //if(Rmax < pRow-1) Rmax = Rmax + 1;     
  //if(Cmin > 0) Cmin = Cmin - 1;
  //if(Cmax < pCol-1) Cmax = Cmax + 1;
  //if(Smin > 0) Smin = Smin - 1;
  //if(Smax < pSli-1) Smax = Smax + 1;

  //std::cout<<Rmin<<" "<<Rmax<<" "<<Cmin<<" "<<Cmax<<" "<<Smin<<" "<<Smax<<std::endl;

  int mode = atoi(argv[4]);
  if(mode == 1)
    {
      Rmin = 0;
      Rmax = pRow-1;
      Cmin = 0;
      Cmax = pCol-1;
      Smin = Smax;
      Smax = pSli-1;
    }
  else if(mode == 2)
    {
      Rmin = 0;
      Rmax = pRow-1;
      Cmin = 0;
      Cmax = pCol-1;
      Smax = Smin;
      Smin = 0;
    }
  else if(mode == 3)
    {
      Rmin = 0;
      Rmax = pRow-1;
      Cmin = 0;
      Cmax = pCol-1;
      Smax = Smax;
      Smin = Smin;
    }
  else if((mode == 4) || (mode == 5))
    {
      float ratio;
      float sizeAdj, sizeAdjR, sizeAdjC, sizeAdjS;
      int upBnd;
      if(mode == 4)
	{
	  ratio = atof(argv[5]);
	  bool sqrFlag = atoi(argv[6]);
	  //std::cout<<sqrFlag<<": ";
	  //Extend the ROI range, limited by the image size
	  int centR, centC, centS;
	  centR = (Rmin+Rmax)/2;
	  centC = (Cmin+Cmax)/2;
	  centS = (Smin+Smax)/2;
	  
	  sizeAdjR = Rmax-Rmin+1;
	  sizeAdjR = sizeAdjR*ratio/2;
	  upBnd = pRow-1-Rmax;
	  if(upBnd > Rmin)
	    upBnd = Rmin;
	  if(sizeAdjR>upBnd)
	    sizeAdjR = upBnd;

	  sizeAdjC = Cmax-Cmin+1;
	  sizeAdjC = sizeAdjC*ratio/2;
	  upBnd = pCol-1-Cmax;
	  if(upBnd > Cmin)
	    upBnd = Cmin;
	  if(sizeAdjC>upBnd)
	    sizeAdjC = upBnd;

	  sizeAdjS = Smax-Smin+1;
	  sizeAdjS = sizeAdjS*ratio/2;
	  upBnd = pSli-1-Smax;
	  if(upBnd > Smin)
	    upBnd = Smin;
	  if(sizeAdjS>upBnd)
	    sizeAdjS = upBnd;
	  if(sqrFlag)
	    {
	      if(sizeAdjR > sizeAdjC)
		sizeAdj = sizeAdjR;
	      else
		sizeAdj = sizeAdjC;
	      if(sizeAdj < sizeAdjS)
		sizeAdj = sizeAdjS;
	      //std::cout<<sizeAdj<<": ";
	      Rmin = centR - sizeAdj;
	      Rmax = Rmin + sizeAdj*2;
	      Cmin = centC - sizeAdj;
	      Cmax = Cmin + sizeAdj*2;
	      Smin = centS - sizeAdj;
	      Smax = Smin + sizeAdj*2;
	    }
	  else
	    {
	      Rmin = centR - sizeAdjR;
	      Rmax = Rmin + sizeAdjR*2;
	      Cmin = centC - sizeAdjC;
	      Cmax = Cmin + sizeAdjC*2;
	      Smin = centS - sizeAdjS;
	      Smax = Smin + sizeAdjS*2;
	    }
	  if(Rmin < 0) Rmin = 0;
	  if(Rmax > pRow-1) Rmax = pRow-1;
	  if(Cmin < 0) Cmin = 0;
	  if(Cmax > pCol-1) Cmax = pCol-1;
	  if(Smin < 0) Smin = 0;
	  if(Smax > pSli-1) Smax = pSli-1;
	}
      else
	{
	  sizeAdj = floor(atof(argv[5])/spa);
	  //std::cout<<argv[5]<<" "<<spa<<" "<<sizeAdj<<std::endl;
	  int centR, centC, centS;
	  centR = (Rmin+Rmax)/2;
	  centC = (Cmin+Cmax)/2;
	  centS = (Smin+Smax)/2;  
	  Rmin = centR - sizeAdj;
	  Rmax = centR + sizeAdj;
	  Cmin = centC - sizeAdj;
	  Cmax = centC + sizeAdj;
	  Smin = centS - sizeAdj;
	  Smax = centS + sizeAdj;
	  if(Rmin < 0) Rmin = 0;
	  if(Rmax > pRow-1) Rmax = pRow-1;
	  if(Cmin < 0) Cmin = 0;
	  if(Cmax > pCol-1) Cmax = pCol-1;
	  if(Smin < 0) Smin = 0;
	  if(Smax > pSli-1) Smax = pSli-1;
	}
  
    }

  else
    {
      std::cerr << "Invalid mode" << std::endl;
      return EXIT_FAILURE;
    }

  //std::cout<<Rmin<<" "<<Rmax<<" "<<Cmin<<" "<<Cmax<<" "<<Smin<<" "<<Smax<<std::endl;

  ImageType::SizeType size;
  size[0] = Rmax - Rmin + 1;
  size[1] = Cmax - Cmin + 1;
  size[2] = Smax - Smin + 1;
  //std::cout<<size<<std::endl;
  outputImage->SetRegions( size );
  outputImage->SetSpacing( inputImage->GetSpacing() );
  outputImage->SetOrigin( inputImage->GetOrigin() );
  outputImage->SetDirection( inputImage->GetDirection() );
  outputImage->Allocate(); 
  outputImage->FillBuffer( 0.0 );

  ImageType::PixelType pixel;
  ImageType::IndexType inputIndex;
  ImageType::IndexType outputIndex;
  
  for( outputIndex[2] = 0, inputIndex[2] = Smin; outputIndex[2] < size[2], inputIndex[2] <= Smax; outputIndex[2]++, inputIndex[2]++ )
    {
      for( outputIndex[1] = 0, inputIndex[1] = Cmin; outputIndex[1] < size[1], inputIndex[1] <= Cmax; outputIndex[1]++, inputIndex[1]++ )
	{
	  for( outputIndex[0] = 0, inputIndex[0] = Rmin; outputIndex[0] < size[0], inputIndex[0] <= Rmax; outputIndex[0]++, inputIndex[0]++ )
	    {
	      pixel = inputImage->GetPixel( inputIndex );
	      outputImage->SetPixel( outputIndex, pixel );
	    }
	}
    }

  writer->SetFileName( argv[3] );
  writer->SetInput( outputImage );
  writer->Update();

  return EXIT_SUCCESS;
}

