#include "itkImage.h"
#include "itkImageFileReader.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFileL inputImageFileS threshold" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  //Filters
  ReaderType::Pointer readerL = ReaderType::New();
  ReaderType::Pointer readerS = ReaderType::New();

  //Parameters
  readerL->SetFileName( argv[1] );
  readerS->SetFileName( argv[2] );
  
  //Pipeline
  readerL->Update();
  readerS->Update();
  ImageType::SizeType size = readerL->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol;
  pRow = size[0];
  pCol = size[1];

  //Correspondance map
  int total = 9;
  float corres[9] = {1, 1, 1.5, 2, 3, 3.5, 4, 5, 5.5}; 
  //float corres[9] = {1, 1.5, 2, 3, 3.5, 4, 5, 5.5, 6}; 
  //float corres[9] = {1, 1, 2, 2.5, 3, 4, 4.5, 5, 6}; 
  //float corres[9] = {1, 2, 2.5, 3, 4, 4.5, 5, 6, 6}; 
  int thre = atoi( argv[3] );

  //Count
  ImageType::IndexType indexL, indexS1, indexS2;
  ImageType::PixelType pixelL, pixelS1, pixelS2;
  int i,j,k;
  double PP, PN, NP, NN;
  PP = 0;
  PN = 0;
  NP = 0;
  NN = 0;
  for(k=0;k<total;k++)
    {
      indexL[2] = k;
      indexS1[2] = floor(corres[k]) - 1;
      indexS2[2] = ceil(corres[k]) - 1;
      std::cout<<indexL[2]+1<<": "<<indexS1[2]+1<<"; "<<indexS2[2]+1<<std::endl;
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {
	    indexL[0] = i;
	    indexL[1] = j;
	    indexS1[0] = i;
	    indexS1[1] = j;
	    indexS2[0] = i;
	    indexS2[1] = j;

	    pixelL = readerL->GetOutput()->GetPixel( indexL );
	    pixelS1 = readerS->GetOutput()->GetPixel( indexS1 );
	    pixelS2 = readerS->GetOutput()->GetPixel( indexS2 );

	    if((pixelL>=thre)&&((pixelS1+pixelS2)/2>=thre))
	      PP++;
	    if((pixelL>=thre)&&((pixelS1+pixelS2)/2<thre))
	      PN++;
	    if((pixelL<thre)&&((pixelS1+pixelS2)/2>=thre))
	      NP++;
	    if((pixelL<thre)&&((pixelS1+pixelS2)/2<thre))
	      NN++;
	  }
    }
  
  std::cout<<"PP "<<PP<<" PN "<<PN<<" NP "<<NP<<" NN "<<NN<<std::endl;

  return EXIT_SUCCESS;
}
