//The program compute multi-seed Iterative Rlative Fuzzy Connectedness
#include "FuzzyCon.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " seedFile outputImageFile inputImageFile sigma "<< std::endl;
      return EXIT_FAILURE;
    }
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  writer->SetFileName( argv[2] );
  reader->SetFileName( argv[3] );

  int sigma = atoi(argv[4]);
 
  //Pipeline
  std::cout<<std::endl<<"Reading images......."<<std::flush;
  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
 
  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  OutImageType::Pointer image = OutImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  

  //Read image
  int *inputImage;
  inputImage = (   int *)malloc(sizeof(   int ) * pRow * pCol * pSli); 

  ImageType::IndexType pixelIndex;
  ImageType::PixelType pixelValue;
  int i, j, k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  inputImage(i,j,k) = reader->GetOutput()->GetPixel(pixelIndex);
	}
  std::cout<<"finished"<<std::endl;
  std::cout<<std::endl;
  
  //Read seed
    int seedCount;
    int x, y, z, index, lable;
  std::cout<<"Reading seed file......."<<std::flush;
  FILE *seedFile = fopen(argv[1], "r");
  if(seedFile == NULL) 
    {
      std::cerr << argv[1] << " Seed file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  seedCount = 0;
  int lab;
  while(fscanf(seedFile, "%d %d %d %d", &x, &y, &z, &lable) == 4)  
    seedCount++;
  fclose(seedFile);

  if(seedCount)
    {
      seedFile = fopen(argv[1], "r");
      int seed[seedCount];
      i = 0;
      while(fscanf(seedFile, "%d %d %d %d", &x, &y, &z, &lable) == 4)  
	{  
	  index = z*pRow*pCol+y*pRow+x;
	  seed[i] = index;
	  i++;
	}
      fclose(seedFile);
      std::cout<<"finished"<<std::endl;
      std::cout<<std::endl;
      /*
      //Seed
      std::cout<<"Seed List:"<<std::endl;
      for(i=0; i<seedCount; i++)
      {
      index = seed[i];
      int xN = index % pRow;
      int yN = ((index-xN) / pRow) % pCol;
      int zN = (index-xN-yN*pRow) / (pRow*pCol);
      std::cout<<xN<<" "<<yN<<" "<<zN<<": "<<inputImage(xN,yN,zN)<<std::endl;
      }
      std::cout<<std::endl;
      */
      //multi-seed single-object kFOEMS
      //Compute fuzzy connectedness
      int *fuzzyConImage;
      fuzzyConImage = (  int *)malloc(sizeof(  int) * pRow * pCol * pSli); 
      memset(fuzzyConImage,0,pRow*pCol*pSli*sizeof(  int));
      //Compute use kFOEMS
      std::cout <<"Computing kFOEMS......."<<std::flush;
      kFOEMS( inputImage, spacing, sigma, seed, seedCount, fuzzyConImage);
      std::cout <<"finished"<<std::endl;
      std::cout<<std::endl;
      //Set output
      std::cout <<"Writing output......."<<std::flush;
      float value;
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  for(k=0;k<pSli;k++)
	    {
	      pixelIndex[0]=i;
	      pixelIndex[1]=j;
	      pixelIndex[2]=k;
	      value = float(fuzzyConImage(i,j,k));
	      image->SetPixel(pixelIndex, value);
	    }
    }
  else
    {
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  for(k=0;k<pSli;k++)
	    {
	      pixelIndex[0]=i;
	      pixelIndex[1]=j;
	      pixelIndex[2]=k;
	      image->SetPixel(pixelIndex, 0);
	    }
    }

  //Write output  
  writer->SetInput( image );
  try
    {
      writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }
  std::cout <<"finished"<<std::endl;
  std::cout<<std::endl;
  return EXIT_SUCCESS;
}
