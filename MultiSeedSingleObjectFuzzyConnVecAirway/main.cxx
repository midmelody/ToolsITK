//The program compute multi-seed Iterative Rlative Fuzzy Connectedness
#include "FuzzyCon.h"

int main( int argc, char * argv[] )
{
  if( argc < 14 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " seedFile outputImageFile intImageFile vesImageFile scaImageFile gryImageFile intSigma vesSigma grySigma intMean vesMean gryMean scaThre [maskAll maskLung]" << std::endl;
      return EXIT_FAILURE;
    }
  
  //Filters
  ReaderType::Pointer intReader = ReaderType::New();
  ReaderType::Pointer vesReader = ReaderType::New();
  ReaderType::Pointer scaReader = ReaderType::New();
  ReaderType::Pointer gryReader = ReaderType::New(); 
  WriterType::Pointer writer = WriterType::New();

  ReaderType::Pointer mskReader = ReaderType::New();
  ReaderType::Pointer lngReader = ReaderType::New(); 

  //Parameters
  writer->SetFileName( argv[2] );
  intReader->SetFileName( argv[3] );
  vesReader->SetFileName( argv[4] );
  scaReader->SetFileName( argv[5] );
  gryReader->SetFileName( argv[6] );

  VectorType sigma;
  sigma.intensity = atoi(argv[7]);
  sigma.vesselness = atoi(argv[8]);
  sigma.scale = 0;
  sigma.grayReconst = atoi(argv[9]);

  VectorType mean;
  mean.intensity = atoi(argv[10]);
  mean.vesselness = atoi(argv[11]);
  mean.scale = atoi(argv[13]);
  mean.grayReconst = atoi(argv[12]);

  if(argc == 16)
    { 
      mskReader->SetFileName( argv[14] );
      lngReader->SetFileName( argv[15] );
    }
  
  //Pipeline
  std::cout<<std::endl<<"Reading images......."<<std::flush;
  try
    {
      intReader->Update();
      vesReader->Update();
      scaReader->Update();
      gryReader->Update();
      if(argc == 16)
	{
	  mskReader->Update();
	  lngReader->Update();
	}
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
 
  //Get image specs
  ImageType::SpacingType spacing = intReader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = intReader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = intReader->GetOutput()->GetDirection();
  ImageType::SizeType  size = intReader->GetOutput()->GetRequestedRegion().GetSize();
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
  VectorType *inputImage;
  inputImage = (VectorType *)malloc(sizeof(VectorType) * pRow * pCol * pSli); 
  bool *maskImage;
  maskImage = (bool *)malloc(sizeof(bool) * pRow * pCol * pSli); 
  bool *lungImage;
  lungImage = (bool *)malloc(sizeof(bool) * pRow * pCol * pSli);

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
	  inputImage(i,j,k).intensity = intReader->GetOutput()->GetPixel(pixelIndex);
	  inputImage(i,j,k).vesselness = vesReader->GetOutput()->GetPixel(pixelIndex); 
	  inputImage(i,j,k).scale = scaReader->GetOutput()->GetPixel(pixelIndex); 
	  inputImage(i,j,k).grayReconst = gryReader->GetOutput()->GetPixel(pixelIndex); 
	  if(argc == 16)
	    {
	      maskImage(i,j,k) = mskReader->GetOutput()->GetPixel(pixelIndex); 
	      lungImage(i,j,k) = lngReader->GetOutput()->GetPixel(pixelIndex); 
	    }
	  else
	    {
	      maskImage(i,j,k) = true; 
	      lungImage(i,j,k) = true;
	    }
	}
  std::cout<<"finished"<<std::endl;
  std::cout<<std::endl;
  
  //Read seed
  unsigned int seedCount;
  int x, y, z, label, index;
  std::cout<<"Reading seed file......."<<std::flush;
  FILE *seedFile = fopen(argv[1], "r");
  if(seedFile == NULL) 
    {
      std::cerr << argv[1] << " Seed file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  seedCount = 0;
  while(fscanf(seedFile, "%d %d %d %d\n", &x, &y, &z, &label) == 4)  
    seedCount++;
  fclose(seedFile);

  seedFile = fopen(argv[1], "r");
  unsigned int seed[seedCount];
  i = 0;
  while(fscanf(seedFile, "%d %d %d %d\n", &x, &y, &z, &label) == 4)  
    {  
      index = z*pRow*pCol+y*pRow+x;
      seed[i] = index;
      i++;
    }
  fclose(seedFile);
  std::cout<<"finished"<<std::endl;
  std::cout<<std::endl;
  
  //Seed
  std::cout<<"Seed List:"<<std::endl;
  for(i=0; i<seedCount; i++)
    {
      index = seed[i];
      int xN = index % pRow;
      int yN = ((index-xN) / pRow) % pCol;
      int zN = (index-xN-yN*pRow) / (pRow*pCol);
      std::cout<<xN<<" "<<yN<<" "<<zN<<std::endl;
    }
  std::cout<<std::endl;
  
  //multi-seed single-object kFOEMS
  //Compute fuzzy connectedness
  unsigned int *fuzzyConImage;
  fuzzyConImage = (unsigned int *)malloc(sizeof(unsigned int) * pRow * pCol * pSli); 
  memset(fuzzyConImage,0,pRow*pCol*pSli*sizeof(unsigned int));
  //Compute use kFOEMS
  std::cout <<"Computing kFOEMS......."<<std::flush;
  kFOEMS( inputImage, spacing, sigma, mean, seed, seedCount, fuzzyConImage, maskImage, lungImage);
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
