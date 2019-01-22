//The program compute multi-seed Iterative Rlative Fuzzy Connectedness
#include "FuzzyCon.h"

int main( int argc, char * argv[] )
{
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile seedFile meanObject sigmaObject confidentBackground " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  float meanObj = atof(argv[4]);
  float sigmaObj = atof(argv[5]);
  float confidentBackground = atof(argv[6]);

  //Pipeline
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
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  

  //Read image
  float *inputImage;
  inputImage = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  ImageType::IndexType pixelIndex;
  int i, j, k, value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  inputImage(i,j,k) = value;
	}

  //Read seed
  IndexLabel *seed;
  int seedCount;
  int x, y, z, label;
  FILE *seedFile = fopen(argv[3], "r");
  if(seedFile == NULL) 
    {
      std::cerr << argv[4] << " Seed file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  seedCount = 0;
  while(fscanf(seedFile, "%d %d %d %d", &x, &y, &z, &label) == 4)  
    {
      seedCount = seedCount + 1;
    }
  fclose(seedFile);
  
  seed = new IndexLabel[seedCount];

  seedFile = fopen(argv[4], "r");
  for(i=0;i<seedCount;i++)  
    {
      fscanf(seedFile, "%d %d %d %d", &x, &y, &z, &label);
      seed[i].index = z*pRow*pCol+y*pRow+x;
      seed[i].label = label;    
      std::cout<<"seed["<<i+1<<"]: "<<x<<" "<<y<<" "<<z<<": "<<label<<std::endl;
    }
  fclose(seedFile);
  
  //Compute fuzzy affinity relation kappa
  typedef std::map< unsigned int, IndexKappa > AffinityMapType;
  AffinityMapType affinityMap;
  computeAffinity(affinityMap, spacing);





  //Compute fuzzy connectedness
  float *fuzzyConImage;
  fuzzyConImage = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  computeFC(inputImage, fuzzyConImage, seed, seedCount, meanObj, sigmaObj, confidentBackground);

  //Set output
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = fuzzyConImage(i,j,k);
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

  return EXIT_SUCCESS;
}
