//The program compute multi-seed Iterative Rlative Fuzzy Connectedness
#include "FuzzyCon.h"

int main( int argc, char * argv[] )
{
  if( argc < 15 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " seedFile outputImageFile intImageFile vesImageFile scaImageFile intSigma vesSigma scaSigma intMean vesMean scaMean intDarkObjFlag vesDarkObjFlag scaDarkObjFlag" << std::endl;
      return EXIT_FAILURE;
    }
  
  //Filters
  ReaderType::Pointer intReader = ReaderType::New();
  ReaderType::Pointer vesReader = ReaderType::New();
  ReaderType::Pointer scaReader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  intReader->SetFileName( argv[3] );
  vesReader->SetFileName( argv[4] );
  scaReader->SetFileName( argv[5] );
  writer->SetFileName( argv[2] );


  VectorType sigma;
  sigma.intensity = atoi(argv[6]);
  sigma.vesselness = atof(argv[7]);
  sigma.scale = atoi(argv[8]);
  VectorType mean;
  mean.intensity = atoi(argv[9]);
  mean.vesselness = atof(argv[10]);
  mean.scale = atoi(argv[11]);
  VectorType darkObj;
  darkObj.intensity = atoi(argv[12]);
  darkObj.vesselness = atof(argv[13]);
  darkObj.scale = atoi(argv[14]);

  //Pipeline
  try
    {
      intReader->Update();
      vesReader->Update();
      scaReader->Update();
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
	}

  //Read seed
  SeedMapType seedMap;
  SeedLabelListType seedLabelList;
  int x, y, z, index, label;
  std::cout<<std::endl<<"Reading seed file......."<<std::flush;
  FILE *seedFile = fopen(argv[1], "r");
  if(seedFile == NULL) 
    {
      std::cerr << argv[1] << " Seed file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(seedFile, "%d %d %d %d", &x, &y, &z, &label) == 4)  
    {
      index = z*pRow*pCol+y*pRow+x;
      //if label hasn't been recorded, insert to list 
      if(!seedLabelList.count(label))
	seedLabelList.insert(label);
      //insert seed to map
      seedMap.insert(std::pair< unsigned short, unsigned int >(label, index));
    }
  fclose(seedFile);
  std::cout<<"finished"<<std::endl;
  std::cout<<std::endl;
  

  //Seed
  SeedLabelListType::iterator seedLabelIt;
  SeedMapType::iterator seedIt;
  std::pair< SeedMapType::iterator, SeedMapType::iterator > seedRange;

  //Output seed list correspondant to label
  std::cout<<"Seed List:"<<std::endl;
  for(seedLabelIt=seedLabelList.begin(); seedLabelIt!=seedLabelList.end(); ++seedLabelIt)
    {
      label = *seedLabelIt;
      std::cout<<"Label "<<label<<std::endl;
      //Get seed for each label
      seedRange = seedMap.equal_range(label);
      for (seedIt=seedRange.first; seedIt!=seedRange.second; ++seedIt)
	{
	  int seedLabel = (*seedIt).first;
	  int seedIndex = (*seedIt).second;
	  int xN = seedIndex % pRow;
	  int yN = ((seedIndex-xN) / pRow) % pCol;
	  int zN = (seedIndex-xN-yN*pRow) / (pRow*pCol);
	  std::cout<<seedLabel<<": \t "<<xN<<" "<<yN<<" "<<zN<<std::endl;
	}
    }
  std::cout<<std::endl;
  

  //Initialize adjacency relationship
  std::cout<<"Initializing affinity map......."<<std::flush;
  AffinityType *affinityMap;
  affinityMap = (AffinityType *)malloc(sizeof(AffinityType) * pRow * pCol * pSli); 
  computeAffinityMap(affinityMap, inputImage, spacing, sigma, mean, darkObj);
  std::cout<<"finished"<<std::endl;
  std::cout<<std::endl;
  
  //Try multi-seed single-object kFOEMS
  SeedListType seedListS;
  seedLabelIt=seedLabelList.begin();
  label = *seedLabelIt;
  seedRange = seedMap.equal_range(label);
  for(seedIt=seedRange.first; seedIt!=seedRange.second; ++seedIt)
    {
      int seedIndex = (*seedIt).second;
      seedListS.push_front(seedIndex);
    }
  //Compute fuzzy connectedness
  unsigned short *fuzzyConImage;
  fuzzyConImage = (unsigned short *)malloc(sizeof(unsigned short) * pRow * pCol * pSli); 
  memset(fuzzyConImage,0,pRow*pCol*pSli*sizeof(unsigned short));
  //Compute use kFOEMS
  kFOEMS( affinityMap, seedListS, fuzzyConImage);
  std::cout<<std::endl;
  //Set output
  unsigned short value;
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
