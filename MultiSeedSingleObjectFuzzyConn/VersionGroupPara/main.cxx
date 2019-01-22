//The program compute multi-seed Iterative Rlative Fuzzy Connectedness
#include "FuzzyCon.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile seedFile parameterFile" << std::endl;
      return EXIT_FAILURE;
    }

  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
 
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
  int i, j, k, index, value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  index = k*pRow*pCol+j*pRow+i;
	  inputImage[index] = value;
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
  seedFile = fopen(argv[3], "r");
  for(i=0;i<seedCount;i++)  
    {
      fscanf(seedFile, "%d %d %d %d", &x, &y, &z, &label);
      seed[i].index = z*pRow*pCol+y*pRow+x;
      seed[i].label = label;    
      std::cout<<"seed["<<i+1<<"]: ("<<x<<" "<<y<<" "<<z<<") Label: "<<label<<std::endl;
    }
  fclose(seedFile);
 
  //Read parameter
  LabelPara *parameters;
  int labelCount;
  float meanObj, sigmaObj, confidentBackground;
  int darkObj;
  FILE *paraFile = fopen(argv[4], "r");
  if(paraFile == NULL) 
    {
      std::cerr << argv[4] << " Para file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  labelCount = 0;
  while(fscanf(paraFile, "%d %f %f %f %d", &label, &meanObj, &sigmaObj, &confidentBackground, &darkObj) == 5)  
    {
      labelCount = labelCount + 1;
    }
  fclose(paraFile);
 
  parameters = new LabelPara[labelCount];
  paraFile = fopen(argv[4], "r");
  for(i=0;i<labelCount;i++)  
    {
      fscanf(paraFile, "%d %f %f %f %d", &label, &meanObj, &sigmaObj, &confidentBackground, &darkObj);
      parameters[i].label = label;
      parameters[i].meanObj = meanObj;
      parameters[i].sigmaObj = sigmaObj;
      parameters[i].confidentBackground = confidentBackground;
      if(darkObj!=0)
	parameters[i].darkObj = true;
      else
	parameters[i].darkObj = false;
 	
      std::cout<<"Group Label "<<label<<" with parameters: meanObj "<<meanObj<<" sigmaObj "<<sigmaObj<<" confidentBackground "<<confidentBackground<<" darkObj "<<darkObj<<std::endl;
    }
  fclose(paraFile);
 
  //Initialize adjacency relation map
  AdjacencyMapType adjacencyMap;
  initializeAdjacencyMap(adjacencyMap);
  
  //Check adjacency map example
  std::pair< AdjacencyMapType::iterator, AdjacencyMapType::iterator > rangeDefine;
  //Locate range
  rangeDefine = adjacencyMap.equal_range(seed[0].index);
  //iterate
  AdjacencyMapType::iterator it;
  int count = 0;
  for (it=rangeDefine.first; it!=rangeDefine.second; ++it)
    {
      count++;
      int seedIndex = (*it).first;
      int xS = seedIndex % pRow;
      int yS = ((seedIndex-xS) / pRow) % pCol;
      int zS = (seedIndex-xS-yS*pRow) / (pRow*pCol);
      int neighborIndex = (*it).second;
      int xN = neighborIndex % pRow;
      int yN = ((neighborIndex-xN) / pRow) % pCol;
      int zN = (neighborIndex-xN-yN*pRow) / (pRow*pCol);
      std::cout<<count<<" \t "<<xS<<" "<<yS<<" "<<zS<<"; "<<inputImage[seedIndex]<<" => "<<xN<<" "<<yN<<" "<<zN<<"; "<<inputImage[neighborIndex] << std::endl;
    }

  //Compute fuzzy connectedness 


  /*
  //Compute fuzzy connectedness
  float *fuzzyConImage;
  fuzzyConImage = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  computeAFC(inputImage, fuzzyConImage, seed, seedCount, meanObj, sigmaObj, confidentBackground);

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
  */
  return EXIT_SUCCESS;
}
