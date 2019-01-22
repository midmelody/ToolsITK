//The program compute multi-seed Iterative Rlative Fuzzy Connectedness
#include "FuzzyCon.h"

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile ROIImageFile ObjectMaskImageFile outputImageFile sigma " << std::endl;
      return EXIT_FAILURE;
    }

  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerROI = ReaderType::New();
  ReaderType::Pointer readerMSK = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  readerROI->SetFileName( argv[2] );
  readerMSK->SetFileName( argv[3] );
  writer->SetFileName( argv[4] );
  float sigma = atof(argv[5]);

  //Pipeline
  try
    {
      reader->Update();
      readerROI->Update();
      readerMSK->Update();
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
  bool *ROIImage, *MSKImage;
  ROIImage = (bool *)malloc(sizeof(bool) * pRow * pCol * pSli); 
  MSKImage = (bool *)malloc(sizeof(bool) * pRow * pCol * pSli);
 
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
	  if(value<-600)
	    inputImage(i,j,k) = -600;
	  value = readerROI->GetOutput()->GetPixel(pixelIndex);
	  ROIImage(i,j,k) = value;
	  value = readerMSK->GetOutput()->GetPixel(pixelIndex);
	  MSKImage(i,j,k) = value;
	}

  //Generate Seed
  SeedMapType seedMap;
  SeedLabelListType seedLabelList;
  seedLabelList.insert(1);
  seedLabelList.insert(2);
  bool ROIValue, MSKValue;
  int index, label;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  index = k*pRow*pCol+j*pRow+i;
	  if(MSKImage(i,j,k)&&(ROIImage(i,j,k)))
	    {
	      label = 1;
	      seedMap.insert(std::pair< unsigned short, unsigned int >(label, index));	      
	    }
	  else if((ROIImage(i,j,k))&&(inputImage(i,j,k)>150)&&(inputImage(i,j,k)<200))
	    {
	      label = 2;
	      seedMap.insert(std::pair< unsigned short, unsigned int >(label, index));	      
	    }	    
	}
 
  //Initialize adjacency relationship
  std::cout<<"Initializing affinity map.......";
  AffinityType *affinityMap;
  affinityMap = (AffinityType *)malloc(sizeof(AffinityType) * pRow * pCol * pSli); 
  computeAffinityMap(affinityMap, inputImage, spacing, sigma, ROIImage);
  std::cout<<"finished"<<std::endl;
  std::cout<<std::endl;

  //Try multi-seed multi-object 
  unsigned short *labelSegImage;
  labelSegImage = (unsigned short *)malloc(sizeof(unsigned short) * pRow * pCol * pSli); 
  memset(labelSegImage,0,pRow*pCol*pSli*sizeof(unsigned short));
  //Compute use kIRMOFC
  kIRMOFC( affinityMap, seedMap, seedLabelList, labelSegImage, ROIImage);

  //Set output
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = labelSegImage(i,j,k);
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
