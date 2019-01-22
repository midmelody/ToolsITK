#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "MIP3D.h"

int main(int argc, char *argv[])
{
  int i,j,k; //counters

  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage outputImage mode(1:max, 2:avg, 3:min, 4:sum) orientation(1:X 2:Y 3:Z) [Mask]" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK Filters
  ImageReaderType::Pointer Reader = ImageReaderType::New();
  ImageWriterType::Pointer Writer = ImageWriterType::New();
  ImageReaderType::Pointer ReaderM = ImageReaderType::New();

  //Parameters
  Reader->SetFileName( argv[1] );
  Writer->SetFileName( argv[2] );
  int mode = atoi(argv[3]);
  int ori = atoi(argv[4]);
  if((mode!=1)&&(mode!=2)&&(mode!=3)&&(mode!=4))
    {
      std::cerr << "Illegal Mode!"<<std::endl;
      return EXIT_FAILURE;
    }
  if((ori!=1)&&(ori!=2)&&(ori!=3))
    {
      std::cerr << "Illegal Orientation!"<<std::endl;
      return EXIT_FAILURE;
    }
  bool maskFlag = false;
  if(argc == 6)
    {
 
      ReaderM->SetFileName( argv[5] );
      ReaderM->Update();
      maskFlag = true;
    }

  //Get image specs
  Reader->Update();
  ImageType::SpacingType spacing = Reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = Reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = Reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = Reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 

  //Set output
  int ctX, ctY, ctZ;
  if(ori == 1)
    {
      ctX = size[1];
      ctY = size[2];
      ctZ = size[0];
    }
  else if(ori == 2)
    { 
      ctX = size[0];
      ctY = size[2];
      ctZ = size[1];
    }
  else if(ori == 3)
    {
      ctX = size[0];
      ctY = size[1];
      ctZ = size[2];
    }

  size[0] = ctX;
  size[1] = ctY;
  size[2] = 1;
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  

  //Process  
  ImageType::IndexType pixelIndex;
  ImageType::IndexType outIndex;
  float value;    
  int valueMsk;  
  float valueMaxGen, valueMaxMsk, valueMax;
  float valueAvgGen, valueAvgMsk, valueAvg;
  float valueMinGen, valueMinMsk, valueMin;
  float valueCtGen, valueCtMsk, valueCt;
  for(i=0;i<ctX;i++)
    for(j=0;j<ctY;j++)
      {
	valueMaxGen = -10000000;
	valueAvgGen = 0;
	valueMinGen = 10000000;
	valueCtGen = 0;
	valueMaxMsk = -10000000;
	valueAvgMsk = 0;
	valueMinMsk = 10000000;
	valueCtMsk = 0;

	for(k=0;k<ctZ;k++)
	  {
	    if(ori == 1)
	      {
		pixelIndex[0]=k;
		pixelIndex[1]=i;
		pixelIndex[2]=j;
	      }
	    else if(ori == 2)
	      {
		pixelIndex[0]=i;
		pixelIndex[1]=k;
		pixelIndex[2]=j;
	      }
	    else if(ori == 3)
	      {
		pixelIndex[0]=i;
		pixelIndex[1]=j;
		pixelIndex[2]=k;
	      }
	    
	    value = Reader->GetOutput()->GetPixel(pixelIndex);
	    //General max, avg, and min
	    if(valueMaxGen<value)
	      valueMaxGen=value;
	    valueAvgGen = valueAvgGen+value;
	    if(valueMinGen>value)
	      valueMinGen=value;
	    valueCtGen = valueCtGen+1;
	    //Masked case, only count masked region
	    if(maskFlag)
	      {
		valueMsk = ReaderM->GetOutput()->GetPixel(pixelIndex);
		if(valueMsk!=0)
		  {
		    if(valueMaxMsk<value)
		      valueMaxMsk=value;
		    valueAvgMsk = valueAvgMsk+value;
		    if(valueMinMsk>value)
		      valueMinMsk=value;
		    valueCtMsk = valueCtMsk+1;		    
		  }
	      }
	  }
	//if mask and value exist, assign mask value
	if((maskFlag)&&(valueCtMsk))
	  {
	    valueMax = valueMaxMsk;
	    valueAvg = valueAvgMsk;
	    valueMin = valueMinMsk;
	    valueCt = valueCtMsk;
	  }
	else
	  {
	    valueMax = valueMaxGen;
	    valueAvg = valueAvgGen;
	    valueMin = valueMinGen;
	    valueCt = valueCtGen;
	  }

	outIndex[0] = i;
	outIndex[1] = j;
	outIndex[2] = 0;

	if(mode==1)
	  image->SetPixel(outIndex, valueMax);
	else if(mode==2)
	  image->SetPixel(outIndex, valueAvg/valueCt);
	else if(mode==3)
	  image->SetPixel(outIndex, valueMin);
	else if(mode==4)
	  image->SetPixel(outIndex, valueAvg);
      }

  //Write output 
 
  Writer->SetInput( image );
  try
    {
      Writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }
 
  return EXIT_SUCCESS;
}


