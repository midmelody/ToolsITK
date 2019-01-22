//The program compute multi-seed Iterative Rlative Fuzzy Connectedness
#include "combine.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " intensityImageFile vesselnessImageFile scaleImageFile outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }

  
  //Filters
  ReaderType::Pointer intReader = ReaderType::New();
  ReaderType::Pointer vesReader = ReaderType::New();
  ReaderType::Pointer scaReader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  intReader->SetFileName( argv[1] );
  vesReader->SetFileName( argv[2] );
  scaReader->SetFileName( argv[3] );
  writer->SetFileName( argv[4] );
 
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
  OutImageType::SpacingType spacing = intReader->GetOutput()->GetSpacing(); 
  OutImageType::PointType origin = intReader->GetOutput()->GetOrigin(); 
  OutImageType::DirectionType direction = intReader->GetOutput()->GetDirection();
  OutImageType::SizeType  size = intReader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  OutImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  OutImageType::Pointer image = OutImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  

  //Read image
  OutImageType::IndexType pixelIndex;
  OutImageType::PixelType pixelValue;
  int i, j, k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  pixelValue[0] = intReader->GetOutput()->GetPixel(pixelIndex);
	  pixelValue[1] = vesReader->GetOutput()->GetPixel(pixelIndex);
	  pixelValue[2] = scaReader->GetOutput()->GetPixel(pixelIndex);
	  image->SetPixel(pixelIndex, pixelValue);
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
