#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//ITK settings
const unsigned int Dimension = 3;
typedef unsigned char PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile" << std::endl;
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
  size[0] = size[0] + 2;
  size[1] = size[1] + 2;
  size[2] = size[2] + 2;
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate image for output
  ImageType::Pointer labelImage = ImageType::New();
  labelImage->SetRegions( region );
  labelImage->SetSpacing( spacing );
  labelImage->SetOrigin( origin );
  labelImage->SetDirection( direction );
  labelImage->Allocate();

  ImageType::IndexType pixelIndex;
  int i, j, k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  labelImage->SetPixel(pixelIndex, 0);
	}
  int value;
  for(i=1;i<pRow-1;i++)
    for(j=1;j<pCol-1;j++)
      for(k=1;k<pSli-1;k++)
	{
	  pixelIndex[0]=i-1;
	  pixelIndex[1]=j-1;
	  pixelIndex[2]=k-1;
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  labelImage->SetPixel(pixelIndex, value);
	}

  //Write output  
  writer->SetInput( labelImage );
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
