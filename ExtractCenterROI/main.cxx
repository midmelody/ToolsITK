#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 9 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile centerX centerY centerZ sizeX sizeY sizeZ" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  
  //Parameters
  reader->SetFileName( argv[1] );
  reader->Update();
  writer->SetFileName( argv[2] );
  int centerX = atoi(argv[3]);
  int centerY = atoi(argv[4]);
  int centerZ = atoi(argv[5]);    
  int sizeX = atoi(argv[6]);
  int sizeY = atoi(argv[7]);
  int sizeZ = atoi(argv[8]);
  
  int startX = floor(centerX - sizeX/2 + 1);
  int startY = floor(centerY - sizeY/2 + 1); 
  int startZ = floor(centerZ - sizeZ/2 + 1); 
  int endX = floor(centerX + sizeX/2);
  int endY = floor(centerY + sizeY/2); 
  int endZ = floor(centerZ + sizeZ/2);

  if(sizeX != (endX-startX+1))
   {
      std::cerr << "ROI Size Error!" << std::endl;
      return EXIT_FAILURE;
    }  
    
  //Set image specs
  ImageType::SizeType  size;
  size[0] = sizeX;
  size[1] = sizeY;
  size[2] = sizeZ;
  ImageType::RegionType region;
  region.SetSize( size ); 
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( reader->GetOutput()->GetSpacing() );
  image->SetOrigin( reader->GetOutput()->GetOrigin() );
  image->SetDirection( reader->GetOutput()->GetDirection() );
  image->Allocate();
  image->FillBuffer( 0.0 );

  //Set values
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::IndexType pixelInIndex, pixelOutIndex;
  int i, j, k, value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelInIndex[0] = i+startX;
	  pixelInIndex[1] = j+startY;
	  pixelInIndex[2] = k+startZ;	  
	  value = reader->GetOutput()->GetPixel(pixelInIndex);
	  pixelOutIndex[0] = i;
	  pixelOutIndex[1] = j;
	  pixelOutIndex[2] = k;
	  image->SetPixel(pixelOutIndex, value);
	}

  writer->SetInput( image );
  writer->Update();
  
  return EXIT_SUCCESS;
}
