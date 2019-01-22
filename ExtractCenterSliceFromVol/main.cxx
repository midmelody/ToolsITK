#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <cstdlib>
#include <ctime>
int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
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
    
  //Set image specs and assign output image
  int pRow, pCol, pSli;
  ImageType::SizeType   size = reader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  size[2] = 1;

  ImageType::RegionType region;
  region.SetSize( size ); 
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( reader->GetOutput()->GetSpacing() );
  image->SetOrigin( reader->GetOutput()->GetOrigin() );
  image->SetDirection( reader->GetOutput()->GetDirection() );
  image->Allocate();
  image->FillBuffer( 0.0 );

  ImageType::IndexType pixelIndex;
  int i, j, k, value;
  k = int(pSli/2);
  ImageType::IndexType pixelOutIndex;
  
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      {
	pixelIndex[0] = i;
	pixelIndex[1] = j;
	pixelIndex[2] = k;	  
	value = reader->GetOutput()->GetPixel(pixelIndex);
	pixelOutIndex[0] = i;
	pixelOutIndex[1] = j;
	pixelOutIndex[2] = 0;	
	image->SetPixel(pixelOutIndex, value);
      }

  writer->SetInput(image);
  writer->Update();
  
  return EXIT_SUCCESS;
}
