#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif
#include <stdlib.h>
#include <stdio.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << "  inputImage outputImage SliceMin"<< std::endl;
      return EXIT_FAILURE;
    }
  
  const unsigned int Dimension = 3;
  typedef short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  int sliceMin = atoi(argv[3]);
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  reader->SetFileName ( argv[1] );
  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem reading image file : " << argv[1] << std::endl;
      std::cerr << excp << std::endl;
      return -1;
    }

  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  size[0] = size[0];
  size[1] = size[1];
  size[2] = size[2] +1 - sliceMin;
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate image for output
  ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions( region );
  outImage->SetSpacing( spacing );
  outImage->SetOrigin( origin );
  outImage->SetDirection( direction );
  outImage->Allocate();
  outImage->FillBuffer( 0.0 );
  //Assign value
  ImageType::PixelType pixel;
  ImageType::IndexType indexIn, indexOut;
  int i,j,k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  indexIn[0] = i;
	  indexIn[1] = j;
	  indexIn[2] = k+sliceMin-1;
	  indexOut[0] = i;
	  indexOut[1] = j;
	  indexOut[2] = k;
	  pixel = reader->GetOutput()->GetPixel( indexIn );
	  outImage->SetPixel( indexOut, pixel );
	}

  writer->SetFileName( argv[2] );
  writer->SetInput( outImage );
  writer->Update();
  
  return EXIT_SUCCESS;
}

