#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile maxSizeX maxSizeY maxSizeZ" << std::endl;
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
  int sizeMaxX = atoi(argv[3]);
  int sizeMaxY = atoi(argv[4]);
  int sizeMaxZ = atoi(argv[5]);
    
  //Get image specs
  ImageType::SizeType size = reader->GetOutput()->GetRequestedRegion().GetSize(); 
  int sizeOriX = size[0];
  int sizeOriY = size[1];
  int sizeOriZ = size[2];
  if((sizeOriX<sizeMaxX)&&(sizeOriZ<sizeMaxZ)&&(sizeOriZ<sizeMaxZ))
    {
      writer->SetInput( reader->GetOutput() );
    }
  else
    {
      int startX, startY, startZ;
      int endX, endY, endZ; 
      int sizeX, sizeY, sizeZ;
      
      if(sizeOriX>sizeMaxX)
	{
	  startX = sizeOriX/2-sizeMaxX/2;
	  endX = startX+sizeMaxX-1;
	  sizeX = sizeMaxX;
	}
      else
	{
	  startX = 0;
	  endX = sizeOriX-1;
	  sizeX = sizeOriX;
	}
      if(sizeOriY>sizeMaxY)
	{
	  startY = sizeOriY/2-sizeMaxY/2;
	  endY = startY+sizeMaxY-1;
	  sizeY = sizeMaxY;
	}
      else
	{
	  startY = 0;
	  endY = sizeOriY-1;
	  sizeY = sizeOriY;
	}
      if(sizeOriZ>sizeMaxZ)
	{
	  startZ = sizeOriZ/2-sizeMaxZ/2;
	  endZ = startZ+sizeMaxZ-1;
	  sizeZ = sizeMaxZ;
	}
      else
	{
	  startZ = 0;
	  endZ = sizeOriZ-1;
	  sizeZ = sizeOriZ;
	}   

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
    }
  writer->Update();
  
  return EXIT_SUCCESS;
}
