#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>        IteratorType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
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
  int pRow = size[0];
  int pCol = size[1];
  int pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  
  /*
  ImageType::IndexType pixelIndex;
  int i, j, k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  if(reader->GetOutput()->GetPixel(pixelIndex))
	    

	  inputImage(i,j,k).intensity = intReader->GetOutput()->GetPixel(pixelIndex);
	  inputImage(i,j,k).vesselness = vesReader->GetOutput()->GetPixel(pixelIndex); 
	  inputImage(i,j,k).scale = scaReader->GetOutput()->GetPixel(pixelIndex); 
	}
  */
  //Iterator
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  NeighborhoodIteratorType it( radius, reader->GetOutput(), reader->GetOutput()->GetRequestedRegion() );
  IteratorType out(image, region);
  IteratorType in(reader->GetOutput(), region);

  //Get sum
  for (it.GoToBegin(), in.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++in, ++out)
    {
      out.Set(0);
      if(in.Get() == 1)
	{
	  unsigned char sum = -1;
	  for (unsigned int i = 0; i < it.Size(); ++i)
	    {
	      if(it.GetPixel(i))
		sum++;
	    }
	  if(sum<0) sum=0;
	  out.Set(sum);
	}
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
