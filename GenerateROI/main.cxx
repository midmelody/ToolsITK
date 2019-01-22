#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile direction(1:x; 2:y; 3:z) location(1: Left/Upper; 2:Right/Lower; 3:Central) widthRatio(max 0.5)" << std::endl;
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
  int dir = atoi(argv[3]);
  int loc = atoi(argv[4]);
  float wid = atof(argv[5]);
  if(wid>0.5)
    {
      std::cerr << "Invalid Width Value"<<std::endl;
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
 
  //Set image specs
  ImageType::RegionType region;
  region.SetSize( size ); 
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();
  image->FillBuffer( 0.0 );

  //Set values
  ImageType::IndexType pixelIndex;
  int i, j, k, value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0] = i;
	  pixelIndex[1] = j;
	  pixelIndex[2] = k;
	  if(dir==1)
	    {
	      if(loc==1)
		{
		  if(i<pRow*wid)
		    value = 1;
		  else
		    value = 0;
		}
	      else if(loc==2)
		{
		  if(i<pRow*(1-wid))
		    value = 0;
		  else
		    value = 1;
		}	
	      else if(loc==3)
		{
		  if((i>pRow*(0.5-wid/2))&&(i<pRow*(0.5+wid/2)))
		    value = 1;
		  else
		    value = 0;
		}
	      else
		{
		  std::cerr << "Invalid Location Value"<<std::endl;
		  return EXIT_FAILURE;
		}
	    }
	  else if(dir==2)	
	    {
	      if(loc==1)
		{
		  if(j<pCol*wid)
		    value = 1;
		  else
		    value = 0;
		}
	      else if(loc==2)
		{
		  if(j<pCol*(1-wid))
		    value = 0;
		  else
		    value = 1;
		}
	      else if(loc==3)
		{
		  if((j>pCol*(0.5-wid/2))&&(j<pCol*(0.5+wid/2)))
		    value = 1;
		  else
		    value = 0;
		}	
	      else
		{
		  std::cerr << "Invalid Location Value"<<std::endl;
		  return EXIT_FAILURE;
		}
	    }
	  else if(dir==3)
	    {
	      if(loc==1)
		{
		  if(k<pSli*wid)
		    value = 1;
		  else
		    value = 0;
		}
	      else if(loc==2)
		{
		  if(k<pSli*(1-wid))
		    value = 0;
		  else
		    value = 1;
		}
	      else if(loc==3)
		{
		  if((k>pSli*(0.5-wid/2))&&(k<pSli*(0.5+wid/2)))
		    value = 1;
		  else
		    value = 0;
		}	
	      else
		{
		  std::cerr << "Invalid Location Value"<<std::endl;
		  return EXIT_FAILURE;
		}
	    }
	  else
	    {
	      std::cerr << "Invalid Direction Value"<<std::endl;
	      return EXIT_FAILURE;
	    }
	  image->SetPixel( pixelIndex, value );
	}

  writer->SetInput( image );
  writer->Update();
  
  return EXIT_SUCCESS;
}
