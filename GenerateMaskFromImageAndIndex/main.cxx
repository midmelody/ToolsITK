#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile inputIndexFile outputImageFile mode(0:xyz, 1:yxz)" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[3] );
  int mode = atoi( argv[4] );
  if(mode==0)
    {
      std::cout<<"xyz mode"<<std::endl;
    }
  else if(mode==1)
    {
      std::cout<<"yxz mode"<<std::endl;
    }
  else
    {
      std::cerr << "invalid mode" << std::endl;
      return EXIT_FAILURE;
    }
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
  ImageType::RegionType region;
  region.SetSize( size );
  
  //Allocate new image
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  
  image->FillBuffer( itk::NumericTraits< PixelType >::Zero );
  
  //Read mask indecies and set to 1
  ImageType::IndexType pixelIndex;
  int x, y, z;
  std::cout<<"Reading mask file.......";
  FILE *maskFile = fopen(argv[2], "r");
  if(maskFile == NULL) 
    {
      std::cerr << " Seed file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(maskFile, "%d %d %d", &x, &y, &z) == 3)  
    {
      if(mode==0)
        {
          pixelIndex[0]=x;
          pixelIndex[1]=y;
          pixelIndex[2]=z+1;
        }
      else if(mode==1)
        {
          pixelIndex[0]=y;
          pixelIndex[1]=x;
          pixelIndex[2]=z+1;
        }
 
	  image->SetPixel(pixelIndex, 1);
	}
  fclose(maskFile);
  std::cout<<"finished"<<std::endl;
  
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
