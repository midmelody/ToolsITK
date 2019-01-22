#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage outputFileNameBase outputFileNameExt trim" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, 3 > ImageType3D;
  typedef itk::ImageFileReader< ImageType3D > ImageReaderType;
  typedef itk::Image< PixelType, 2 > ImageType2D; 
  typedef itk::ImageFileWriter< ImageType2D > ImageWriterType;
 
  //Read image
  ImageReaderType::Pointer Reader = ImageReaderType::New();
  Reader->SetFileName( argv[1] );
  Reader->Update();
  ImageWriterType::Pointer Writer = ImageWriterType::New();
 
  //Get image specs
  ImageType3D::SizeType size3D = Reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = size3D[0];
  int pCol = size3D[1];
  int pSli = size3D[2];

  //Trim parameter
  int trim = atoi(argv[4]);
  
  //Set output image
  ImageType2D::SizeType size2D;
  size2D[0] = pRow-trim;
  size2D[1] = pCol-trim;
  std::cout<<"Size: "<<size2D[0]<<" "<<size2D[1]<<std::endl;
  ImageType2D::RegionType region;
  region.SetSize( size2D ); 
  ImageType2D::Pointer imageOut = ImageType2D::New();
  imageOut->SetRegions( region );
  imageOut->Allocate();
  imageOut->FillBuffer( 0.0 );

  //Iterate to extract all slices
  ImageType3D::IndexType pixelInIndex;
  int i, j, k;
  unsigned char value;
  ImageType2D::IndexType pixelOutIndex;
  char filename[100];
  
  for(k=0;k<pSli;k++)
    {
      //Set value
      for(i=trim;i<pRow;i++)
	for(j=trim;j<pCol;j++)
	  {
	    pixelInIndex[0] = i;
	    pixelInIndex[1] = j;
	    pixelInIndex[2] = k;	  
	    value = Reader->GetOutput()->GetPixel(pixelInIndex);
	    pixelOutIndex[0] = i-trim;
	    pixelOutIndex[1] = j-trim;
	    imageOut->SetPixel(pixelOutIndex, value);
	    
	  }    
      //Output
      sprintf(filename, "%s_%d.%s", argv[2], k, argv[3]);
      //std::cout<<"Output: "<<filename<<std::endl;
      Writer->SetFileName( filename );
      Writer->SetInput(imageOut);
      Writer->Update();  
      imageOut->DisconnectPipeline();
    }


  
  return EXIT_SUCCESS;
}


