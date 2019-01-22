#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile inputMaskFile inputROIImageFile inputROIMaskFile outputImageFile outputMaskFile" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer readerIMG = ReaderType::New();
  ReaderType::Pointer readerMSK = ReaderType::New();
  ReaderType::Pointer readerROI = ReaderType::New();
  ReaderType::Pointer readerMSR = ReaderType::New();  
  WriterType::Pointer writerIMG = WriterType::New(); 
  WriterType::Pointer writerMSK = WriterType::New();
  
  //Parameters
  readerIMG->SetFileName( argv[1] );
  readerIMG->Update();
  readerMSK->SetFileName( argv[2] );
  readerMSK->Update();
  readerROI->SetFileName( argv[3] );
  readerROI->Update();
  readerMSR->SetFileName( argv[4] );
  readerMSR->Update();  
  writerIMG->SetFileName( argv[5] );
  writerMSK->SetFileName( argv[6] );
  
  //Parse for center location
  char *filenameRef;
  int centerX, centerY, centerZ;
  int caseID, nodID, rad;
  char *lastdot;
  filenameRef = basename(argv[4]);
  lastdot = strrchr(filenameRef, '.');
  if (lastdot != NULL)
    *lastdot = '\0'; 
  sscanf(filenameRef, "%d_%d_%d_%d_%d_%d", &caseID, &nodID, &rad, &centerX, &centerY, &centerZ);
  std::cout<<"Center: "<<centerX<<" "<<centerY<<" "<<centerZ<<std::endl;

  //Get image specs
  ImageType::SizeType sizeIMG = readerIMG->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = sizeIMG[0];
  int pCol = sizeIMG[1];
  int pSli = sizeIMG[2];
  ImageType::SizeType sizeROI = readerROI->GetOutput()->GetRequestedRegion().GetSize();
  
  //Set output
  ImageType::SizeType sizeOut;
  sizeOut[0] = pRow;
  sizeOut[1] = pCol;  
  sizeOut[2] = sizeROI[2];
  std::cout<<"Size: "<<sizeOut<<std::endl; 
  ImageType::RegionType region;
  region.SetSize( sizeOut );
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( readerIMG->GetOutput()->GetSpacing() );
  image->SetOrigin( readerIMG->GetOutput()->GetOrigin() );
  image->SetDirection( readerIMG->GetOutput()->GetDirection() );
  image->Allocate();
  image->FillBuffer( -1000.0 );
  ImageType::Pointer mask = ImageType::New();
  mask->SetRegions( region );
  mask->SetSpacing( readerIMG->GetOutput()->GetSpacing() );
  mask->SetOrigin( readerIMG->GetOutput()->GetOrigin() );
  mask->SetDirection( readerIMG->GetOutput()->GetDirection() );
  mask->Allocate();
  mask->FillBuffer( 0.0 );
  
  //Copy intensity from original image and replace around the center point
  int minX, maxX, minY, maxY, minZ, maxZ;
  minX = centerX - sizeROI[0]/2;
  minY = centerY - sizeROI[1]/2;
  minZ = centerZ - sizeROI[2]/2;
  if(minX < 0)
    minX = 0;
  if(minY < 0)
    minY = 0;
  if(minZ < 0)
    minZ = 0;
  maxX = minX + sizeROI[0];
  maxY = minY + sizeROI[1];
  maxZ = minZ + sizeROI[2];
  if(maxX > pRow)
    maxX = pRow;
  if(maxY > pCol)
    maxY = pCol;
  if(maxZ > pSli)
    maxZ = pSli;
 
  //Set output image value around selected center
  ImageType::IndexType pixelIndexIMG, pixelIndexROI, pixelIndexOut;
  int i, j, k, valueIMG, valueROI, valueMSK, valueMSR;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=minZ;k<maxZ;k++)
	{
	  pixelIndexIMG[0] = i;
	  pixelIndexIMG[1] = j;
	  pixelIndexIMG[2] = k;
	  valueIMG = readerIMG->GetOutput()->GetPixel(pixelIndexIMG);
	  valueMSK = readerMSK->GetOutput()->GetPixel(pixelIndexIMG);
	  
	  pixelIndexOut[0] = i;
	  pixelIndexOut[1] = j;
	  pixelIndexOut[2] = k-minZ;
	  mask->SetPixel(pixelIndexOut, valueMSK);
	  image->SetPixel(pixelIndexOut, valueIMG);

	  if((i>=minX)&&(i<maxX-1)&&(j>=minY)&&(j<maxY-1)&&(k<maxZ-1))
	    {
	      pixelIndexROI[0] = i-minX;
	      pixelIndexROI[1] = j-minY;	  
	      pixelIndexROI[2] = k-minZ;
	      valueMSR = readerMSR->GetOutput()->GetPixel(pixelIndexROI);
	      valueROI = readerROI->GetOutput()->GetPixel(pixelIndexROI);
	      if(valueMSR)
		image->SetPixel(pixelIndexOut, valueROI);	      
	    }  
	}
  
  writerIMG->SetInput( image );
  writerIMG->Update();
  writerMSK->SetInput( mask );
  writerMSK->Update();

  return EXIT_SUCCESS;
}
