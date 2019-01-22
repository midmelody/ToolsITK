#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <cstdlib>
#include <ctime>
int main( int argc, char * argv[] )
{
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFileL inputImageFileH outputImageFolder sizeROI totalSample [mask]" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer readerL = ReaderType::New();
  ReaderType::Pointer readerH = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New();
  WriterType::Pointer writerL = WriterType::New();
  WriterType::Pointer writerH = WriterType::New();
  
  //Parameters
  readerL->SetFileName( argv[1] );
  readerL->Update();
  readerH->SetFileName( argv[2] );
  readerH->Update();  
  int sizeROI = atoi(argv[4]);
  int totalSample = atoi(argv[5]);
  if(argc==7)
    {
      readerM->SetFileName( argv[6] );
      readerM->Update();  
    }
    
  //Set image specs
  ImageType::SizeType  size;
  size[0] = sizeROI;
  size[1] = sizeROI;
  size[2] = sizeROI;
  ImageType::SizeType  sizeOri = readerL->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize( size ); 
  ImageType::Pointer imageL = ImageType::New();
  imageL->SetRegions( region );
  imageL->SetSpacing( readerL->GetOutput()->GetSpacing() );
  imageL->SetOrigin( readerL->GetOutput()->GetOrigin() );
  imageL->SetDirection( readerL->GetOutput()->GetDirection() );
  imageL->Allocate();
  imageL->FillBuffer( 0.0 );
  ImageType::Pointer imageH = ImageType::New();
  imageH->SetRegions( region );
  imageH->SetSpacing( readerL->GetOutput()->GetSpacing() );
  imageH->SetOrigin( readerL->GetOutput()->GetOrigin() );
  imageH->SetDirection( readerL->GetOutput()->GetDirection() );
  imageH->Allocate();
  imageH->FillBuffer( 0.0 );

  // If mask given, find the range
  ImageType::IndexType pixelIndex;
  int i, j, k, valueM;
  int minX, maxX, minY, maxY, minZ, maxZ;
  minX = 10000;
  minY = 10000;
  minZ = 10000;
  maxX = 0;
  maxY = 0;
  maxZ = 0;

  //std::cout<<sizeOri<<std::endl;
  
  if(argc==7)
    {
      for(i=0;i<sizeOri[0];i++)
	for(j=0;j<sizeOri[1];j++)
	  for(k=0;k<sizeOri[2];k++)
	    {
	      pixelIndex[0] = i;
	      pixelIndex[1] = j;
	      pixelIndex[2] = k;	  
	      valueM = readerM->GetOutput()->GetPixel(pixelIndex);
	      if(valueM)
		{
		  if(minX > i)
		    minX = i;
		  if(maxX < i)
		    maxX = i;
		  if(minY > j)
		    minY = j;
		  if(maxY < j)
		    maxY = j;
		  if(minZ > k)
		    minZ = k;
		  if(maxZ < k)
		    maxZ = k;		  
		}
	    }
      minX = minX + sizeROI;
      maxX = maxX - sizeROI;
      minY = minY + sizeROI;
      maxY = maxY - sizeROI;
      minZ = minZ + sizeROI;
      maxZ = maxZ - sizeROI;
      //std::cout<<"From ("<<minX<<" "<<minY<<" "<<minZ<<") To ("<<maxX<<" "<<maxY<<" "<<maxZ<<")"<<std::endl;
    }

 
  
  int ctROI;
  char fileNameL[100];
  char fileNameH[100];  
  srand((unsigned)time(0));
  ctROI=0;
  while(ctROI<totalSample)
    {
      //Get ROI
      ImageType::IndexType pixelInIndex, pixelOutIndex;
      int valueL, valueH;
      //Generate random center point within [sizeROI/2, sizeImage-sizeROI/2]
      //If given mask, generate within [min max] & check until mask=1
      int rangeX, rangeY, rangeZ, centerX, centerY, centerZ;
      if(argc==7)
	{
	  //valueM = 0;
	  //while(valueM!=1)
	  //{ 
	  centerX = minX + (rand() % static_cast<int>(maxX - minX + 1));
	  centerY = minY + (rand() % static_cast<int>(maxY - minY + 1));
	  centerZ = minZ + (rand() % static_cast<int>(maxZ - minZ + 1));
	  pixelInIndex[0] = centerX;
	  pixelInIndex[1] = centerY;
	  pixelInIndex[2] = centerZ;	  
	  valueM = readerM->GetOutput()->GetPixel(pixelInIndex);
	      //}
	  //std::cout<<ctROI+1<<": "<<valueM<<"; "<<centerX<<" "<<centerY<<" "<<centerZ<<std::endl;
	}
      else
	{
	  rangeX = sizeOri[0]-sizeROI/2;
	  rangeY = sizeOri[1]-sizeROI/2;
	  rangeZ = sizeOri[2]-sizeROI/2;
	  centerX = rand() % rangeX + sizeROI/2;
	  centerY = rand() % rangeY + sizeROI/2;
	  centerZ = rand() % rangeZ + sizeROI/2;
	  std::cout<<ctROI+1<<": "<<centerX<<" "<<centerY<<" "<<centerZ<<std::endl;	  
	}
      
      ctROI++;    
      for(i=0;i<sizeROI;i++)
	for(j=0;j<sizeROI;j++)
	  for(k=0;k<sizeROI;k++)
	    {
	      pixelInIndex[0] = i+centerX-sizeROI/2;
	      pixelInIndex[1] = j+centerY-sizeROI/2;
	      pixelInIndex[2] = k+centerZ-sizeROI/2;	  
	      valueL = readerL->GetOutput()->GetPixel(pixelInIndex);
	      valueH = readerH->GetOutput()->GetPixel(pixelInIndex);
	      //std::cout<<pixelInIndex[0]<<" "<<pixelInIndex[1]<<" "<<pixelInIndex[2]<<": "<<valueL<<" "<<valueH<<std::endl;
	      pixelOutIndex[0] = i;
	      pixelOutIndex[1] = j;
	      pixelOutIndex[2] = k;
	      imageL->SetPixel(pixelOutIndex, valueL);
	      imageH->SetPixel(pixelOutIndex, valueH);	      
	    }
      sprintf(fileNameL, "%s/Low/ROI_%04d.hdr", argv[3],ctROI);
      sprintf(fileNameH, "%s/High/ROI_%04d.hdr", argv[3],ctROI);     
      
      writerL->SetFileName( fileNameL );
      writerL->SetInput(imageL);
      writerL->Update();
      imageL->DisconnectPipeline();
      writerH->SetFileName( fileNameH );
      writerH->SetInput(imageH);     
      writerH->Update();
      imageH->DisconnectPipeline();
    }

  return EXIT_SUCCESS;
}
