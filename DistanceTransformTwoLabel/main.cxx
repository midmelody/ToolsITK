#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " maskImage inputLabel1 inputLabel2 outputLabel outputDist" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::SignedMaurerDistanceMapImageFilter < ImageType, ImageType > DistanceTransformFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer readerM = ReaderType::New();
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  DistanceTransformFilterType::Pointer DTFilter1 = DistanceTransformFilterType::New();
  DistanceTransformFilterType::Pointer DTFilter2 = DistanceTransformFilterType::New();
  WriterType::Pointer writerL = WriterType::New();
  WriterType::Pointer writerD = WriterType::New();

  //Parameters
  readerM->SetFileName( argv[1] );
  reader1->SetFileName( argv[2] );
  reader2->SetFileName( argv[3] );
  writerL->SetFileName( argv[4] );
  writerD->SetFileName( argv[5] );

  //Pipeline
  try
    {
      readerM->Update();
      reader1->Update();
      reader2->Update(); 
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  //Get image specs
  ImageType::SpacingType spacing = readerM->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = readerM->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = readerM->GetOutput()->GetDirection();
  ImageType::SizeType  size = readerM->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = size[0];
  int pCol = size[1];
  int pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  ImageType::Pointer imageL = ImageType::New();
  imageL->SetRegions( region );
  imageL->SetSpacing( spacing );
  imageL->SetOrigin( origin );
  imageL->SetDirection( direction );
  imageL->Allocate();  
  ImageType::Pointer imageD = ImageType::New();
  imageD->SetRegions( region );
  imageD->SetSpacing( spacing );
  imageD->SetOrigin( origin );
  imageD->SetDirection( direction );
  imageD->Allocate();  

  //Assign label and dist infor
  DTFilter1->SetInput(reader1->GetOutput());
  DTFilter1->SetBackgroundValue( 0 );
  DTFilter1->UseImageSpacingOn();

  DTFilter2->SetInput(reader2->GetOutput());
  DTFilter2->SetBackgroundValue( 0 );
  DTFilter2->UseImageSpacingOn();

  DTFilter1->Update(); 
  DTFilter2->Update(); 

  ImageType::IndexType pixelIndex;
  int i, j, k;
  float DT1, DT2;

  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  if(readerM->GetOutput()->GetPixel(pixelIndex))
	    {   
	      DT1 = DTFilter1->GetOutput()->GetPixel(pixelIndex);
	      DT2 = DTFilter2->GetOutput()->GetPixel(pixelIndex);
	      if(DT1<DT2)
		{
		  imageL->SetPixel(pixelIndex, 1);
		  imageD->SetPixel(pixelIndex, DT1);
		}
	      else if(DT1>DT2)
		{
		  imageL->SetPixel(pixelIndex, 2);
		  imageD->SetPixel(pixelIndex, DT2);
		}
	    }
	  else
	    {
	      imageL->SetPixel(pixelIndex, 0);
	      imageD->SetPixel(pixelIndex, 0);
	    }
	}

  //Write output  
  writerL->SetInput( imageL );
  writerD->SetInput( imageD );
  try
    {
      writerL->Update();
      writerD->Update();    
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
