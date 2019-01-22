#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageFileWriter.h"
#include <iostream>

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " LabelImage1 LabelImage2 maskImage originalImage outputImage " << std::endl;
      return EXIT_FAILURE;
    }
 
  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef unsigned char OutPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< OutPixelType, Dimension > OutImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< OutImageType > WriterType;
    
  //Filters
  ReaderType::Pointer readerL1 = ReaderType::New();
  ReaderType::Pointer readerL2 = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New();
  ReaderType::Pointer readerO = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  WriterType::Pointer writerF = WriterType::New();   

  //Parameters
  readerL1->SetFileName( argv[1] );
  readerL2->SetFileName( argv[2] );
  readerM->SetFileName( argv[3] );
  readerO->SetFileName( argv[4] );  
  writer->SetFileName( argv[5] );
  writerF->SetFileName( "tempCC.hdr");
  //Pipeline
  try
    {
      readerL1->Update();
      readerL2->Update();
      readerM->Update();
      readerO->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  //Get image specs
  ImageType::SpacingType spacing = readerL1->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = readerL1->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = readerL1->GetOutput()->GetDirection();
  ImageType::SizeType  size = readerL1->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  OutImageType::Pointer image = OutImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  
  OutImageType::Pointer imageF = OutImageType::New();
  imageF->SetRegions( region );
  imageF->SetSpacing( spacing );
  imageF->SetOrigin( origin );
  imageF->SetDirection( direction );
  imageF->Allocate();
 
  //Read image and set absolute part and fuzzy part
  ImageType::IndexType pixelIndex;
  int i, j, k, valueL1, valueL2, valueM, valueO;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  valueL1 = readerL1->GetOutput()->GetPixel(pixelIndex);
	  valueL2 = readerL2->GetOutput()->GetPixel(pixelIndex);
	  valueM = readerM->GetOutput()->GetPixel(pixelIndex);
	  valueO = readerO->GetOutput()->GetPixel(pixelIndex);

	  if(valueL1>valueL2)
	    image->SetPixel(pixelIndex, 2);
	  else
	    image->SetPixel(pixelIndex, 1);
	  
	  //if((i==183)&&(j==264)&&(k==243))
	  //std::cout<<valueL1<<" "<<valueL2<<" "<<valueM<<" "<<valueO<<std::endl;
	  if((fabs(valueL1-valueL2)<20)&&valueM&&(valueO<-700))
	    imageF->SetPixel(pixelIndex, 0);
	  else
	    imageF->SetPixel(pixelIndex, 0);
	}	  
  
  typedef itk::ConnectedComponentImageFilter< OutImageType, OutImageType > ConnectedComponentFilterType;
  typedef itk::RelabelComponentImageFilter< OutImageType, OutImageType >  RelabelFilterType;
  ConnectedComponentFilterType::Pointer connectComponentFilter = ConnectedComponentFilterType::New();
  RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  connectComponentFilter->SetInput( imageF );
  connectComponentFilter->SetFullyConnected(0);
  relabelFilter->SetInput( connectComponentFilter->GetOutput() );
  relabelFilter->SetMinimumObjectSize(1000);
  relabelFilter->Update();

  typedef std::vector< itk::SizeValueType > SizesInPixelsType;
  const SizesInPixelsType &  sizesInPixels = relabelFilter->GetSizeOfObjectsInPixels();
  SizesInPixelsType::const_iterator sizeItr = sizesInPixels.begin();
  SizesInPixelsType::const_iterator sizeEnd = sizesInPixels.end();
  //std::cout << "Number of pixels per class " << std::endl;
  unsigned int kclass = 1;
  while( sizeItr != sizeEnd )
    {
      //std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
      ++kclass;
      ++sizeItr;
    }
  writerF->SetInput( relabelFilter->GetOutput() );
  writerF->Update();

  kclass = kclass-1;
  int Accu[kclass];
  int Count[kclass];
  int ct=0;
  sizeItr = sizesInPixels.begin();
  sizeEnd = sizesInPixels.end();
  while( sizeItr != sizeEnd )
    {
      Count[ct] = *sizeItr;
      Accu[ct] = 0;
      ++ct;
      ++sizeItr;
    }
  
  int valueSep, valueCC;
  int ct4 = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  valueSep = image->GetPixel(pixelIndex);
	  valueCC = relabelFilter->GetOutput()->GetPixel(pixelIndex);
	  if(valueCC!=0)
	    Accu[valueCC-1] = Accu[valueCC-1] + valueSep;
	}

  for(ct=0;ct<kclass;ct++)
    {
      if(float(Accu[ct])/float(Count[ct])>1.5)
	Accu[ct] = 2;
      else
	Accu[ct] = 1;
    }
  
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  valueCC = relabelFilter->GetOutput()->GetPixel(pixelIndex);
	  if(valueCC!=0)
	    image->SetPixel(pixelIndex, Accu[valueCC-1]);
	}

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
