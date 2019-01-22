#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>


//ITK settings
const unsigned int Dimension = 3;
typedef unsigned char PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::Image< PixelType, Dimension+1 > MaskType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::BinaryThresholdImageFilter <ImageType, ImageType> BinThresholdImageFilterType;
typedef itk::ImageFileWriter< MaskType > WriterType;
typedef itk::CastImageFilter<ImageType, MaskType> CastType;
typedef itk::PasteImageFilter<MaskType, MaskType> PasteType;

// Store a 3D image inside a 4D image at a particular index
MaskType::Pointer copy_3d_to_4d(ImageType::Pointer image, MaskType::Pointer image4d, int idx)
{
    //this works by first casting the 3d image to a 4d one (with a single "frame")
    //then "pasting" this 4d image inside the image4d at the correct place


    CastType::Pointer castFilter = CastType::New();
    castFilter->SetInput(image);
    castFilter->Update();
    PasteType::Pointer pasteFilter = PasteType::New();
    pasteFilter->SetSourceImage(castFilter->GetOutput());
    pasteFilter->SetDestinationImage(image4d);
    pasteFilter->SetSourceRegion(castFilter->GetOutput()->GetLargestPossibleRegion());
    MaskType::IndexType destinationIndex;
    destinationIndex[0] = 0;
    destinationIndex[1] = 0;
    destinationIndex[2] = 0;
    destinationIndex[3] = idx;
    pasteFilter->SetDestinationIndex(destinationIndex);
    pasteFilter->Update();
    MaskType::Pointer newimage4d = pasteFilter->GetOutput();
    return newimage4d;
}

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile totalObject outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  BinThresholdImageFilterType::Pointer thresholdFilter = BinThresholdImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  int totalObj = atoi( argv[2] );
  writer->SetFileName( argv[3] );
 
  //Get image specs
  reader->Update();
  ImageType::SpacingType spacing0 = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin0 = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction0 = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::Pointer image = ImageType::New();
    
  //Set 4D image
  MaskType::Pointer mask  = MaskType::New();
  MaskType::IndexType start;
  for (int i=0; i<4; ++i) start[i] = 0;
  MaskType::SizeType masksize;
  for (int i=0; i<3; ++i) masksize[i] = size[i];
  masksize[3] = totalObj + 1; // background + objects
  MaskType::RegionType region;
  region.SetSize( masksize );
  region.SetIndex( start );
  mask->SetRegions( region );
  MaskType::DirectionType direction;
  for (int i=0; i<3; ++i) 
    for (int j=0; j<3; ++j)
      {
	direction[i][j]=direction0[i][j];
	direction[3][i]=0;
	direction[i][3]=0;
      }
  direction[3][3]=1;
  mask->SetDirection(direction);
  MaskType::SpacingType spacing;
  for (int i=0; i<3; ++i) spacing[i]=spacing0[i];
  spacing[3]=1;
  mask->SetSpacing(spacing);
  MaskType::PointType origin;
  for (int i=0; i<3; ++i) origin[i]=origin0[i];
  origin[3]=0;
  mask->SetOrigin(origin);
  mask->Allocate(); 

  for(int i=0; i<masksize[3]; i++)
    {
      thresholdFilter->SetInput( reader->GetOutput() );
      thresholdFilter->SetUpperThreshold( i );
      thresholdFilter->SetLowerThreshold( i );
      thresholdFilter->SetInsideValue( 1 );
      thresholdFilter->SetOutsideValue( 0 );
      thresholdFilter->Update();
      image = thresholdFilter->GetOutput();
      // copy to mask
      mask=copy_3d_to_4d(image, mask.GetPointer(), i);
      image->DisconnectPipeline();
    }

  
  //Write output  
  writer->SetInput( mask );
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
