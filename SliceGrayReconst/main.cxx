#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkMaximumImageFilter.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile radiusMin radiusMax direction " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension3 = 3;
  const unsigned int Dimension2 = 2;  
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageSliceIteratorWithIndex < ImageType > SliceIteratorType;
  typedef itk::Image< PixelType, Dimension2 > SliceImageType;
  typedef itk::ImageLinearIteratorWithIndex< SliceImageType > LinearIteratorType;
  typedef itk::BinaryBallStructuringElement< PixelType, Dimension2 > StructuringElementType;
  typedef itk::GrayscaleMorphologicalClosingImageFilter< SliceImageType, SliceImageType, StructuringElementType >  ClosingFilterType;
  typedef itk::GrayscaleErodeImageFilter< SliceImageType, SliceImageType, StructuringElementType> ErodeFilterType;
  typedef itk::MaximumImageFilter< SliceImageType, SliceImageType, SliceImageType> MaxFilterType;
  typedef itk::AbsoluteValueDifferenceImageFilter< SliceImageType, SliceImageType, SliceImageType> DiffFilterType;
  typedef itk::StatisticsImageFilter< SliceImageType > StatFilterType;
  typedef itk::SubtractImageFilter< SliceImageType, SliceImageType, SliceImageType> SubtractFilterType;
  typedef itk::ImageDuplicator< SliceImageType > DuplicateFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
  ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
  MaxFilterType::Pointer maxFilter = MaxFilterType::New();
  DiffFilterType::Pointer diffFilter = DiffFilterType::New();
  StatFilterType::Pointer statFilter = StatFilterType::New();
  SubtractFilterType::Pointer subtractFilter = SubtractFilterType::New();
  DuplicateFilterType::Pointer duplicateFilter = DuplicateFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  unsigned int radiusSEMin = atoi(argv[3]);
  unsigned int radiusSEMax = atoi(argv[4]);
  unsigned int sliceDirection = atoi(argv[5])-1;

  //Setting
  StructuringElementType  erodeSE;
  erodeSE.SetRadius( 1 );
  erodeSE.CreateStructuringElement();
  erodeFilter->SetKernel( erodeSE );
  unsigned int i, j;
  unsigned int direction[2];
  for(i=0, j=0; i<3; ++i)
    {
      if(i!=sliceDirection)
	{
	  direction[j] = i;
	  j++;
	}
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
  ImageType::DirectionType direction3 = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize( size );
  
  //Allocate output image
  ImageType::Pointer outputImage = ImageType::New();
  outputImage->SetRegions( region );
  outputImage->SetSpacing( spacing );
  outputImage->SetOrigin( origin );
  outputImage->SetDirection( direction3 );
  outputImage->Allocate();  

  //Slice specs
  SliceImageType::Pointer oneSliceIn = SliceImageType::New();
  SliceImageType::Pointer oneSliceOut = SliceImageType::New();
  SliceImageType::Pointer oneSlice = SliceImageType::New();
  SliceImageType::RegionType regionSlice;
  SliceImageType::RegionType::SizeType sizeSlice;
  SliceImageType::RegionType::IndexType indexSlice;

  //Allocate one slice
  ImageType::RegionType requestedRegion = reader->GetOutput()->GetRequestedRegion();
  indexSlice[ direction[0] ]    = requestedRegion.GetIndex()[ direction[0] ];
  indexSlice[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  sizeSlice[ direction[0] ]     = requestedRegion.GetSize()[  direction[0] ];
  sizeSlice[ 1- direction[0] ]  = requestedRegion.GetSize()[  direction[1] ];
  regionSlice.SetSize( sizeSlice );
  regionSlice.SetIndex( indexSlice );
  oneSliceIn->SetRegions( regionSlice );
  oneSliceIn->Allocate();
  oneSliceOut->SetRegions( regionSlice );
  oneSliceOut->Allocate();
  oneSliceOut->FillBuffer(0);

  //Iterator
  SliceIteratorType sliceItIn( reader->GetOutput(), reader->GetOutput()->GetRequestedRegion() );
  sliceItIn.SetFirstDirection( direction[1] );
  sliceItIn.SetSecondDirection( direction[0] );
  SliceIteratorType sliceItOut( outputImage, outputImage->GetRequestedRegion() );
  sliceItOut.SetFirstDirection( direction[1] );
  sliceItOut.SetSecondDirection( direction[0] );
  sliceItIn.GoToBegin();
  sliceItOut.GoToBegin();

  /*
  for(int a = 0;a<186;a++)
    {
      sliceItIn.NextSlice();
      sliceItOut.NextSlice();
    }
  */
  int a = 0;
  //Iterate through slices
  while( !sliceItIn.IsAtEnd() )
    //  if( !sliceItIn.IsAtEnd() )
    {
      std::cout<<"\r";
      std::cout<<int((a+1)*100/size[sliceDirection])<<"%"<<std::flush;
      //Extract slice
      LinearIteratorType lineItIn( oneSliceIn, oneSliceIn->GetRequestedRegion() );
      lineItIn.SetDirection( 1 - direction[0] );    
      lineItIn.GoToBegin(); 
      while ( !sliceItIn.IsAtEndOfSlice() )
	{
	  while ( !sliceItIn.IsAtEndOfLine() )
	    {
	      lineItIn.Set( sliceItIn.Get());
	      ++sliceItIn;
	      ++lineItIn;
	    }
	  lineItIn.NextLine();
	  sliceItIn.NextLine();
	}

      //Multi radius
      for(int radiusSE=radiusSEMin;radiusSE<=radiusSEMax;radiusSE++)
	{
	  StructuringElementType closingSE;
	  closingSE.SetRadius( radiusSE );
	  closingSE.CreateStructuringElement();
	  closingFilter->SetKernel( closingSE );

	  //Close the input slice
	  closingFilter->SetInput( oneSliceIn );
	  closingFilter->Update();
	  oneSlice = closingFilter->GetOutput();
	  //iterative grayscale reconstruction
	  int statSum = 1000;
	  while(statSum > 1)
	    { 
	      erodeFilter->SetInput( oneSlice );
	      maxFilter->SetInput1(erodeFilter->GetOutput());
	      maxFilter->SetInput2(oneSliceIn);
	      //test difference
	      diffFilter->SetInput1(oneSlice);
	      diffFilter->SetInput2(maxFilter->GetOutput());
	      statFilter->SetInput(diffFilter->GetOutput());
	      statFilter->Update();
	      statSum = statFilter->GetSum();
	      //assign for next iteration
	      oneSlice = maxFilter->GetOutput();
	      oneSlice->DisconnectPipeline();
	      //std::cout<<statSum<<std::endl;
	    }
	  //Final responseresult for current radius
	  subtractFilter->SetInput1(oneSlice);
	  subtractFilter->SetInput2(oneSliceIn);
	  subtractFilter->Update();
	  //record current max
	  duplicateFilter->SetInputImage(oneSliceOut);
	  duplicateFilter->Update();
	  maxFilter->SetInput1(subtractFilter->GetOutput());
	  maxFilter->SetInput2(duplicateFilter->GetOutput());
	  maxFilter->Update();
	  oneSliceOut = maxFilter->GetOutput();
	  oneSliceOut->DisconnectPipeline();
	}

      LinearIteratorType lineItOut( oneSliceOut, oneSliceOut->GetRequestedRegion() );
      lineItOut.SetDirection( 1 - direction[0] );
      lineItOut.GoToBegin();
      //Write slice
      while ( !sliceItOut.IsAtEndOfSlice() )
	{
	  while ( !sliceItOut.IsAtEndOfLine() )
	    {
	      sliceItOut.Set( lineItOut.Get());
	      ++sliceItOut;
	      ++lineItOut;
	    }
	  lineItOut.NextLine();
	  sliceItOut.NextLine();
	}

      //Prepare to extract next slice
      sliceItIn.NextSlice();
      sliceItOut.NextSlice();
      oneSliceIn->DisconnectPipeline();
      oneSliceOut->DisconnectPipeline();
      oneSliceIn->FillBuffer(0);
      oneSliceOut->FillBuffer(0);
      a++;
    }
  writer->SetInput(outputImage);
  std::cout<<std::endl;
  try
    {
      writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout << "ExceptionObject caught !" << std::endl;
      std::cout << err << std::endl;
      return -1;
    }

  return EXIT_SUCCESS;
}
