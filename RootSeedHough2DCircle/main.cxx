#include "itkHoughTransform2DCirclesImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"

#include <list>

#include "vnl/vnl_math.h"

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Missing Parameters " << std::endl;
      std::cerr << "Usage: " << argv[0] << " inputImage outputFile" << std::endl;
      return EXIT_FAILURE;
    }

  const unsigned int Dimension3 = 3;
  const unsigned int Dimension2 = 2;  
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageSliceIteratorWithIndex < ImageType > SliceIteratorType;
  typedef itk::Image< PixelType, Dimension2 > SliceImageType;
  typedef itk::ImageLinearIteratorWithIndex< SliceImageType > LinearIteratorType;
  typedef float AccumulatorPixelType;
  SliceImageType::IndexType localIndex;
  typedef itk::Image< AccumulatorPixelType, Dimension2 > AccumulatorImageType;


  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction3 = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize( size );

  //Slice specs
  SliceImageType::Pointer oneSlice = SliceImageType::New();
  SliceImageType::RegionType regionSlice;
  SliceImageType::RegionType::SizeType sizeSlice;
  SliceImageType::RegionType::IndexType indexSlice;
  ImageType::RegionType requestedRegion = reader->GetOutput()->GetRequestedRegion();
  indexSlice[ 0 ] = requestedRegion.GetIndex()[ 0 ];
  indexSlice[ 1 ] = requestedRegion.GetIndex()[ 1 ];
  sizeSlice[ 0 ]  = requestedRegion.GetSize()[ 0 ];
  sizeSlice[ 1 ]  = requestedRegion.GetSize()[ 1 ];
  regionSlice.SetSize( sizeSlice );
  regionSlice.SetIndex( indexSlice );
  oneSlice->SetRegions( regionSlice );
  oneSlice->Allocate();
  oneSlice->FillBuffer(0);
  //Iterator
  SliceIteratorType sliceIt( reader->GetOutput(), reader->GetOutput()->GetRequestedRegion() );
  sliceIt.SetFirstDirection( 0 );
  sliceIt.SetSecondDirection( 1 );
  sliceIt.GoToBegin();

  int upper = 20;
  int lower = 0;

  int a = 0;
  //skip first 0 slices
  while(a<lower)
    {
      sliceIt.NextSlice();
      a++;
    }


  int maxX=0, maxY=0, maxZ=0, maxAccu=0;

  //Iterate through slices
  while( !sliceIt.IsAtEnd() && (a<upper) )
    //  if( !sliceItIn.IsAtEnd() )
    {
      //Extract slice
      LinearIteratorType lineIt( oneSlice, oneSlice->GetRequestedRegion() );
      lineIt.SetDirection( 0 );    
      lineIt.GoToBegin(); 
      while ( !sliceIt.IsAtEndOfSlice() )
	{
	  while ( !sliceIt.IsAtEndOfLine() )
	    {
	      lineIt.Set( sliceIt.Get());
	      ++sliceIt;
	      ++lineIt;
	    }
	  lineIt.NextLine();
	  sliceIt.NextLine();
	}
      
      typedef itk::HoughTransform2DCirclesImageFilter<PixelType, AccumulatorPixelType> HoughTransformFilterType;
      HoughTransformFilterType::Pointer houghFilter = HoughTransformFilterType::New();
      houghFilter->SetInput( oneSlice );
      houghFilter->SetNumberOfCircles( 1 );
      houghFilter->SetMinimumRadius( 5 );
      houghFilter->SetMaximumRadius( 10 );
      houghFilter->Update();
      AccumulatorImageType::Pointer localAccumulator = houghFilter->GetOutput();
      AccumulatorImageType::IndexType accuIndex;
      HoughTransformFilterType::CirclesListType circles;
      circles = houghFilter->GetCircles( 1 );
      //std::cout << "Found " << circles.size() << " circle(s)." << std::endl;

      typedef HoughTransformFilterType::CirclesListType CirclesListType;
      CirclesListType::const_iterator itCircles = circles.begin();

      while( itCircles != circles.end() )
	{
	  //std::cout << "Center: ";
	  accuIndex[0] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
	  accuIndex[1] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
	  //std::cout << accuIndex[0]+1 << " "<< accuIndex[1]+1 << " "<< a+1 <<" "<< localAccumulator->GetPixel(accuIndex)<< std::endl;
	  if(maxAccu<localAccumulator->GetPixel(accuIndex))
	    {
	      maxAccu = localAccumulator->GetPixel(accuIndex);
	      maxX = accuIndex[0];
	      maxY = accuIndex[1];
	      maxZ = a;
	    }
	  //std::cout << "Radius: " << (*itCircles)->GetRadius()[0] << std::endl;
	  itCircles++;
	}

      //Prepare to extract next slice
      sliceIt.NextSlice();
      oneSlice->DisconnectPipeline();
      oneSlice->FillBuffer(0);
      a++;
    }

  FILE *seedFile = fopen(argv[2], "w");
  fprintf(seedFile, "%d %d %d 1\n", maxX, maxY, maxZ);
  fclose(seedFile);
  std::cout << "Seed: "<< maxX+1 << " "<< maxY+1 << " "<< maxZ+1 << std::endl;;
  return EXIT_SUCCESS;
}
