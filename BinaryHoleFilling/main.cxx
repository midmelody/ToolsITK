#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile radius iterations" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::VotingBinaryIterativeHoleFillingImageFilter< ImageType >  HoleFillingFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  HoleFillingFilterType::Pointer holeFillingFilter = HoleFillingFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  ImageType::SizeType indexRadius;
  indexRadius[0] = atoi( argv[3] );
  indexRadius[1] = atoi( argv[3] );
  indexRadius[2] = atoi( argv[3] );
  holeFillingFilter->SetRadius( indexRadius );
  holeFillingFilter->SetBackgroundValue(   0 );
  holeFillingFilter->SetForegroundValue( 255 );
  holeFillingFilter->SetMaximumNumberOfIterations( atoi( argv[4] ) );

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
  
  holeFillingFilter->SetInput( reader->GetOutput() ); 

  //Write output  
  writer->SetInput( holeFillingFilter->GetOutput() );
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

  const unsigned int iterationsUsed = holeFillingFilter->GetCurrentNumberOfIterations();
  std::cout << "The filter used " << iterationsUsed << " iterations " << std::endl;

  const unsigned int numberOfPixelsChanged = holeFillingFilter->GetNumberOfPixelsChanged();
  std::cout << "and changed a total of " << numberOfPixelsChanged << " voxels " << std::endl;
  
  return EXIT_SUCCESS;
}
