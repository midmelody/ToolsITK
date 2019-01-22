#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkConstrainedValueDifferenceImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"

//The program enhance the airway tubular structure with grayscale morphological operations as mentioned in "Segmentation and analysis of the human airway tree from three-demensional X-ray CT images"
//Black top hat: returns an image, containing the "objects" or "elements" that:
//Are "smaller" than the structuring element, and
//Are darker than their surroundings.
//Grayscale closing and take difference by original   

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile radius" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::BinaryBallStructuringElement< PixelType, Dimension > StructuringElementType;
  typedef itk::GrayscaleMorphologicalClosingImageFilter< ImageType, ImageType, StructuringElementType >  ClosingFilterType;
  typedef itk::ConstrainedValueDifferenceImageFilter< ImageType, ImageType, ImageType > SubtractionFilterType;
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType> RescaleFilterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  SubtractionFilterType::Pointer topHat = SubtractionFilterType::New();
  ClosingFilterType::Pointer closing = ClosingFilterType::New();
  
  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  StructuringElementType  structuringElement;
  structuringElement.SetRadius( atoi(argv[3]) ); //(argv[3]+1)x(argv[3]+1) structuring element
  structuringElement.CreateStructuringElement();
  
  //Pipeline
  //Read image
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
  //Setup closing
  closing->SetInput( reader->GetOutput() );
  closing->SetKernel( structuringElement );
  //Top hat
  topHat->SetInput1( closing->GetOutput() );
  topHat->SetInput2( reader->GetOutput() );
  //Rescale
  rescaleFilter->SetInput( topHat->GetOutput() );
  rescaleFilter->SetOutputMinimum( 0 );
  rescaleFilter->SetOutputMaximum( 255 );
  //Write Image
  writer->SetInput( rescaleFilter->GetOutput() );
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
