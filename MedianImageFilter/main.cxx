#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMedianImageFilter.h"
#include "itkImageFileWriter.h"
 
int main(int argc, char * argv[])
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " InputImageFile OutputImageFile [radius]" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::MedianImageFilter<ImageType, ImageType > FilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
 
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  FilterType::Pointer medianFilter = FilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  FilterType::InputSizeType radius;
  radius.Fill(2);
  if (argc > 3)
    {
      radius.Fill(atoi(argv[3]));
    }
  medianFilter->SetRadius(radius);
  medianFilter->SetInput( reader->GetOutput() );
  writer->SetFileName( argv[2] );

  writer->SetInput( medianFilter->GetOutput() );
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
