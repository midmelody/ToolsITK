#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkAndImageFilter.h" 
#include "itkInvertIntensityImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program 

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " oriImageFile dropImageFile outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::InvertIntensityImageFilter< ImageType, ImageType > InvFilterType;
  typedef itk::AndImageFilter< ImageType, ImageType, ImageType > AndFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer readerOri = ReaderType::New();
  ReaderType::Pointer readerDrp = ReaderType::New();  
  InvFilterType::Pointer invFilter = InvFilterType::New();
  AndFilterType::Pointer andFilter = AndFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  readerOri->SetFileName( argv[1] );
  readerDrp->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
 
  //Pipeline
  try
    {
      readerOri->Update();
      readerDrp->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  invFilter->SetInput( readerDrp->GetOutput() );
  invFilter->SetMaximum( 1 );
  andFilter->SetInput1( invFilter->GetOutput() );
  andFilter->SetInput2( readerOri->GetOutput() ); 


  //Write output  
  writer->SetInput( andFilter->GetOutput() );
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
