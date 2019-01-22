#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkImage.h"

#include "itkAnalyzeObjectLabelMapImageIO.h"
#include "itkAnalyzeObjectMap.h"
#include "itkAnalyzeObjectLabelMapImageIOFactory.h"

int main( int argc, char ** argv )
{
  if ( argc < 3 )
    {
      std::cerr << "USAGE: " << argv[0] << "inputObjectFileName outputImageFileName" << std::endl;
    }
  const char *InputObjectFileName = argv[1];
  const char *OuptputObjectFileName = argv[2];

  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
 
  itk::ObjectFactoryBase::RegisterFactory( itk::AnalyzeObjectLabelMapImageIOFactory::New() );

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( InputObjectFileName);
  try
    {
      reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl
		<< err << std::endl;
      return EXIT_FAILURE;
    }

  writer->SetFileName(OuptputObjectFileName);
  writer->SetInput(reader->GetOutput());
  try
    {
      writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl
		<< err << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
