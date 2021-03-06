#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkMetaDataObject.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 
#include "gdcmGlobal.h"

//The program converts Dicom series to Analyze format with capability of multiple series read and get patient name as output  

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " DicomDirectory outputFileNameHead" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef double PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageSeriesReader< ImageType > ReaderType;
  typedef itk::GDCMImageIO ImageIOType;
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  nameGenerator->SetDirectory( argv[1] );
  char outFileName[50];

  //Pipeline
  reader->SetImageIO( dicomIO );
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  writer->SetInput( reader->GetOutput() );

  typedef std::vector< std::string >    SeriesIdContainer;
  const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
  SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
  SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

  std::cout << std::endl << "The directory: ";
  std::cout << std::endl << argv[1] << std::endl << std::endl;
  std::cout << "Contains the following DICOM Series: ";
  std::cout << std::endl ;
  while( seriesItr != seriesEnd )
    {
      std::cout << seriesItr->c_str() << std::endl;
      ++seriesItr;
    } 
     
  seriesItr = seriesUID.begin();

  while( seriesItr != seriesEnd )
    {
      std::strcpy(outFileName,argv[2]);
      std::string seriesIdentifier;
      seriesIdentifier = seriesItr->c_str();
      std::cout << std::endl;
      std::cout << "Now reading series: " << std::endl;
      std::cout << seriesIdentifier << std::endl;
      typedef std::vector< std::string > FileNamesContainer;
      FileNamesContainer fileNames;
      fileNames = nameGenerator->GetFileNames( seriesIdentifier ); 
      reader->SetFileNames( fileNames );
      reader->Update();
      typedef itk::MetaDataDictionary   DictionaryType;
      const  DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
      typedef itk::MetaDataObject< std::string > MetaDataStringType;
      DictionaryType::ConstIterator itr = dictionary.Begin();
      DictionaryType::ConstIterator end = dictionary.End();
      std::string entryId = "0010|0020";
      DictionaryType::ConstIterator tagItr = dictionary.Find( entryId );
      MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>( tagItr->second.GetPointer() );
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      std::cout << "Patient's Name is: " << tagvalue.c_str() << std::endl;
      
      //std::strcat(outFileName,"_");
      //std::strcat(outFileName,tagvalue.c_str());
      //std::strcat(outFileName,".hdr");

      std::string entryId2 = "0028|1053";
      DictionaryType::ConstIterator tagItr2 = dictionary.Find( entryId2 );
      MetaDataStringType::ConstPointer entryvalue2 = dynamic_cast<const MetaDataStringType *>( tagItr2->second.GetPointer() );
      std::cout << "Rescale slope is: " << std::setprecision (std::numeric_limits<double>::digits10 + 1) << entryvalue2->GetMetaDataObjectValue() << std::endl;

      std::string entryId3 = "0028|1052";
      DictionaryType::ConstIterator tagItr3 = dictionary.Find( entryId3 );
      MetaDataStringType::ConstPointer entryvalue3 = dynamic_cast<const MetaDataStringType *>( tagItr3->second.GetPointer() );
      std::cout << "Rescale intercept is: "<< std::setprecision (std::numeric_limits<double>::digits10 + 1) << entryvalue3->GetMetaDataObjectValue() << std::endl;     

      std::cout << "Writing the image as " << std::endl;
      std::cout << outFileName << std::endl<< std::endl;

      writer->SetFileName( outFileName );
      writer->Update();

      ++seriesItr;
    } 
     
  return EXIT_SUCCESS;
}
