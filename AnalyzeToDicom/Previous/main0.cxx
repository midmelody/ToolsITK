#include "itkImage.h"
#include "itkMinimumMaximumImageFilter.h" 
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "gdcmUIDGenerator.h"
#include <string>
#include <sstream>


//The program converts Analyze images to Dicom series according to a Dicom reference image
static void CopyDictionary (itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict);

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " InputImage ReferenceDicomDirectory OutputDicomDirectory" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef signed short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageSeriesReader< ImageType > RefReaderType;
  typedef itk::ImageFileReader< ImageType > InputReaderType;
  typedef itk::GDCMImageIO ImageIOType;
  typedef itk::GDCMSeriesFileNames InputNamesGeneratorType;
  typedef itk::NumericSeriesFileNames OutputNamesGeneratorType;
  typedef itk::ImageSeriesWriter< ImageType, ImageType > SeriesWriterType;

  //Read input and reference images
  InputReaderType::Pointer readerInput = InputReaderType::New();
  readerInput->SetFileName( argv[1] );
  readerInput->Update();
  ImageType::SizeType size = readerInput->GetOutput()->GetRequestedRegion().GetSize();

  ImageIOType::Pointer dicomIO = ImageIOType::New();
  InputNamesGeneratorType::Pointer refNames = InputNamesGeneratorType::New();
  refNames->SetInputDirectory( argv[2] );
  const RefReaderType::FileNamesContainer & filenames = refNames->GetInputFileNames();
  RefReaderType::Pointer readerRef = RefReaderType::New();
  readerRef->SetImageIO( dicomIO );
  readerRef->SetFileNames( filenames );
  readerRef->Update();

  //Copy the dictionary from the first image  
  RefReaderType::DictionaryRawPointer inputDict = (*(readerRef->GetMetaDataDictionaryArray()))[0];
  RefReaderType::DictionaryArrayType outputArray;
  gdcm::UIDGenerator suid;
  std::string seriesUID = suid.Generate();
  gdcm::UIDGenerator fuid;
  std::string frameOfReferenceUID = fuid.Generate();
  std::string studyUID;
  std::string sopClassUID;
  itk::ExposeMetaData<std::string>(*inputDict, "0020|000d", studyUID);
  itk::ExposeMetaData<std::string>(*inputDict, "0008|0016", sopClassUID);
  dicomIO->KeepOriginalUIDOn();
  
  //Write DICOM slices
  int i;
  for (i=0;i<size[2];i++)
    {
      // Create a new dictionary for this slice
      RefReaderType::DictionaryRawPointer dict = new RefReaderType::DictionaryType;
      // Copy the dictionary from the first slice
      CopyDictionary (*inputDict, *dict);
      // Set the UID's for the study, series, SOP  and frame of reference
      itk::EncapsulateMetaData<std::string>(*dict,"0020|000d", studyUID);
      itk::EncapsulateMetaData<std::string>(*dict,"0020|000e", seriesUID);
      itk::EncapsulateMetaData<std::string>(*dict,"0020|0052", frameOfReferenceUID);
      gdcm::UIDGenerator sopuid;
      std::string sopInstanceUID = sopuid.Generate();
      itk::EncapsulateMetaData<std::string>(*dict,"0008|0018", sopInstanceUID);
      itk::EncapsulateMetaData<std::string>(*dict,"0002|0003", sopInstanceUID);
      outputArray.push_back(dict);
    }
  
  // Make the output directory and generate the file names.
  itksys::SystemTools::MakeDirectory( argv[3] );
  
  // Generate the file names
  OutputNamesGeneratorType::Pointer outputNames = OutputNamesGeneratorType::New();
  std::string seriesFormat( argv[3] );
  seriesFormat = seriesFormat + "/" + "IM%d.dcm";
  outputNames->SetSeriesFormat (seriesFormat.c_str());
  outputNames->SetStartIndex (1);
  outputNames->SetEndIndex (size[2]);
  
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput( readerInput->GetOutput() );
  seriesWriter->SetImageIO( dicomIO );
  seriesWriter->SetFileNames( outputNames->GetFileNames() );
  seriesWriter->SetMetaDataDictionaryArray( &outputArray );
  seriesWriter->Update();
  return EXIT_SUCCESS;
}
 
void CopyDictionary (itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict)
{
  typedef itk::MetaDataDictionary DictionaryType;
 
  DictionaryType::ConstIterator itr = fromDict.Begin();
  DictionaryType::ConstIterator end = fromDict.End();
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
 
  while( itr != end )
    {
    itk::MetaDataObjectBase::Pointer  entry = itr->second;
 
    MetaDataStringType::Pointer entryvalue =
      dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;
    if( entryvalue )
      {
      std::string tagkey   = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      itk::EncapsulateMetaData<std::string>(toDict, tagkey, tagvalue);
      }
    ++itr;
    }
}
