// The program reads a 3D image from a non DICOM file and writes it as a series of DICOM slices
// Output header comes from a Dicom reference image

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "gdcmUIDGenerator.h"
#include "itkImageSeriesReader.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"
#include <vector>
#include "itksys/SystemTools.hxx"

static void CopyDictionary (itk::MetaDataDictionary &fromDict,
                     itk::MetaDataDictionary &toDict);

int main( int argc, char* argv[] )
{   
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " InputImage OutputDicomDirectory ReferenceDicomDirectory" << std::endl;
      return EXIT_FAILURE;
    }

  typedef int    PixelType;
  const unsigned int      Dimension = 3;
  typedef itk::Image< PixelType, Dimension >      ImageType;

  //Image input
  typedef itk::ImageFileReader< ImageType >       ReaderInType;
  ReaderInType::Pointer readerIn = ReaderInType::New();
  readerIn->SetFileName( argv[1] );
  readerIn->Update();
  ImageType::RegionType region = readerIn->GetOutput()->GetLargestPossibleRegion();
  ImageType::IndexType start = region.GetIndex();
  ImageType::SizeType  size  = region.GetSize();
  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer gdcmIO = ImageIOType::New();

  //Reference input
  typedef itk::ImageSeriesReader< ImageType >     ReaderRefType;
  typedef itk::GDCMSeriesFileNames                InNamesGeneratorType;
  InNamesGeneratorType::Pointer namesGeneratorIn = InNamesGeneratorType::New();
  namesGeneratorIn->SetInputDirectory( argv[3] );
  const ReaderRefType::FileNamesContainer & filenames = namesGeneratorIn->GetInputFileNames();
  const unsigned int numberOfFileNames =  filenames.size();
  if(size[2] != numberOfFileNames)
    {
      std::cerr << "Input and reference images do not match!" << std::endl;
      return EXIT_FAILURE;
    }
  else
    std::cout << "Total slice number: "<< numberOfFileNames << std::endl;
  ReaderRefType::Pointer readerRef = ReaderRefType::New();
  readerRef->SetImageIO( gdcmIO );
  readerRef->SetFileNames( filenames );
  readerRef->Update();
  
  //Write output
  typedef itk::NumericSeriesFileNames             OutNamesGeneratorType;
  const char * outputDirectory = argv[2];
  itksys::SystemTools::MakeDirectory( outputDirectory );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;
  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;
  typedef itk::ImageSeriesWriter< ImageType, Image2DType >  SeriesWriterType;
  OutNamesGeneratorType::Pointer namesGeneratorOut = OutNamesGeneratorType::New();

  //Copy image dictionary from reference 
  ReaderRefType::DictionaryRawPointer inputDict = (*(readerRef->GetMetaDataDictionaryArray()))[0];
  ReaderRefType::DictionaryArrayType outputArray;
  gdcm::UIDGenerator suid;
  std::string seriesUID = suid.Generate();
  gdcm::UIDGenerator fuid;
  //std::string frameOfReferenceUID = fuid.Generate();
  std::string studyUID;
  std::string sopClassUID;
  itk::ExposeMetaData<std::string>(*inputDict, "0020|000d", studyUID);
  itk::ExposeMetaData<std::string>(*inputDict, "0008|0016", sopClassUID);
  gdcmIO->KeepOriginalUIDOn();

  for (unsigned int f = 0; f < size[2]; f++)
    {
      // Create a new dictionary for this slice
      ReaderRefType::DictionaryRawPointer dict = new ReaderRefType::DictionaryType;
      // Copy the dictionary from the first slice
      CopyDictionary (*inputDict, *dict);
      // Set the UID's for the study, series, SOP  and frame of reference
      itk::EncapsulateMetaData<std::string>(*dict,"0020|000d", studyUID);
      itk::EncapsulateMetaData<std::string>(*dict,"0020|000e", seriesUID);
      //itk::EncapsulateMetaData<std::string>(*dict,"0020|0052", frameOfReferenceUID);
      gdcm::UIDGenerator sopuid;
      std::string sopInstanceUID = sopuid.Generate();
      itk::EncapsulateMetaData<std::string>(*dict,"0008|0018", sopInstanceUID);
      itk::EncapsulateMetaData<std::string>(*dict,"0002|0003", sopInstanceUID);
      itk::EncapsulateMetaData<std::string>(*dict,"0008|0008", "DERIVED\\SECONDARY");
      itk::EncapsulateMetaData<std::string>(*dict,"0008|103e", "SEGMENTATION");
      std::ostringstream value;
      value.str("");
      value << f + 1;
      itk::EncapsulateMetaData<std::string>(*dict,"0020|0013", value.str());
      value.str("");
      value << 1001;
      itk::EncapsulateMetaData<std::string>(*dict,"0020|0011", value.str());
      ImageType::PointType position;
      ImageType::IndexType index;
      index[0] = 0;
      index[1] = 0;
      index[2] = f;
      readerRef->GetOutput()->TransformIndexToPhysicalPoint(index, position);
      value.str("");
      value << position[0] << "\\" << position[1] << "\\" << position[2];
      itk::EncapsulateMetaData<std::string>(*dict,"0020|0032", value.str());
      value.str("");
      value << position[2];
      itk::EncapsulateMetaData<std::string>(*dict,"0020|1041", value.str());
      value.str("");
      value << "0";
      itk::EncapsulateMetaData<std::string>(*dict,"0028|1052", value.str());
      value.str("");
      value << "1";
      itk::EncapsulateMetaData<std::string>(*dict,"0028|1053", value.str());
      outputArray.push_back(dict);
    }

  //Write output
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput( readerIn->GetOutput() );
  seriesWriter->SetImageIO( gdcmIO );
  std::string format = outputDirectory;
  format += "/image%03d.dcm";
  namesGeneratorOut->SetSeriesFormat( format.c_str() );
  namesGeneratorOut->SetStartIndex( start[2] );
  namesGeneratorOut->SetEndIndex( start[2] + size[2] - 1 );
  namesGeneratorOut->SetIncrementIndex( 1 );
  seriesWriter->SetFileNames( namesGeneratorOut->GetFileNames() );
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
