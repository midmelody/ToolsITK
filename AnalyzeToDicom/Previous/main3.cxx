// The program reads a 3D image from a non DICOM file and writes it as a series of DICOM slices
// Output header comes from a Dicom reference image

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"
#include <vector>
#include "itksys/SystemTools.hxx"

int main( int argc, char* argv[] )
{   
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " InputImage OutputDicomDirectory [ReferenceDicomDirectory]" << std::endl;
      return EXIT_FAILURE;
    }

  typedef signed short    PixelType;
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
  ImageIOType::Pointer gdcmIOIn = ImageIOType::New();

  //Reference input
  if(argc == 4)
    {
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
      readerRef->SetImageIO( gdcmIOIn );
      readerRef->SetFileNames( filenames );
      readerRef->Update();
    }

  //Write output
  typedef itk::NumericSeriesFileNames             OutNamesGeneratorType;
  ImageIOType::Pointer gdcmIOOut = ImageIOType::New();
  const char * outputDirectory = argv[2];
  itksys::SystemTools::MakeDirectory( outputDirectory );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;
  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;
  typedef itk::ImageSeriesWriter< ImageType, Image2DType >  SeriesWriterType;
  OutNamesGeneratorType::Pointer namesGeneratorOut = OutNamesGeneratorType::New();

  typedef itk::MetaDataDictionary   DictionaryType;
  DictionaryType & dict = gdcmIOOut->GetMetaDataDictionary();

  //Copy image dictionary from reference if given
  if(argc == 4)
    {
      const  DictionaryType & dictRef = gdcmIOIn->GetMetaDataDictionary();
      DictionaryType::ConstIterator itr = dictRef.Begin();
      DictionaryType::ConstIterator end = dictRef.End();
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
	      
	      if(tagkey == "0008|0008")
		tagvalue = "DERIVED\\SECONDARY";
	      else if(tagkey == "0008|103e")
		tagvalue = "SEGMENTATION";
	      else if(tagkey == "0020|000e")
		tagvalue = tagvalue + ".222";
	      else if(tagkey == "0020|0011")
		tagvalue = tagvalue + ".111";
	      else if(tagkey == "0020|0012")
		tagvalue = tagvalue + ".111";
	      else if(tagkey == "0028|1052")
		tagvalue = "0";
	      else if(tagkey == "0028|1053")
	      	tagvalue = "1";
	      else if((tagkey == "0020|0013")||(tagkey == "0020|0032")||(tagkey == "0020|1041"))
	        break;		
	      itk::EncapsulateMetaData<std::string>(dict, tagkey, tagvalue);
	      //std::cout << tagkey <<  " = " << tagvalue << std::endl;
	      /*
	      tagkey = "0008|0008"; // Image Type
	      tagvalue = "DERIVED\\SECONDARY";
	      itk::EncapsulateMetaData<std::string>(dict, tagkey, tagvalue);
	      tagkey = "0008|0064"; // Conversion Type
	      tagvalue = "DV";
	      itk::EncapsulateMetaData<std::string>(dict, tagkey, tagvalue);
	      tagkey = "0008|0060"; // Modality
	      tagvalue = "MR";
	      itk::EncapsulateMetaData<std::string>(dict, tagkey, tagvalue );
	      */
	    }
	  ++itr;
	}
    }

  //Write output
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput( readerIn->GetOutput() );
  seriesWriter->SetImageIO( gdcmIOOut );
  std::string format = outputDirectory;
  format += "/image%03d.dcm";
  namesGeneratorOut->SetSeriesFormat( format.c_str() );
  namesGeneratorOut->SetStartIndex( start[2] );
  namesGeneratorOut->SetEndIndex( start[2] + size[2] - 1 );
  namesGeneratorOut->SetIncrementIndex( 1 );
  seriesWriter->SetFileNames( namesGeneratorOut->GetFileNames() );
  seriesWriter->Update();
  
  return EXIT_SUCCESS;
}

