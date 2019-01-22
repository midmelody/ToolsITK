// The program reads a 3D image from a non DICOM file and writes it as a series of DICOM slices
// Output header comes from a Dicom reference image

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"
#include <vector>
#include "itksys/SystemTools.hxx"

int main( int argc, char* argv[] )
{   
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " InputImage ReferenceDicomDirectory OutputDicomDirectory" << std::endl;
      return EXIT_FAILURE;
    }

  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;
  typedef itk::Image< PixelType, Dimension >      ImageType;
  typedef itk::ImageFileReader< ImageType >       InReaderType;
  typedef itk::ImageSeriesReader< ImageType >     RefReaderType;
  typedef itk::GDCMImageIO                        ImageIOType;
  typedef itk::GDCMSeriesFileNames                NamesGeneratorType;

  //Image input
  InReaderType::Pointer readerIn = InReaderType::New();
  readerIn->SetFileName( argv[1] );
  readerIn->Update();
  ImageType::RegionType region = readerIn->GetOutput()->GetLargestPossibleRegion();
  ImageType::SizeType  size  = region.GetSize();

  //Reference input
  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
  namesGenerator->SetInputDirectory( argv[2] );
  const RefReaderType::FileNamesContainer & filenames = namesGenerator->GetInputFileNames();
  const unsigned int numberOfFileNames =  filenames.size();
  if(size[2] != numberOfFileNames)
    {
      std::cerr << "Input and reference images do not match!" << std::endl;
      return EXIT_FAILURE;
    }
  else
    std::cout << "Total slice number: "<< numberOfFileNames << std::endl;
  RefReaderType::Pointer readerRef = RefReaderType::New();
  readerRef->SetImageIO( gdcmIO );
  readerRef->SetFileNames( filenames );
  readerRef->Update();

  //Write output
  const char * outputDirectory = argv[3];
  itksys::SystemTools::MakeDirectory( outputDirectory );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;
  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;
  typedef itk::ImageSeriesWriter< ImageType, Image2DType >  SeriesWriterType;
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput( readerIn->GetOutput() );
  seriesWriter->SetImageIO( gdcmIO );
  namesGenerator->SetOutputDirectory( outputDirectory );
  seriesWriter->SetFileNames( namesGenerator->GetOutputFileNames() );
  seriesWriter->SetMetaDataDictionaryArray( readerRef->GetMetaDataDictionaryArray() );
  seriesWriter->Update();
  
  return EXIT_SUCCESS;
}
