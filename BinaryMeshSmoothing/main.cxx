#include "itkImage.h"
#include "itkMesh.h"
#include "itkImageFileReader.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkTriangleMeshToBinaryImageFilter.h"
#include "itkVTKPolyDataWriter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;

  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Mesh< double > MeshType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::BinaryMask3DMeshSource< ImageType, MeshType > MeshSourceType;
  
  typedef itk::TriangleMeshToBinaryImageFilter< MeshType, ImageType > MeshToBinaryFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::VTKPolyDataWriter<MeshType> MeshWriterType;
 
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  MeshSourceType::Pointer meshSource = MeshSourceType::New();
  MeshToBinaryFilterType::Pointer meshBinFilter = MeshToBinaryFilterType::New();
  WriterType::Pointer writer = WriterType::New();
  MeshWriterType::Pointer meshWriter = MeshWriterType::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

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

  meshSource->SetInput( reader->GetOutput() ); 
  meshSource->SetObjectValue( 1 );
  meshSource->Update();



  meshBinFilter->SetTolerance (1.0);
  meshBinFilter->SetSpacing(reader->GetOutput()->GetSpacing());
  meshBinFilter->SetOrigin(reader->GetOutput()->GetOrigin());
  meshBinFilter->SetSize(reader->GetOutput()->GetRequestedRegion().GetSize());
  meshBinFilter->SetDirection(reader->GetOutput()->GetDirection());
  meshBinFilter->SetInput(meshSource->GetOutput());
  meshBinFilter->UpdateLargestPossibleRegion();

  //Write output  
  meshWriter->SetInput(meshSource->GetOutput());
  meshWriter->SetFileName("test.vtk");
  meshWriter->Update();


  writer->SetInput( meshBinFilter->GetOutput() );
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
