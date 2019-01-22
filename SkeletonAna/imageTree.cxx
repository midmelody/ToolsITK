#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//ITK settings
const unsigned int Dimension = 3;
typedef unsigned char PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputImageFile outputLabelImageFile outputParentLabelImageFile" << std::endl;
      return EXIT_FAILURE;
    }
 
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writerLab = WriterType::New();
  WriterType::Pointer writerPar = WriterType::New();

  //Parameters
  reader->SetFileName( argv[2] );
  writerLab->SetFileName( argv[3] );
  writerPar->SetFileName( argv[4] );
 
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

  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate image for output
  ImageType::Pointer labelImage = ImageType::New();
  labelImage->SetRegions( region );
  labelImage->SetSpacing( spacing );
  labelImage->SetOrigin( origin );
  labelImage->SetDirection( direction );
  labelImage->Allocate();  
  ImageType::Pointer parentLabelImage = ImageType::New();
  parentLabelImage->SetRegions( region );
  parentLabelImage->SetSpacing( spacing );
  parentLabelImage->SetOrigin( origin );
  parentLabelImage->SetDirection( direction );
  parentLabelImage->Allocate(); 

  ImageType::IndexType pixelIndex;
  int i, j, k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  labelImage->SetPixel(pixelIndex, 0);
	  parentLabelImage->SetPixel(pixelIndex, 0);
	}

  FILE *treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  int label, parLabel;
  while(fscanf(treeFile, "%d %d %d %d %d", &i, &j, &k, &label, &parLabel) == 5)
    {
      pixelIndex[0]=i;
      pixelIndex[1]=j;
      pixelIndex[2]=k;
      labelImage->SetPixel(pixelIndex, label);
      parentLabelImage->SetPixel(pixelIndex, parLabel);
    }
  fclose(treeFile);



 
  //Write output  
  writerLab->SetInput( labelImage );
  try
    {
      writerLab->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }

  writerPar->SetInput( parentLabelImage );
  try
    {
      writerPar->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
