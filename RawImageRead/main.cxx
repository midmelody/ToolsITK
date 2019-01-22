#include "itkImage.h"
#include "itkRawImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

float SearchWordInFile(char *fname, char *str) 
{
  FILE *fp;
  int lineNum = 1;
  bool findResult = false;
  char temp[512];
  char tempT[20];
  float value;
  if((fp = fopen(fname, "r")) == NULL) 
    {
      std::cerr<<"No File Found"<<std::endl;
      return(-1);
    }

  while((fgets(temp, 512, fp) != NULL) && !findResult)
    {
      if((strstr(temp, str)) != NULL) 
	{
	  sscanf(temp, "%s %f", tempT, &value);
	  findResult = true;
	}
      lineNum++;
    }
  
  if(!findResult) 
    {
      std::cout<<"No match found"<<std::endl;
    }	
  fclose(fp);
  return(value);
}

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << "HeaderFile RawImageFile OutputImageFile";
      std::cerr << " " << std::endl;
      return EXIT_FAILURE;
    }
  char keyword[20];
  
  int dimX, dimY, dimZ;
  strcpy(keyword, "x_dimension");
  dimX = SearchWordInFile(argv[1], keyword);
  strcpy(keyword, "y_dimension");
  dimY = SearchWordInFile(argv[1], keyword);
  strcpy(keyword, "z_dimension");
  dimZ = SearchWordInFile(argv[1], keyword);
  std::cout<<"Size "<<dimX<<" "<<dimY<<" "<<dimZ<<std::endl;

  float oriX, oriY, oriZ;
  strcpy(keyword, "volume_origin_x");
  oriX = SearchWordInFile(argv[1], keyword);
  strcpy(keyword, "volume_origin_y");
  oriY = SearchWordInFile(argv[1], keyword);
  strcpy(keyword, "volume_origin_z");
  oriZ = SearchWordInFile(argv[1], keyword);
  std::cout<<"Origin "<<oriX<<" "<<oriY<<" "<<oriZ<<std::endl;

  float spacX, spacY, spacZ;
  strcpy(keyword, "pixel_size_x");
  spacX = SearchWordInFile(argv[1], keyword);
  strcpy(keyword, "pixel_size_y");
  spacY = SearchWordInFile(argv[1], keyword);
  strcpy(keyword, "pixel_size_z");
  spacZ = SearchWordInFile(argv[1], keyword);
  std::cout<<"Spacing "<<spacX<<" "<<spacY<<" "<<spacZ<<std::endl;

  typedef float PixelType;
  typedef itk::RawImageIO< PixelType, 3 > RawIOType;
  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ImageType::SpacingType spacing; 
  ImageType::PointType origin; 
  ImageType::SizeType size;
  size[0] = dimX;
  size[1] = dimY;
  size[2] = dimZ;
  spacing[0] = spacX;
  spacing[1] = spacY;
  spacing[2] = spacZ;
  origin[0] = oriX;
  origin[1] = oriY;
  origin[2] = oriZ;

  RawIOType::Pointer rawIO = RawIOType::New();
  rawIO->SetFileName( argv[2] );
  for(int i=0;i<3; i++)
    {
      rawIO->SetDimensions(i, size[i]);
      rawIO->SetSpacing(i, spacing[i]);
      rawIO->SetOrigin(i, origin[i]);
    }
  rawIO->SetHeaderSize( 0 );
  rawIO->SetByteOrderToLittleEndian();
  rawIO->SetPixelType(itk::ImageIOBase::SCALAR);
  rawIO->SetNumberOfComponents( 1 );

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->SetImageIO(rawIO);

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

