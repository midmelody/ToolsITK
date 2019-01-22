#include "itkImage.h"
#include "itkImageFileReader.h"

#include <iostream>

//The program counts the voxel in a label image for volume of all label groups

int main( int argc, char * argv[] )
{
  if( argc < 2 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();

  //Parameters
  reader->SetFileName( argv[1] );
 
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
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  size[0] = size[0];
  size[1] = size[1];
  size[2] = size[2];
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
  float voxelSize = spacing[0]*spacing[1]*spacing[2];

  std::map<int, int> labelCt;
  ImageType::IndexType index;
  ImageType::PixelType pixel;
  int i,j,k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  index[0] = i;
	  index[1] = j;
	  index[2] = k;
	  pixel = reader->GetOutput()->GetPixel( index );
	  if(pixel!=0)
	    labelCt[pixel]++;
	}

  std::cout << "Label Group Count: " << labelCt.size() << std::endl;
  std::cout << "Image Spacing: "<<spacing[0]<<" "<<spacing[1]<<" "<<spacing[2]<<std::endl;
  std::map<int,int>::iterator ii;
  for( ii=labelCt.begin(); ii!=labelCt.end(); ++ii)
    {
      std::cout<<"Group "<<(*ii).first<<": "<<(*ii).second<<" voxels; "<<float((*ii).second)*voxelSize/1000<<" mL"<<std::endl;
    }
  return EXIT_SUCCESS;
}
