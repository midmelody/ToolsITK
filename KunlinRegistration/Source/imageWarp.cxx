#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <time.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkBinaryThresholdImageFilter.h"

#include "ImageType.h"
#include "IMGIO.h"
#include "BasicProcess.h"
#include "ArrayProcess.h"

using namespace std;

#define bgValue -1000
#define dbgValue 0
#define debugTag 0

int main( int argc, char* argv[] )
{
  if( argc < 6 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile dispalcementX dispalcementY dispalcementZ outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  string imgName = argv[1];
  string dpXName = argv[2];
  string dpYName = argv[3];
  string dpZName = argv[4];
  string outName = argv[5];
  
  ImageType::Pointer iimageIn = ImageIO.ReadImage< ImageType >( imgName );
  DispImageType::Pointer iimageDX = ImageIO.ReadImage< DispImageType >( dpXName );
  DispImageType::Pointer iimageDY = ImageIO.ReadImage< DispImageType >( dpYName );
  DispImageType::Pointer iimageDZ = ImageIO.ReadImage< DispImageType >( dpZName );

  ImageType::DirectionType direction = iimageIn->GetDirection();
  ImageType::SizeType size = iimageIn->GetRequestedRegion().GetSize();
  ImageType::SpacingType spacing = iimageIn->GetSpacing();
  ImageType::IndexType iStart;
  iStart.Fill(0);
  ImageType::PointType origin = iimageIn->GetOrigin();
  ImageType::RegionType region;
  region.SetIndex(iStart);
  region.SetSize(size);

  int xdimo=size[0];
  int ydimo=size[1];
  int zdimo=size[2];
  int imgSize=xdimo*ydimo*zdimo;

  float spacing_x = static_cast<float>(spacing[0]);
  float spacing_y = static_cast<float>(spacing[1]);
  float spacing_z = static_cast<float>(spacing[2]);
    
  //initialize arrays
  std::cout<<"Initializing arrays ..."<<std::endl; 
  float* imageIn = new float[imgSize];
  float* imageDX = new float[imgSize];
  float* imageDY = new float[imgSize];
  float* imageDZ = new float[imgSize];
  float* imageOut= new float[imgSize];
  //convert input images to arrays
  ITKImage2Array< ImageType >( iimageIn, imageIn );
  DispITKImage2Array< DispImageType >( iimageDX, iimageDY, iimageDZ, imageDX, imageDY, imageDZ);
  
  //generate deformed image
  std::cout<<"Calculating and saving deformed image"<<std::endl;
  deformed2_inner(imageOut, imageIn, xdimo, ydimo, zdimo, imageDX, imageDY, imageDZ, bgValue);
  ImageType::Pointer iimageOut = ImageIO.initImage< ImageType >( iStart, size, direction, origin, spacing );
  Array2ITKImage< ImageType >( imageOut, iimageOut );
  ImageIO.WriteImage< ImageType >( iimageOut, outName );

  return 0;
}


