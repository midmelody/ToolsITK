#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNotImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 
int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " refCenterFile seedCenterFile inputImageFile outputImageFile "<< std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > BinaryThresholdFilterType; 
  typedef itk::NotImageFilter< ImageType, ImageType > InvFilterType;
  typedef itk::AndImageFilter< ImageType, ImageType, ImageType > AndFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  BinaryThresholdFilterType::Pointer maskThresholdWhole = BinaryThresholdFilterType::New();
  BinaryThresholdFilterType::Pointer maskThresholdTarget = BinaryThresholdFilterType::New();
  InvFilterType::Pointer invFilter = InvFilterType::New();
  AndFilterType::Pointer andFilter = AndFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[3] );
  writer->SetFileName( argv[4] );
 
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
  
  //Read reference
  int refX, refY, refZ;
  FILE *refFile = fopen(argv[1], "r");
  fscanf(refFile, "%d %d %d\n", &refX, &refY, &refZ);
  fclose(refFile);

  //Read output
  int distance;
  float minDistance = 100000000000;
  int minLoc;
  FILE *candFile = fopen(argv[2], "r");
  int groupNum = 1;
  int candX, candY, candZ;
  while(fscanf(candFile, "%d %d %d\n", &candX, &candY, &candZ) == 3)  
    {  
      distance = (refX-candX)*(refX-candX);
      distance = distance + (refY-candY)*(refY-candY);
      distance = distance + (refZ-candZ)*(refZ-candZ);
      if(minDistance>distance)
	{
	  minDistance = distance;
	  minLoc = groupNum;
	}
      groupNum++;
    }
  fclose(candFile);

  if(minDistance<2500)
    {
      std::cout<<"Group "<<minLoc<<" is the rest part of lung1"<<std::endl;
      //Get rid of the false positive
      maskThresholdTarget->SetInput(reader->GetOutput());
      maskThresholdTarget->SetOutsideValue( 0 );
      maskThresholdTarget->SetInsideValue( 1 );
      maskThresholdTarget->SetLowerThreshold( minLoc );
      maskThresholdTarget->SetUpperThreshold( minLoc ); 
      invFilter->SetInput(maskThresholdTarget->GetOutput());
      maskThresholdWhole->SetInput(reader->GetOutput());
      maskThresholdWhole->SetOutsideValue( 0 );
      maskThresholdWhole->SetInsideValue( 1 );
      maskThresholdWhole->SetLowerThreshold( 1 );
      maskThresholdWhole->SetUpperThreshold( 100 ); 
      andFilter->SetInput1(maskThresholdWhole->GetOutput());
      andFilter->SetInput2(invFilter->GetOutput());
      writer->SetInput( andFilter->GetOutput() );
    }
  else
    {
      minLoc = 0;
      std::cout<<"No rest part of lung1"<<std::endl;
      writer->SetInput( reader->GetOutput() );
    }


 
  //Write output  
 
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
