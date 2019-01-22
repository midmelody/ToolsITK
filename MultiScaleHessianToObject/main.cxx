#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

//The program computes multiscale objectiveness for 3D image

int main( int argc, char * argv[] )
{
  if( argc < 9 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile enhancedImageFile scalesImageFile brightObjFlag minSigma maxSigma numOfSteps objectDimension (0-blob 1-vessel 2-plate)" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef unsigned char OutputPixelType;
  typedef float FloatPixelType;
  typedef itk::Image< FloatPixelType, Dimension > FloatImageType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::NumericTraits< PixelType >::RealType RealPixelType;
  typedef itk::SymmetricSecondRankTensor< RealPixelType, Dimension > HessianPixelType;
  typedef itk::Image< HessianPixelType, Dimension > HessianImageType;
  typedef itk::HessianToObjectnessMeasureImageFilter< HessianImageType,ImageType > ObjectnessFilterType;
  typedef itk::MultiScaleHessianBasedMeasureImageFilter< ImageType, HessianImageType, ImageType > MultiScaleEnhancementFilterType;
  typedef itk::RescaleIntensityImageFilter< ImageType, OutputImageType > RescaleFilterType;

  typedef itk::RescaleIntensityImageFilter< FloatImageType, OutputImageType > RescaleScaleFilterType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
  MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  RescaleScaleFilterType::Pointer rescaleScaleFilter = RescaleScaleFilterType::New();
  WriterType::Pointer writer1 = WriterType::New();
  WriterType::Pointer writer2 = WriterType::New();
  //Parameters
  bool brightObjFlag = atoi(argv[4]);
  reader->SetFileName( argv[1] );
  if(brightObjFlag)
    objectnessFilter->SetBrightObject(true);
  else
    objectnessFilter->SetBrightObject(false);  
  objectnessFilter->SetAlpha(0.5);
  objectnessFilter->SetBeta(0.5);
  objectnessFilter->SetGamma(5.0);
  objectnessFilter->SetObjectDimension(atoi(argv[8]));
  multiScaleEnhancementFilter->SetSigmaStepMethodToLogarithmic();
  multiScaleEnhancementFilter->SetSigmaMinimum( atof(argv[5]) );
  multiScaleEnhancementFilter->SetSigmaMaximum( atof(argv[6]) );
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps( atoi(argv[7]) );
  writer1->SetFileName( argv[2] );
 
  //Pipeline
  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout << "Problems reading input image" << std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }

  multiScaleEnhancementFilter->SetInput( reader->GetOutput() );
  multiScaleEnhancementFilter->SetHessianToMeasureFilter( objectnessFilter );


  multiScaleEnhancementFilter->GenerateScalesOutputOn();
  if ( !multiScaleEnhancementFilter->GetGenerateScalesOutput() ) 
    {
    std::cerr << "Error in Set/GetGenerateScalesOutput()" << std::endl;
    return EXIT_FAILURE;
    }
  multiScaleEnhancementFilter->SetGenerateScalesOutput( false );
  if ( multiScaleEnhancementFilter->GetGenerateScalesOutput() ) 
    {
    std::cerr << "Error in Set/GetGenerateScalesOutput()" << std::endl;
    return EXIT_FAILURE;
    }
  multiScaleEnhancementFilter->SetGenerateScalesOutput( true );
  
  multiScaleEnhancementFilter->GenerateHessianOutputOn();
  if ( !multiScaleEnhancementFilter->GetGenerateHessianOutput() ) 
    {
    std::cerr << "Error in Set/GetGenerateHessianOutput()" << std::endl;
    return EXIT_FAILURE;
    }
  multiScaleEnhancementFilter->SetGenerateHessianOutput( false );
  if ( multiScaleEnhancementFilter->GetGenerateHessianOutput() ) 
    {
    std::cerr << "Error in Set/GetGenerateHessianOutput()" << std::endl;
    return EXIT_FAILURE;
    }
  multiScaleEnhancementFilter->SetGenerateHessianOutput( true );


  multiScaleEnhancementFilter->NonNegativeHessianBasedMeasureOff(); 






  try
    {
      multiScaleEnhancementFilter->Update();
    }
  catch (itk::ExceptionObject &err)
    {
      std::cout << "Problems multiscale processing" << std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
    }


  rescaleFilter->SetInput(multiScaleEnhancementFilter->GetOutput());
  rescaleFilter->SetOutputMinimum(   0 );
  rescaleFilter->SetOutputMaximum( 255 );
  writer1->SetInput(rescaleFilter->GetOutput());

  try
    {
      writer1->Update();
    }  
  catch( itk::ExceptionObject & err )
    {
      std::cout << "Problems writing image" << std::endl;     
      std::cout << "ExceptionObject caught !" << std::endl;
      std::cout << err << std::endl;
      return EXIT_FAILURE;
    }
 

  writer2->SetFileName(argv[3]);
  //writer->SetInput(multiScaleEnhancementFilter->GetScalesOutput()); 
  rescaleScaleFilter->SetInput(multiScaleEnhancementFilter->GetScalesOutput()); 
  rescaleScaleFilter->SetOutputMinimum(  5 );
  rescaleScaleFilter->SetOutputMaximum( 250 );
  writer2->SetInput(rescaleScaleFilter->GetOutput());

  try
    {
      writer2->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }



  return EXIT_SUCCESS;
}
