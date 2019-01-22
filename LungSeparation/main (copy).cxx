#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkNotImageFilter.h" 
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkOrImageFilter.h"

#include "itkImageFileWriter.h"

#include <iostream>

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile lungMaskFile airwayMaskFile outputImageFile [GapEnhanceImage]" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float FloatPixelType;
  typedef int CharPixelType; 
  typedef itk::Image< FloatPixelType, Dimension > FloatImageType;
  typedef itk::Image< CharPixelType, Dimension > CharImageType;
  typedef itk::ImageFileReader< FloatImageType > FloatReaderType;
  typedef itk::ImageFileReader< CharImageType > CharReaderType;
  typedef itk::BinaryBallStructuringElement< CharPixelType, Dimension > StructuringElementType;
  typedef itk::BinaryDilateImageFilter <CharImageType, CharImageType, StructuringElementType> BinaryDilateImageFilterType;
  typedef itk::BinaryErodeImageFilter <CharImageType, CharImageType, StructuringElementType> BinaryErodeImageFilterType;
  typedef itk::NotImageFilter< CharImageType, CharImageType > InvFilterType;
  typedef itk::NumericTraits< FloatPixelType >::RealType RealPixelType;
  typedef itk::SymmetricSecondRankTensor< RealPixelType, Dimension > HessianPixelType;
  typedef itk::Image< HessianPixelType, Dimension > HessianImageType;
  typedef itk::HessianToObjectnessMeasureImageFilter< HessianImageType,FloatImageType > ObjectnessFilterType;
  typedef itk::MultiScaleHessianBasedMeasureImageFilter< FloatImageType, HessianImageType, FloatImageType > MultiScaleEnhancementFilterType;
  typedef itk::RescaleIntensityImageFilter< FloatImageType, CharImageType > RescalerType;
  typedef itk::BinaryThresholdImageFilter< CharImageType, CharImageType > BinaryThresholdFilterType; 
  typedef itk::AndImageFilter< CharImageType, CharImageType, CharImageType > AndFilterType;
  typedef itk::ConnectedComponentImageFilter< CharImageType, CharImageType > ConnectedComponentFilterType;
  typedef itk::RelabelComponentImageFilter< CharImageType, CharImageType >  RelabelFilterType;
  typedef itk::SignedMaurerDistanceMapImageFilter < CharImageType, FloatImageType > DistanceTransformFilterType;
  typedef itk::OrImageFilter< CharImageType, CharImageType, CharImageType > OrFilterType;
  typedef itk::ImageFileWriter< CharImageType > WriterType;

  //Filters
  FloatReaderType::Pointer readerI = FloatReaderType::New();
  CharReaderType::Pointer readerL = CharReaderType::New();
  CharReaderType::Pointer readerA = CharReaderType::New();
  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
  InvFilterType::Pointer invFilter = InvFilterType::New();
  ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
  MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
  RescalerType::Pointer rescaleFilter = RescalerType::New();
  BinaryThresholdFilterType::Pointer thresholdFilter = BinaryThresholdFilterType::New();
  AndFilterType::Pointer andFilter = AndFilterType::New();
  ConnectedComponentFilterType::Pointer connectComponentFilter = ConnectedComponentFilterType::New();
  RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  DistanceTransformFilterType::Pointer DTFilter = DistanceTransformFilterType::New();
  OrFilterType::Pointer orFilter = OrFilterType::New();
  WriterType::Pointer writer = WriterType::New();
  WriterType::Pointer writerTemp = WriterType::New();

  //Parameters
  readerI->SetFileName( argv[1] );
  readerL->SetFileName( argv[2] );
  readerA->SetFileName( argv[3] );
  writer->SetFileName( argv[4] );
  if(argc == 6)
    writerTemp->SetFileName( argv[5] );
  //Read images - original CT image, lung mask, airway mask
  readerI->Update();
  readerL->Update();  
  readerA->Update(); 

  //Subtract airway - dilate for ensurance
  int radius = 2;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();
  dilateFilter->SetInput(readerA->GetOutput());
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->SetBackgroundValue(0);
  dilateFilter->SetForegroundValue(1);
  invFilter->SetInput(dilateFilter->GetOutput());
  andFilter->SetInput1(invFilter->GetOutput());
  andFilter->SetInput2(readerL->GetOutput());
  andFilter->Update();
  CharImageType::Pointer lungWOairway = andFilter->GetOutput();
  lungWOairway->DisconnectPipeline();
  
  //Calculate gap enhancement
  objectnessFilter->SetBrightObject(true);
  objectnessFilter->SetAlpha(0.5);
  objectnessFilter->SetBeta(0.5);
  objectnessFilter->SetGamma(5.0);
  objectnessFilter->SetObjectDimension(2);
  multiScaleEnhancementFilter->SetSigmaStepMethodToLogarithmic();
  multiScaleEnhancementFilter->SetSigmaMinimum( 0.1 );
  multiScaleEnhancementFilter->SetSigmaMaximum( 1 );
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps( 10 );
  multiScaleEnhancementFilter->SetInput( readerI->GetOutput() );
  multiScaleEnhancementFilter->SetHessianToMeasureFilter( objectnessFilter );
  multiScaleEnhancementFilter->GenerateScalesOutputOn();
  multiScaleEnhancementFilter->SetGenerateScalesOutput( true );
  multiScaleEnhancementFilter->SetGenerateScalesOutput( false );
  multiScaleEnhancementFilter->NonNegativeHessianBasedMeasureOff(); 
  rescaleFilter->SetInput(multiScaleEnhancementFilter->GetOutput());
  rescaleFilter->SetOutputMinimum(   0 );
  rescaleFilter->SetOutputMaximum( 255 );
  rescaleFilter->Update();
  CharImageType::Pointer gapEnh = rescaleFilter->GetOutput();
  if(argc == 6)
    {
      writerTemp->SetInput( gapEnh );
      writerTemp->Update();
    }
  gapEnh->DisconnectPipeline();

  //Subtract gap until 2 groups
  int groupCt = 1;
  int thre = 20;
  int dilateSize = 1;
  CharImageType::Pointer lungGroup = CharImageType::New();
  StructuringElementType gapSE;
  gapSE.SetRadius(dilateSize);
  gapSE.CreateStructuringElement();
  std::cout<<"First"<<std::endl;
  while((groupCt==1)&&(thre>2))
    {
      std::cout<<thre<<std::endl;
      thresholdFilter->SetInput(gapEnh);
      thresholdFilter->SetOutsideValue( 0 );
      thresholdFilter->SetInsideValue( 1 );
      thresholdFilter->SetLowerThreshold( thre );
      thresholdFilter->SetUpperThreshold( 255 );
      dilateFilter->SetKernel(gapSE);
      dilateFilter->SetBackgroundValue(0);
      dilateFilter->SetForegroundValue(1);
      dilateFilter->SetInput(thresholdFilter->GetOutput());   
      invFilter->SetInput(dilateFilter->GetOutput());
      andFilter->SetInput1(invFilter->GetOutput());
      andFilter->SetInput2(lungWOairway);
      connectComponentFilter->SetInput( andFilter->GetOutput() );
      connectComponentFilter->SetFullyConnected(0);
      relabelFilter->SetInput( connectComponentFilter->GetOutput() );
      relabelFilter->SetMinimumObjectSize(1000);
      relabelFilter->Update();
      lungGroup = relabelFilter->GetOutput();
      lungGroup->DisconnectPipeline();
      typedef std::vector< itk::SizeValueType > SizesInPixelsType;
      const SizesInPixelsType &  sizesInPixels = relabelFilter->GetSizeOfObjectsInPixels();
      SizesInPixelsType::const_iterator sizeItr = sizesInPixels.begin();
      SizesInPixelsType::const_iterator sizeBegin = sizesInPixels.begin();
      SizesInPixelsType::const_iterator sizeEnd = sizesInPixels.end();

      std::cout << "Number of pixels per class " << std::endl;
      unsigned int kclass = 1;
      while( sizeItr != sizeEnd )
	{
	  std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
	  ++kclass;
	  ++sizeItr;
	}
      groupCt = kclass-1;
      thre = thre-2;
      /*
      if(groupCt > 1)
	{
	  float ratio;
	  float i1, i2;
	  i1 = *sizeBegin;
	  i2 = *(++sizeBegin);
	  ratio = i1/i2;
	  if(ratio>5)
	    groupCt = 1;
	}
      */
    }

  std::cout<<"Second"<<std::endl;
  thre = 10;
  dilateSize = 2;
  gapSE.SetRadius(dilateSize);
  gapSE.CreateStructuringElement();
  while((groupCt==1)&&(thre>2))
    {
      std::cout<<thre<<std::endl;
      thresholdFilter->SetInput(gapEnh);
      thresholdFilter->SetOutsideValue( 0 );
      thresholdFilter->SetInsideValue( 1 );
      thresholdFilter->SetLowerThreshold( thre );
      thresholdFilter->SetUpperThreshold( 255 );
      dilateFilter->SetKernel(gapSE);
      dilateFilter->SetBackgroundValue(0);
      dilateFilter->SetForegroundValue(1);
      dilateFilter->SetInput(thresholdFilter->GetOutput());   
      invFilter->SetInput(dilateFilter->GetOutput());
      andFilter->SetInput1(invFilter->GetOutput());
      andFilter->SetInput2(lungWOairway);
      connectComponentFilter->SetInput( andFilter->GetOutput() );
      connectComponentFilter->SetFullyConnected(0);
      relabelFilter->SetInput( connectComponentFilter->GetOutput() );
      relabelFilter->SetMinimumObjectSize(1000);
      relabelFilter->Update();
      lungGroup = relabelFilter->GetOutput();
      lungGroup->DisconnectPipeline();
      typedef std::vector< itk::SizeValueType > SizesInPixelsType;
      const SizesInPixelsType &  sizesInPixels = relabelFilter->GetSizeOfObjectsInPixels();
      SizesInPixelsType::const_iterator sizeItr = sizesInPixels.begin();
      SizesInPixelsType::const_iterator sizeBegin = sizesInPixels.begin();
      SizesInPixelsType::const_iterator sizeEnd = sizesInPixels.end();
      std::cout << "Number of pixels per class " << std::endl;
      unsigned int kclass = 1;
      while( sizeItr != sizeEnd )
	{
	  std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
	  ++kclass;
	  ++sizeItr;
	}
      groupCt = kclass-1;
      thre--;
      if(groupCt >1 )
	{
	  float ratio;
	  float i1, i2;
	  i1 = *sizeBegin;
	  i2 = *(++sizeBegin);
	  ratio = i1/i2;
	  if(ratio>5)
	    groupCt = 1;
	}
    } 

  std::cout<<"Third"<<std::endl;
  thre = 5;
  dilateSize = 3;
  gapSE.SetRadius(dilateSize);
  gapSE.CreateStructuringElement();
  while((groupCt==1)&&(thre>2))
    {
      std::cout<<thre<<std::endl;
      thresholdFilter->SetInput(gapEnh);
      thresholdFilter->SetOutsideValue( 0 );
      thresholdFilter->SetInsideValue( 1 );
      thresholdFilter->SetLowerThreshold( thre );
      thresholdFilter->SetUpperThreshold( 255 );
      dilateFilter->SetKernel(gapSE);
      dilateFilter->SetBackgroundValue(0);
      dilateFilter->SetForegroundValue(1);
      dilateFilter->SetInput(thresholdFilter->GetOutput());   
      invFilter->SetInput(dilateFilter->GetOutput());
      andFilter->SetInput1(invFilter->GetOutput());
      andFilter->SetInput2(lungWOairway);
      connectComponentFilter->SetInput( andFilter->GetOutput() );
      connectComponentFilter->SetFullyConnected(0);
      relabelFilter->SetInput( connectComponentFilter->GetOutput() );
      relabelFilter->SetMinimumObjectSize(1000);
      relabelFilter->Update();
      lungGroup = relabelFilter->GetOutput();
      lungGroup->DisconnectPipeline();
      typedef std::vector< itk::SizeValueType > SizesInPixelsType;
      const SizesInPixelsType &  sizesInPixels = relabelFilter->GetSizeOfObjectsInPixels();
      SizesInPixelsType::const_iterator sizeItr = sizesInPixels.begin();
      SizesInPixelsType::const_iterator sizeBegin = sizesInPixels.begin();
      SizesInPixelsType::const_iterator sizeEnd = sizesInPixels.end();
      std::cout << "Number of pixels per class " << std::endl;
      unsigned int kclass = 1;
      while( sizeItr != sizeEnd )
	{
	  std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
	  ++kclass;
	  ++sizeItr;
	}
      groupCt = kclass-1;
      thre--;
      if(groupCt >1)
	{
	  float ratio;
	  float i1, i2;
	  i1 = *sizeBegin;
	  i2 = *(++sizeBegin);
	  ratio = i1/i2;
	  if(ratio>5)
	    groupCt = 1;
	}
    } 

  CharImageType::Pointer lungL = CharImageType::New();
  CharImageType::Pointer lungR = CharImageType::New();
  FloatImageType::Pointer lungLDT = FloatImageType::New();
  FloatImageType::Pointer lungRDT = FloatImageType::New();
  
  if(groupCt!=1)
    {
      //Write output  
      writer->SetInput( lungGroup );
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
    }
  else
    {
      std::cerr << "Failed to separate " << std::endl;
      return EXIT_FAILURE;
    }     
 

  return EXIT_SUCCESS;
}
