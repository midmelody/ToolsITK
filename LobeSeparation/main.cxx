#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkNotImageFilter.h" 
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
      std::cerr << argv[0] << " probImageFile lungMaskFile outputImageFile expectedCT" << std::endl;
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
  typedef itk::RescaleIntensityImageFilter< FloatImageType, CharImageType > RescalerType;
  typedef itk::BinaryThresholdImageFilter< CharImageType, CharImageType > BinaryThresholdFilterType;
  typedef itk::BinaryThresholdImageFilter< FloatImageType, CharImageType > FloatThresholdFilterType; 

  typedef itk::AndImageFilter< CharImageType, CharImageType, CharImageType > AndFilterType;
  typedef itk::ConnectedComponentImageFilter< CharImageType, CharImageType > ConnectedComponentFilterType;
  typedef itk::RelabelComponentImageFilter< CharImageType, CharImageType >  RelabelFilterType;
  typedef itk::SignedMaurerDistanceMapImageFilter < CharImageType, FloatImageType > DistanceTransformFilterType;
  typedef itk::OrImageFilter< CharImageType, CharImageType, CharImageType > OrFilterType;
  typedef itk::ImageFileWriter< CharImageType > WriterType;

  //Filters
  FloatReaderType::Pointer readerPro = FloatReaderType::New();
  CharReaderType::Pointer readerLng = CharReaderType::New();
  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
  InvFilterType::Pointer invFilter = InvFilterType::New();
  RescalerType::Pointer rescaleFilter = RescalerType::New();
  BinaryThresholdFilterType::Pointer thresholdFilter = BinaryThresholdFilterType::New();
  FloatThresholdFilterType::Pointer probThresholdFilter = FloatThresholdFilterType::New();
  AndFilterType::Pointer andFilter = AndFilterType::New();
  ConnectedComponentFilterType::Pointer connectComponentFilter = ConnectedComponentFilterType::New();
  RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  DistanceTransformFilterType::Pointer DTFilter = DistanceTransformFilterType::New();
  OrFilterType::Pointer orFilter = OrFilterType::New();
  WriterType::Pointer writer = WriterType::New();
  WriterType::Pointer writerTemp = WriterType::New();

  //Parameters
  readerPro->SetFileName( argv[1] );
  readerLng->SetFileName( argv[2] );
  writer->SetFileName( argv[3] );
  int expCT;
  expCT = atoi( argv[4] );
  //Read images - fissure probability image, lung mask
  readerPro->Update();
  readerLng->Update();  
 
  //Subtract gap until 2 groups
  int groupCt = 1;
  float thre = 30;
  int dilateSize = 1;
  CharImageType::Pointer lungGroup = CharImageType::New();
  StructuringElementType gapSE;
  gapSE.SetRadius(dilateSize);
  gapSE.CreateStructuringElement();
  std::cout<<"First"<<std::endl;
  while((groupCt<expCT)&&(thre>=10))
    {
      std::cout<<thre<<std::endl;
      probThresholdFilter->SetInput(readerPro->GetOutput());
      probThresholdFilter->SetOutsideValue( 0 );
      probThresholdFilter->SetInsideValue( 1 );
      probThresholdFilter->SetLowerThreshold( thre );
      probThresholdFilter->SetUpperThreshold( 255 );
      
      dilateFilter->SetKernel(gapSE);
      dilateFilter->SetBackgroundValue(0);
      dilateFilter->SetForegroundValue(1);
      dilateFilter->SetInput(probThresholdFilter->GetOutput());   
      invFilter->SetInput(dilateFilter->GetOutput());
      andFilter->SetInput1(invFilter->GetOutput());
      andFilter->SetInput2(readerLng->GetOutput());
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
      thre = thre-10;
      if(groupCt >= expCT)
	{
	  float ratio;
	  float i1, i2;
	  i1 = *sizeBegin;
	  i2 = *(sizeBegin+expCT-1);
	  ratio = i1/i2;
	  if(ratio>20)
	    groupCt = 1;
	}
      
    }

  std::cout<<"Second"<<std::endl;
  thre = 30;
  dilateSize = 2;
  gapSE.SetRadius(dilateSize);
  gapSE.CreateStructuringElement();
  while((groupCt<expCT)&&(thre>=5))
    {
      std::cout<<thre<<std::endl;
      probThresholdFilter->SetInput(readerPro->GetOutput());
      probThresholdFilter->SetOutsideValue( 0 );
      probThresholdFilter->SetInsideValue( 1 );
      probThresholdFilter->SetLowerThreshold( thre );
      probThresholdFilter->SetUpperThreshold( 255 );
      dilateFilter->SetKernel(gapSE);
      dilateFilter->SetBackgroundValue(0);
      dilateFilter->SetForegroundValue(1);
      dilateFilter->SetInput(probThresholdFilter->GetOutput());   
      invFilter->SetInput(dilateFilter->GetOutput());
      andFilter->SetInput1(invFilter->GetOutput());
      andFilter->SetInput2(readerLng->GetOutput());
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
      thre = thre-5;
      if(groupCt >= expCT)
	{
	  float ratio;
	  float i1, i2;
	  i1 = *sizeBegin;
	  i2 = *(sizeBegin+expCT-1);
	  ratio = i1/i2;
	  if(ratio>20)
	    groupCt = 1;
	}
    } 

  std::cout<<"Third"<<std::endl;
  thre = 30;
  dilateSize = 3;
  gapSE.SetRadius(dilateSize);
  gapSE.CreateStructuringElement();
  while((groupCt<expCT)&&(thre>=5))
    {
      std::cout<<thre<<std::endl;
      probThresholdFilter->SetInput(readerPro->GetOutput());
      probThresholdFilter->SetOutsideValue( 0 );
      probThresholdFilter->SetInsideValue( 1 );
      probThresholdFilter->SetLowerThreshold( thre );
      probThresholdFilter->SetUpperThreshold( 255 );
      dilateFilter->SetKernel(gapSE);
      dilateFilter->SetBackgroundValue(0);
      dilateFilter->SetForegroundValue(1);
      dilateFilter->SetInput(probThresholdFilter->GetOutput());   
      invFilter->SetInput(dilateFilter->GetOutput());
      andFilter->SetInput1(invFilter->GetOutput());
      andFilter->SetInput2(readerLng->GetOutput());
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
      thre = thre-5;
      
      if(groupCt >= expCT)
	{
	  float ratio;
	  float i1, i2;
	  i1 = *sizeBegin;
	  i2 = *(sizeBegin+expCT-1);
	  ratio = i1/i2;
	  if(ratio>20)
	    groupCt = 1;
	}
    } 
  
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
