#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkCoocurrenceTextureFeaturesImageFilter.h"

#include <iostream>

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile maskImageFile outputTxtFile" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 2;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::Statistics::CoocurrenceTextureFeaturesImageFilter< ImageType, ImageType > TextureFilterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerM = ReaderType::New();
  TextureFilterType::Pointer textureFilter = TextureFilterType::New();
  
  //Parameters
  reader->SetFileName( argv[1] );
  readerM->SetFileName( argv[2] );
  FILE *featureFile = fopen(argv[3], "w");
  
  //Pipeline
  reader->Update();
  readerM->Update();
  textureFilter->SetInput(reader->GetOutput());
  textureFilter->SetMaskImage(readerM->GetOutput());
  textureFilter->SetInsidePixelValue(255);
  textureFilter->SetNumberOfBinsPerAxis(256);
  textureFilter->Update();

  //Write output
  const TextureFilterType::FeatureValueVector* output = textureFilter->GetFeatureMeans();
  for(unsigned int i = 0; i < output->size(); ++i)
    {
      std::cout<<(*output)[i]<<std::endl;
      fprintf(featureFile, "%f \n", (*output)[i]);
    }

  fclose(featureFile);
  return EXIT_SUCCESS;
}
