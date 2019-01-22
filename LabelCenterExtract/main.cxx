#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkLabelGeometryImageFilter.h"
 
#include <iostream>
#include "math.h"

int main(int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputFile" << std::endl;
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
 
  typedef itk::LabelGeometryImageFilter< ImageType > LabelGeometryImageFilterType;
  LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
  labelGeometryImageFilter->SetInput( reader->GetOutput() );
 
  // These generate optional outputs.
  labelGeometryImageFilter->CalculatePixelIndicesOn();
  labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();
  labelGeometryImageFilter->CalculateOrientedLabelRegionsOn();
  labelGeometryImageFilter->Update();
 
  LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
  LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
  allLabelsIt = allLabels.begin();
  allLabelsIt++;

  FILE *outFile = fopen(argv[2], "w");

  int ct = 1;
  int i,j,k;
  while(allLabelsIt != allLabels.end())
    { 
      LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;

      i = labelGeometryImageFilter->GetCentroid(labelValue)[0];
      j = labelGeometryImageFilter->GetCentroid(labelValue)[1];
      k = labelGeometryImageFilter->GetCentroid(labelValue)[2];

      std::cout << ct <<" " << i << " " << j << " " << k;
      std::cout << std::endl;

      fprintf(outFile, " %d %d %d\n", i, j, k) ;


      allLabelsIt++;   
      ct++;
    }
 
  return EXIT_SUCCESS;
}
