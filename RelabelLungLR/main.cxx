#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkChangeLabelImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"

int main(int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile option(1:Whole L/R; 2: Left U/L; 3: Right U/M/L)" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::LabelGeometryImageFilter< ImageType > LabelGeometryImageFilterType;
  typedef itk::ChangeLabelImageFilter< ImageType, ImageType > ChangeLabelImageFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
  ChangeLabelImageFilterType::Pointer changeLabelImageFilter = ChangeLabelImageFilterType::New();
  WriterType::Pointer writer = WriterType::New();
  
  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  int option = atoi( argv[3] );
   
  //Pipeline
  labelGeometryImageFilter->SetInput( reader->GetOutput() );
  labelGeometryImageFilter->Update();
  //Get all labels
  LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
  LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
  allLabelsIt = allLabels.begin();
  allLabelsIt++;
  //Reassign label
  int ct;
  int labelMap[3];
  int i,j,k;
  int loc[3];
  LabelGeometryImageFilterType::LabelPixelType labelValue;
  if(option == 1)
    {
      //Relabel right:1/left:2 lung, right has smaller X value
      for(ct=0;ct<2;ct++)
	{ 
	  labelValue = *allLabelsIt;
	  i = labelGeometryImageFilter->GetCentroid(labelValue)[0];
	  j = labelGeometryImageFilter->GetCentroid(labelValue)[1];
	  k = labelGeometryImageFilter->GetCentroid(labelValue)[2];
	  std::cout << ct+1 <<" " << i << " " << j << " " << k;
	  std::cout << std::endl;
	  loc[ct] = i;
	  allLabelsIt++;   
	}
      if(loc[0]<loc[1])
	{
	  labelMap[0] = 1;
	  labelMap[1] = 2;
	}
      else
	{
	  labelMap[0] = 2;
	  labelMap[1] = 1;	  
	}
      labelMap[2] = 0;
    }
  else if(option == 2)
    {
      //Relabel left lower:1/upper:2, lower has smaller Z value
      for(ct=0;ct<2;ct++)
	{ 
	  labelValue = *allLabelsIt;
	  i = labelGeometryImageFilter->GetCentroid(labelValue)[0];
	  j = labelGeometryImageFilter->GetCentroid(labelValue)[1];
	  k = labelGeometryImageFilter->GetCentroid(labelValue)[2];
	  std::cout << ct+1 <<" " << i << " " << j << " " << k;
	  std::cout << std::endl;
	  loc[ct] = k;
	  allLabelsIt++;   
	}
      if(loc[0]<loc[1])
	{
	  labelMap[0] = 1;
	  labelMap[1] = 2;
	}
      else
	{
	  labelMap[0] = 2;
	  labelMap[1] = 1;	  
	}
      labelMap[2] = 0;
    }
  else if(option == 3)
    {
      //Relabel right lower:3/middle:4/upper:5, lower has smaller Z value
      for(ct=0;ct<3;ct++)
	{ 
	  labelValue = *allLabelsIt;
	  i = labelGeometryImageFilter->GetCentroid(labelValue)[0];
	  j = labelGeometryImageFilter->GetCentroid(labelValue)[1];
	  k = labelGeometryImageFilter->GetCentroid(labelValue)[2];
	  std::cout << ct+1 <<" " << i << " " << j << " " << k;
	  std::cout << std::endl;
	  loc[ct] = k;
	  allLabelsIt++;   
	}
      if(loc[0]<loc[1])
	{
	  if(loc[1]<loc[2])
	    {
	      labelMap[0] = 1;
	      labelMap[1] = 2;
	      labelMap[2] = 3;
	    }
	  else if(loc[2]<loc[0])
	    {
	      labelMap[0] = 2;
	      labelMap[1] = 3;
	      labelMap[2] = 1;
	    }
	  else
	    {
	      labelMap[0] = 1;
	      labelMap[1] = 3;
	      labelMap[2] = 2;
	    }
	}
      else
	{
	  if(loc[2]>loc[0])
	    {
	      labelMap[0] = 2;
	      labelMap[1] = 1;
	      labelMap[2] = 3;
	    }
	  else if(loc[2]<loc[1])
	    {
	      labelMap[0] = 3;
	      labelMap[1] = 2;
	      labelMap[2] = 1;
	    }
	  else
	    {
	      labelMap[0] = 3;
	      labelMap[1] = 1;
	      labelMap[2] = 2;
	    }
	}
    }
 
  //Get image specs
  ImageType::Pointer image = ImageType::New();
  ImageType::SpacingType spacing =  reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin =  reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction =  reader->GetOutput()->GetDirection();
  ImageType::SizeType  size =  reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = size[0];
  int pCol = size[1];
  int pSli = size[2];
  image->SetRegions( size );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();
  image->FillBuffer( 0 );
  
  int value;
  ImageType::IndexType pixelIndex;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = reader->GetOutput()->GetPixel(pixelIndex);
	  if(value>0)
	    image->SetPixel( pixelIndex, labelMap[value-1] );
	}
  
  writer->SetInput( image );
  writer->Update();
  
  return EXIT_SUCCESS;
}
