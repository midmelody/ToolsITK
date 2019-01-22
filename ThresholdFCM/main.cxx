#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

void calculateCenter(double **degreeMembership, int totalCt, int clusterCt, double *dataPoint, double *clusterCenter, double fuzziness) 
{
  double numerator, denominator;
  double **temp;
  int i, j, k;
  temp = new double*[totalCt];
  for(i=0;i<totalCt;i++)
    {
      temp[i] = new double[clusterCt];
    }
  
  for(i=0;i<totalCt;i++) 
    {
      for(j=0;j<clusterCt;j++) 
	{
	  temp[i][j] = pow(degreeMembership[i][j], fuzziness);
        }
    }

  for(j=0;j<clusterCt;j++)
    {
      numerator = 0.0;
      denominator = 0.0;
      for(i=0;i<totalCt;i++) 
	{
	  numerator +=temp[i][j]*dataPoint[i];
	  denominator += temp[i][j];
	}
      clusterCenter[j] = numerator/denominator;
    }
  
  for(i=0;i<totalCt;++i) 
    {
      delete [] temp[i];
    }
  delete [] temp;
}

double updateDegreeMembership(double **degreeMembership, int totalCt, int clusterCt, double *dataPoint, double *clusterCenter, double fuzziness) 
{
  int i, j, k;
  double t, p, sum;
  double normJ, normK;
  double newUij;
  double maxDiff=0.0, diff;
  for(j=0;j<clusterCt;j++) 
    {
      for(i=0;i<totalCt;i++) 
	{
	  //Get new value
	  sum = 0.0;
	  p = 2 / (fuzziness - 1);
	  for (k=0;k<clusterCt;k++) 
	    {
	      normJ = fabs(dataPoint[i] - clusterCenter[j]);
	      normK = fabs(dataPoint[i] - clusterCenter[k]);			   
	      t = normJ / normK;
	      t = pow(t, p);
	      sum += t;
	    }
	  newUij = 1.0/sum;
	  diff = newUij - degreeMembership[i][j];
	  if (diff > maxDiff)
	    maxDiff = diff;
	  degreeMembership[i][j] = newUij;
        }
    }
  return maxDiff;
}

int comp (const void * elem1, const void * elem2) 
{
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    return (f-s);
}

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile clusterNum epsilon" << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef unsigned char OutPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< OutPixelType, Dimension > OutImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< OutImageType > WriterType;
    
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
   
  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

  //FCM Settings
  int clusterCt = atoi( argv[3] );
  float fuzziness = 2.0;
  float epsilon = atof( argv[4] );

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
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];   


  //Read data points
  int i,j,k;
  double low = 1000000;
  double high = -1000000;
  int totalCt;

  ImageType::IndexType index;
  int t = 0;
  float value;

  /*
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  index[0] = i;
	  index[1] = j;
	  index[2] = k;
	  value = reader->GetOutput()->GetPixel( index );
	  if(value>0)
	    t++;
	}
  totalCt = t;
  */

  totalCt = pRow*pCol*pSli;

  double *dataPoint;
  dataPoint = new double[totalCt];
  double **degreeMembership;
  degreeMembership = new double*[totalCt];
  for(i=0;i<totalCt;i++)
    {
      degreeMembership[i] = new double[clusterCt];
    }
  double *clusterCenter;
  clusterCenter = new double[clusterCt];

  t = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  index[0] = i;
	  index[1] = j;
	  index[2] = k;
	  value = reader->GetOutput()->GetPixel( index );

	  dataPoint[t] = value;
	  if(low>value)
	    low = value;
	  if(high<value)
	    high = value;
	  t++;
	  
	  /*
	  if(value>0)
	    {
	      dataPoint[t] = value;
	      if(low>value)
		low = value;
	      if(high<value)
		high = value;
	      t++;
	    }
	  */
	}

  std::cout<<"Min:"<<low<<" Max:"<<high<<std::endl;

  //Initialize membership as random
  float s, rval;
  int r;
  for(i=0;i<totalCt;i++) 
    {
      s = 0.0;
      r = 100;
      for(j=1;j<clusterCt;j++) 
	{
	  rval = rand() % (r + 1);
	  r -= rval;
	  degreeMembership[i][j] = rval/100.0;
	  s += degreeMembership[i][j];
        }
      degreeMembership[i][0] = 1.0 - s;
    }

  double maxDiff = 10;
  while(maxDiff>epsilon)
    {
      calculateCenter(degreeMembership, totalCt, clusterCt, dataPoint, clusterCenter, fuzziness);
      maxDiff = updateDegreeMembership(degreeMembership, totalCt, clusterCt, dataPoint, clusterCenter, fuzziness);
      std::cout<<"Center: ( ";
      for(i=0;i<clusterCt;i++)
	std::cout<<std::fixed<<std::setprecision(2)<<clusterCenter[i]<<" ";
      std::cout<<"): "<<std::setprecision(7)<<maxDiff<<std::endl;
    }

  std::qsort(clusterCenter, clusterCt, sizeof(double), comp);
  std::cout<<std::endl<<"Final Center:"<<std::endl; 
  for(i=0;i<clusterCt;i++)
    std::cout<<std::fixed<<std::setprecision(2)<<clusterCenter[i]<<" ";
  std::cout<<std::endl;

  delete [] dataPoint;
  for(i=0;i<totalCt;++i) 
    {
      delete [] degreeMembership[i];
    }
  delete [] degreeMembership;

  //Allocate new image
  OutImageType::RegionType region;
  region.SetSize( size );
  OutImageType::Pointer image = OutImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();
  
  double *clusterThre;
  clusterThre = new double[clusterCt-1];
  std::cout<<std::endl<<"Final Thresholds:"<<std::endl;
  for(i=0;i<clusterCt-1;i++)
    {
      clusterThre[i] = (clusterCenter[i]+clusterCenter[i+1])/2;
      std::cout<<std::fixed<<std::setprecision(2)<<clusterThre[i]<<" ";
    }
  std::cout<<std::endl;

  OutImageType::IndexType indexOut;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  index[0] = i;
	  index[1] = j;
	  index[2] = k;
	  value = reader->GetOutput()->GetPixel( index );
	  //std::cout<<value<<std::endl;
	  if(value<clusterThre[0])
	    {
	      indexOut[0] = i;
	      indexOut[1] = j;
	      indexOut[2] = k;	      
	      image->SetPixel(indexOut, 1);
	    }
	  else
	    {
	      t = 0;
	      while((t<clusterCt-1)&&(value>clusterThre[t]))
		t++;
	      t = t+1;
	      indexOut[0] = i;
	      indexOut[1] = j;
	      indexOut[2] = k;
	      image->SetPixel(indexOut, t);		
	    }
	}

  delete [] clusterCenter;
  delete [] clusterThre;

  writer->SetInput( image );
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
