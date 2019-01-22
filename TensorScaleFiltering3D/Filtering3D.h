#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

struct VECTOR
{
  float x;
  float y;
  float z;
  float length;
};

struct DIFFUSIZE
{
  float upD;
  float upN;
  float upS;
  float upW;
  float upE;
  float upNW;
  float upNE;
  float upSW;
  float upSE;
  float cuN;
  float cuE;
  float cuNW;
  float cuNE;
};

//------------------ITK Part-----------------------------------------------
const unsigned int Dimension = 3;
typedef float PixelType; 
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ImageReaderType;
typedef itk::ImageFileWriter< ImageType > ImageWriterType;
typedef itk::MedianImageFilter< ImageType, ImageType >  FilterType;
typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType >  GaussFilterType;
//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image
int spaRow, spaCol, spaSli; //Spacing

//------------------Coordinate convert-------------------------------------
#define vector1(i,j,k) vector1[(k)*pRow*pCol+(j)*pRow+(i)]
#define vector2(i,j,k) vector2[(k)*pRow*pCol+(j)*pRow+(i)]
#define vector3(i,j,k) vector3[(k)*pRow*pCol+(j)*pRow+(i)]
#define diffuSize(i,j,k) diffuSize[(k)*pRow*pCol+(j)*pRow+(i)]
#define image(i,j,k) image[(k)*pRow*pCol+(j)*pRow+(i)]
#define imageFiltered(i,j,k) imageFiltered[(k)*pRow*pCol+(j)*pRow+(i)]
//------------------Functions ---------------------------------------------
#include "shapeTensorFiltering.cxx"
#include "elliRadius.cxx"

