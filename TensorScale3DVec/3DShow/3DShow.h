#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkPNGImageIO.h"

#include "itkRGBPixel.h"
#include "HSVRGB.cxx"

#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>


struct TENSOR
{
  float x;
  float y;
  float lengthMajor;
  float lengthMinor;
};

//------------------ITK Part-----------------------------------------------
const unsigned int Dimension = 3;
typedef float PixelType; 
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ImageReaderType;
typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType >  GaussFilterType;
typedef itk::ImageFileWriter< ImageType > ImageWriterType;
typedef itk::RGBPixel< unsigned char > RGBPixelType;
typedef itk::Image< RGBPixelType, 3 > RGB3DImageType;
typedef itk::Image< RGBPixelType, 2 > RGB2DImageType;
typedef itk::ImageSeriesWriter< RGB3DImageType, RGB2DImageType > SeriesWriterType;
//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image: row, column, slice
int spaRow, spaCol, spaSli; //Spacing: row, column, slice

//------------------Coordinate convert-------------------------------------
#define tensorXY(i,j,k) tensorXY[(k)*pRow*pCol+(j)*pRow+(i)]
#define tensorYZ(i,j,k) tensorYZ[(k)*pRow*pCol+(j)*pRow+(i)]
#define tensorXZ(i,j,k) tensorXZ[(k)*pRow*pCol+(j)*pRow+(i)]
#define tensorSmtXY(i,j,k) tensorSmtXY[(k)*pRow*pCol+(j)*pRow+(i)]
#define tensorSmtYZ(i,j,k) tensorSmtYZ[(k)*pRow*pCol+(j)*pRow+(i)]
#define tensorSmtXZ(i,j,k) tensorSmtXZ[(k)*pRow*pCol+(j)*pRow+(i)]
//------------------Functions ---------------------------------------------
#include "funcs.cxx"


