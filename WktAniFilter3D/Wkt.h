 #include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDerivativeImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>
//------------------ITK Part-----------------------------------------------
const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ImageReaderType;
typedef itk::ImageFileWriter< ImageType > ImageWriterType;
typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType >  GaussFilterType;
typedef itk::DerivativeImageFilter< ImageType, ImageType > DevFilterType;

//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image: row, column, slice
int spaRow, spaCol, spaSli; //Spacing: row, column, slice

//------------------Coordinate convert-------------------------------------
#define oriImage(i,j,k) oriImage[(k)*pRow*pCol+(j)*pRow+(i)]
