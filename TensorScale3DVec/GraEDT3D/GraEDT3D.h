#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"

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
typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType >  FilterType;
typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType >  GraMagType;
//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image: row, column, slice
int spaRow, spaCol, spaSli; //Spacing: row, column, slice

//------------------Coordinate convert-------------------------------------
#define DT(i,j,k) DT[(k)*pRow*pCol+(j)*pRow+(i)]
#define GradientDTX(i,j,k) GradientDTX[(k)*pRow*pCol+(j)*pRow+(i)]
#define GradientDTY(i,j,k) GradientDTY[(k)*pRow*pCol+(j)*pRow+(i)]
#define GradientDTZ(i,j,k) GradientDTZ[(k)*pRow*pCol+(j)*pRow+(i)]
#define vecOriGraZenith(i,j,k) vecOriGraZenith[(k)*pRow*pCol+(j)*pRow+(i)]
#define vecOriGraAzimuth(i,j,k) vecOriGraAzimuth[(k)*pRow*pCol+(j)*pRow+(i)]
#define SmoothGraDTX(i,j,k) SmoothGraDTX[(k)*pRow*pCol+(j)*pRow+(i)]
#define SmoothGraDTY(i,j,k) SmoothGraDTY[(k)*pRow*pCol+(j)*pRow+(i)]
#define SmoothGraDTZ(i,j,k) SmoothGraDTZ[(k)*pRow*pCol+(j)*pRow+(i)]

//---------------------Functions --------------------------------------------
#include "funcs.cxx"
