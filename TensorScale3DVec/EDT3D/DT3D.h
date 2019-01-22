#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkHConvexImageFilter.h"

struct VOXEL
{
  float x;
  float y;
  float z;
};

//------------------ITK Part-----------------------------------------------
const unsigned int Dimension = 3;
typedef float PixelType; 
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ImageReaderType;
typedef itk::ImageFileWriter< ImageType > ImageWriterType;
typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType, ImageType>  LoGFilterType; 
typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType>  GraFilterType;
typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType> DevFilterType;
typedef itk::HConvexImageFilter<ImageType, ImageType> LocMaxFilterType;

//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image: row, column, slice
float spaRow, spaCol, spaSli; //Spacing: row, column, slice

//------------------Coordinate convert-------------------------------------
#define LoGImageHigh(i,j,k) LoGImageHigh[(k)*pRow*pCol+(j)*pRow+(i)]
#define LoGImageLow(i,j,k) LoGImageLow[(k)*pRow*pCol+(j)*pRow+(i)]
#define DoGImageHigh(i,j,k) DoGImageHigh[(k)*pRow*pCol+(j)*pRow+(i)]
#define DoGImageLow(i,j,k) DoGImageLow[(k)*pRow*pCol+(j)*pRow+(i)]
#define DT(i,j,k) DT[(k)*pRow*pCol+(j)*pRow+(i)]
#define edgeCount(i,j,k) edgeCount[(k)*pRow*pCol+(j)*pRow+(i)]
#define edgePosition(i,j,k) edgePosition[(k)*pRow*pCol+(j)*pRow+(i)]
#define nearestEdgePoint(i,j,k) nearestEdgePoint[(k)*pRow*pCol+(j)*pRow+(i)]
#define NEP(i,j,k) NEP[(k)*pRow*pCol+(j)*pRow+(i)]
#define inputData(i,j,k) inputData[(k)*pRow*pCol+(j)*pRow+(i)]
#define localSize(i,j,k) localSize[(k)*pRow*pCol+(j)*pRow+(i)]
//------------------Functions ---------------------------------------------
#include "funcs.cxx"
#include "ComputeZeroCrossing.cxx"
#include "RasterScanEDT.cxx"

