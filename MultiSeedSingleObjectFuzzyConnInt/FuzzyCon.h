#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <list>
#include "assert.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

//Coordinate convert
int pRow, pCol, pSli;
#define inputImage(i,j,k) inputImage[(k)*pRow*pCol+(j)*pRow+(i)]
#define affinityMap(i,j,k) affinityMap[(k)*pRow*pCol+(j)*pRow+(i)]
#define fuzzyConImage(i,j,k) fuzzyConImage[(k)*pRow*pCol+(j)*pRow+(i)]
#define labelSegImage(i,j,k) labelSegImage[(k)*pRow*pCol+(j)*pRow+(i)]


const int maxAffinity = 10000;

const unsigned int Dimension = 3;
typedef float OutPixelType;
typedef int PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::Image< OutPixelType, Dimension > OutImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< OutImageType > WriterType;

//Function
#include "computeFC.cxx"
