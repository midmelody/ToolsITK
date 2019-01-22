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
#include "malloc.h" 

//Coordinate convert
int pRow, pCol, pSli;
#define inputImage(i,j,k) inputImage[(k)*pRow*pCol+(j)*pRow+(i)]
#define affinityMap(i,j,k) affinityMap[(k)*pRow*pCol+(j)*pRow+(i)]
#define fuzzyConImage(i,j,k) fuzzyConImage[(k)*pRow*pCol+(j)*pRow+(i)]
#define labelSegImage(i,j,k) labelSegImage[(k)*pRow*pCol+(j)*pRow+(i)]

//struct for adjacency relationship 
struct AffinityType
{
  unsigned short affinity[27];
};

//struct for affinity relationship 
struct IndexKappa
{
  unsigned int index;
  unsigned short kappa;
};

const int maxAffinity = 10000;

const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;
typedef std::multimap< unsigned short, unsigned int > SeedMapType;
typedef std::set< unsigned short > SeedLabelListType;
typedef std::list< unsigned int > SeedListType;

//Function
#include "computeFC.cxx"
