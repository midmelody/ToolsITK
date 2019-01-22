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

//struct for seed and parameter 
struct IndexLabel
{
  unsigned int index;
  unsigned short label;
};

struct LabelPara
{
  unsigned short label;
  float meanObj;
  float sigmaObj;
  float confidentBackground;
  bool darkObj; //dark obj (true) or bright obj (false)
};

const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;
typedef std::multimap< unsigned int, unsigned int > AdjacencyMapType;

//Function
#include "computeFC.cxx"
