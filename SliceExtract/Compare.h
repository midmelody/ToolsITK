#include <iostream>
#include <math.h>
#include <string>
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

using namespace std;

//------------------ITK Part-----------------------------------------------
const unsigned int Dimension = 3;
typedef float PixelType; 
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ImageReaderType;

//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image: row, column, slice
int spaRow, spaCol, spaSli; //Spacing: row, column, slice

