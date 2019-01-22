#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileWriter.h"

//------------------ITK Part-----------------------------------------------
const unsigned int Dimension = 3;
typedef float PixelType; 
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileWriter< ImageType > ImageWriterType;

//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image: row, column, slice
int spaRow, spaCol, spaSli; //Spacing: row, column, slice




