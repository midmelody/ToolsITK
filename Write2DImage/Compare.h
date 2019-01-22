#include <iostream>
#include <math.h>
#include <string>
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileWriter.h"

using namespace std;

//------------------ITK Part-----------------------------------------------
const unsigned int Dimension = 2;
typedef unsigned char PixelType; 
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileWriter< ImageType > WriterType;

//-----------------Global Variables----------------------------------------
int pRow,pCol;
