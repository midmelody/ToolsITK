#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>

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

//-----------------Global Variables----------------------------------------
int pRow,pCol,pSli; //Size of Image: row, column, slice
int spaRow, spaCol, spaSli; //Spacing: row, column, slice

//------------------Coordinate convert-------------------------------------
#define DT(i,j,k) DT[(k)*pRow*pCol+(j)*pRow+(i)]
#define vecOri(i,j,k) vecOri[(k)*pRow*pCol+(j)*pRow+(i)]
#define tangentOri1(i,j,k) tangentOri1[(k)*pRow*pCol+(j)*pRow+(i)]
#define tangentOri2(i,j,k) tangentOri2[(k)*pRow*pCol+(j)*pRow+(i)]
#define smoothTangOri1(i,j,k) smoothTangOri1[(k)*pRow*pCol+(j)*pRow+(i)]
//------------------Functions ---------------------------------------------
#include "funcs.cxx"
#include "FindTangentOri.cxx"


