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
#define fuzzyConImage(i,j,k) fuzzyConImage[(k)*pRow*pCol+(j)*pRow+(i)]
#define fuzzyConFlag(i,j,k) fuzzyConFlag[(k)*pRow*pCol+(j)*pRow+(i)]

//struct for seed and affinity relationship 
struct IndexLabel
{
  unsigned int index;
  unsigned short label;
};

struct IndexKappa
{
  unsigned int index;
  unsigned short kappa;
};

//Function
#include "computeFC.cxx"
