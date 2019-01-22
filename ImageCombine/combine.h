#include "itkImage.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Vector< float, 3 > VectorPixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::Image< VectorPixelType, Dimension > OutImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< OutImageType > WriterType;
