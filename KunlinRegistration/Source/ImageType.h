/*
 * =====================================================================================
 *
 *       Filename:  ImageType.h
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  10/25/2009 02:43:17 PM CDT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */

#ifndef IMGAGETYPE_H
#define IMGAGETYPE_H


#include "itkImage.h"

#include "itkRGBPixel.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

const unsigned char Dimension = 3;
typedef float PixelType;
typedef float DispPixelType;
typedef short int HUPixelType;
typedef unsigned char MaskPixelType;

typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::Image< DispPixelType, Dimension > DispImageType;
typedef itk::Image< HUPixelType, Dimension > HUImageType;
typedef itk::Image< MaskPixelType, Dimension > MaskImageType;

typedef itk::Vector< DispPixelType, Dimension > DispVectorType;   
typedef itk::Image< DispVectorType, Dimension > DispFieldType;

typedef itk::RGBPixel< unsigned char > RGBPixelType;
typedef itk::Image< RGBPixelType, Dimension > RGBImageType;
////////////////////////////////////////////////////////////////////////

typedef itk::ImageRegionIterator< ImageType > IteratorType;
typedef itk::ImageRegionIterator< DispImageType > DispIteratorType;
typedef itk::ImageRegionIterator< HUImageType > HUIteratorType;
typedef itk::ImageRegionConstIterator< MaskImageType > MaskIteratorType;


#endif
