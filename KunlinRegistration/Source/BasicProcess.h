/*
 * =====================================================================================
 *
 *       Filename:  BasicProcess.h
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  04/29/2010 04:38:15 PM CDT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "cubic_BSpline.h"
#include "ImageType.h"
#include "IMGIO.h"

#include "itkCastImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkThresholdImageFilter.h"

#include "itkBinaryBallStructuringElement.h" 
#include "itkBinaryDilateImageFilter.h"

#include "itkRecursiveGaussianImageFilter.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkBSplineResampleImageFilterBase.h"
#include "itkBSplineDownsampleImageFilter.h"
#include "itkBSplineUpsampleImageFilter.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"

#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

int round(float a) 
{
    return int(a + 0.5);
}


int findMaskCT( unsigned char* mask, int imageSize )
{
    int ct=0;
    for (int i=0; i<imageSize; i++)
    {
        if (mask[i]!=0)
            ct++;
    }

    return ct;
}


template < typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer CastImage( const typename TInputImage::Pointer inputImage, bool bFillZero )
{
    typedef itk::CastImageFilter< TInputImage, TOutputImage > CastFilterType;
    typename CastFilterType::Pointer caster = CastFilterType::New();
    caster->SetInput( inputImage );

    try
    {
        caster->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "CastImage exception caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    if ( bFillZero )
    {
        caster->GetOutput()->FillBuffer( 0 );
    }

    return caster->GetOutput();
}


template < typename TImage >
bool CopyImage( const typename TImage::Pointer inputImage, typename TImage::Pointer& outputImage )
{
    typedef itk::ImageRegionIterator< TImage > IteratorType;
    IteratorType it1( inputImage, inputImage->GetRequestedRegion() );
    IteratorType it2( outputImage, outputImage->GetRequestedRegion() );

    for ( it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2 )
    {
        it2.Set( it1.Get() );
    }

    return true;
}

template < typename TImage >
bool NormalizeImage( const typename TImage::Pointer inputImage, typename TImage::Pointer& outputImage, float normalizeFactor )
{
    typedef itk::ImageRegionIterator< TImage > IteratorType;
    IteratorType it1( inputImage, inputImage->GetRequestedRegion() );
    IteratorType it2( outputImage, outputImage->GetRequestedRegion() );

    for ( it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2 )
    {
        it2.Set( it1.Get()*normalizeFactor );
    }

    return true;
}

template < typename TImage >
typename TImage::Pointer ThresholdImage( const typename TImage::Pointer inputImage, const typename TImage::PixelType lowerThreshold, const typename TImage::PixelType upperThreshold, const typename TImage::PixelType insideValue, const typename TImage::PixelType outsideValue )
{
    typedef itk::BinaryThresholdImageFilter< TImage, TImage >  ThresholdFilterType;
    typename ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();
    thresholder->SetInput( inputImage );
    thresholder->SetLowerThreshold( lowerThreshold );
    thresholder->SetUpperThreshold( upperThreshold );
    thresholder->SetInsideValue( insideValue );
    thresholder->SetOutsideValue( outsideValue );

    try
    {
        thresholder->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "ThresholdImage exception caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return thresholder->GetOutput();
}

template < typename TImage, typename TOutputImage >
typename TImage::Pointer ThresholdCastImage( const typename TImage::Pointer inputImage, const typename TImage::PixelType lowerThreshold, const typename TImage::PixelType upperThreshold, const typename TImage::PixelType insideValue, const typename TImage::PixelType outsideValue )
{
    typedef itk::BinaryThresholdImageFilter< TImage, TImage >  ThresholdFilterType;
    typename ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();
    thresholder->SetInput( inputImage );
    thresholder->SetLowerThreshold( lowerThreshold );
    thresholder->SetUpperThreshold( upperThreshold );
    thresholder->SetInsideValue( insideValue );
    thresholder->SetOutsideValue( outsideValue );

    try
    {
        thresholder->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "ThresholdImage exception caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    typedef itk::CastImageFilter< TImage, TOutputImage > CastFilterType;
    typename CastFilterType::Pointer caster = CastFilterType::New();
    caster->SetInput( thresholder->GetOutput() );

    try
    {
        caster->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "CastImage exception caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return caster->GetOutput();
}


template < typename TImage, typename TMask >
typename TImage::Pointer MaskImage( const typename TImage::Pointer image, const typename TMask::Pointer mask, const typename TImage::PixelType bkValue, const typename TMask::PixelType unmarkValue )
{
    if ( image->GetRequestedRegion().GetSize() != mask->GetRequestedRegion().GetSize() )
    {
        std::cout<<"image and mask have different sizes!"<<std::endl;
        return NULL;
    }

    typedef itk::ImageRegionConstIterator< TImage > ConstIteratorType;
    typedef itk::ImageRegionConstIterator< TMask > MaskConstIteratorType;
    ConstIteratorType imageIt( image, image->GetRequestedRegion() );
    MaskConstIteratorType maskIt( mask, mask->GetRequestedRegion() );

    typename TImage::Pointer imageM = CastImage< TImage, TImage >( image, false );
    typedef itk::ImageRegionIterator< TImage > IteratorType;
    IteratorType outIt( imageM, imageM->GetRequestedRegion() );

    for ( outIt.GoToBegin(), maskIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt, ++maskIt )
    {
        if ( maskIt.Get() == unmarkValue )
        {
            outIt.Set( bkValue );
        }
    }

    return imageM;
}


template < typename TImage, typename TMask >
bool MaskImageBool( typename TImage::Pointer& image, const typename TMask::Pointer mask, const typename TImage::PixelType bkValue, const typename TMask::PixelType unmarkValue )
{
    typedef itk::ImageRegionIterator< TImage > IteratorType;
    typedef itk::ImageRegionConstIterator< TMask > MaskIteratorType;
    IteratorType imageIt( image, image->GetRequestedRegion() );
    MaskIteratorType maskIt( mask, mask->GetRequestedRegion() );

    for ( imageIt.GoToBegin(), maskIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, ++maskIt )
    {
        if ( maskIt.Get() == unmarkValue )
        {
            imageIt.Set( bkValue );
        }
    }

    return true;
}

template < typename TImage, typename TMask >
typename TImage::Pointer HistMatching( const typename TImage::Pointer image1, const typename TImage::Pointer image2, const typename TMask::Pointer mask1, const typename TMask::Pointer mask2, const typename TMask::PixelType unmarkValue )
{
    ImageType::Pointer img1 = MaskImage< TImage, TMask >( image1, mask1, -2000, unmarkValue );
    ImageType::Pointer img2 = MaskImage< TImage, TMask >( image2, mask2, -2000, unmarkValue );

    typedef itk::HistogramMatchingImageFilter< TImage, TImage > MatchingFilterType;
    typename MatchingFilterType::Pointer matcher = MatchingFilterType::New();

    matcher->SetInput( img1 );
    matcher->SetReferenceImage( img2 );

    matcher->SetNumberOfHistogramLevels( 1024 );
    matcher->SetNumberOfMatchPoints( 7 );
    matcher->ThresholdAtMeanIntensityOn();

    try
    {
        matcher->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "HistMatching exception caught !" << std::endl;
        std::cerr << err << std::endl;
        exit(2);
    }

    typename TImage::Pointer outImage1 = matcher->GetOutput();
    MaskImageBool< TImage, TMask >( outImage1, mask1, 0, unmarkValue );
    //CopyImage< TImage >( outImage1, image1 );

    //MaskImageBool< TImage, TMask >( image2, mask2, 0, unmarkValue );    
    

    /*typedef itk::ImageRegionIterator< TImage > IteratorType;
    IteratorType it( outImage1, outImage1->GetRequestedRegion() );
    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
        it.Set( static_cast<short int>(it.Get()) );
    }*/


    return outImage1;
}

template < typename TImage, typename TMask >
typename TImage::Pointer HistMatching_bg( const typename TImage::Pointer image1, const typename TImage::Pointer image2, const typename TMask::Pointer mask1, const typename TMask::Pointer mask2, const typename TImage::PixelType bgValue, const typename TMask::PixelType unmarkValue )
{
    ImageType::Pointer img1 = MaskImage< TImage, TMask >( image1, mask1, -2000, unmarkValue );
    ImageType::Pointer img2 = MaskImage< TImage, TMask >( image2, mask2, -2000, unmarkValue );

    typedef itk::HistogramMatchingImageFilter< TImage, TImage > MatchingFilterType;
    typename MatchingFilterType::Pointer matcher = MatchingFilterType::New();

    matcher->SetInput( img1 );
    matcher->SetReferenceImage( img2 );

    matcher->SetNumberOfHistogramLevels( 1024 );
    matcher->SetNumberOfMatchPoints( 7 );
    matcher->ThresholdAtMeanIntensityOn();

    try
    {
        matcher->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "HistMatching exception caught !" << std::endl;
        std::cerr << err << std::endl;
        exit(2);
    }

    typename TImage::Pointer outImage1 = matcher->GetOutput();
    MaskImageBool< TImage, TMask >( outImage1, mask1, bgValue, unmarkValue );
    //CopyImage< TImage >( outImage1, image1 );

    //MaskImageBool< TImage, TMask >( image2, mask2, 0, unmarkValue );    
    

    /*typedef itk::ImageRegionIterator< TImage > IteratorType;
    IteratorType it( outImage1, outImage1->GetRequestedRegion() );
    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
        it.Set( static_cast<short int>(it.Get()) );
    }*/


    return outImage1;
}


template < typename TMask >
typename TMask::Pointer UnionMasks( const typename TMask::Pointer mask1, const typename TMask::Pointer mask2 )
{
    if ( mask1->GetRequestedRegion().GetSize() != mask2->GetRequestedRegion().GetSize() )
    {
        std::cout<<"mask1 and mask2 have different sizes!"<<std::endl;
        return NULL;
    }

    typedef itk::ImageRegionConstIterator< TMask > MaskConstIteratorType;
    MaskConstIteratorType mask1It( mask1, mask1->GetRequestedRegion() );
    MaskConstIteratorType mask2It( mask2, mask2->GetRequestedRegion() );

    typename TMask::Pointer maskUnion = CastImage< TMask, TMask >( mask1, true );
    typedef itk::ImageRegionIterator< TMask > MaskIteratorType;
    MaskIteratorType maskUnionIt( maskUnion, maskUnion->GetRequestedRegion() );

    std::cout<<"Union Mask Images..."<<std::endl;
    typedef typename TMask::PixelType PixelType;
    PixelType PixelValue;
    int counter = 0;
    for ( mask1It.GoToBegin(), mask2It.GoToBegin(), maskUnionIt.GoToBegin(); !maskUnionIt.IsAtEnd(); ++mask1It, ++mask2It, ++maskUnionIt )
    {
        PixelValue = static_cast<PixelType>( (mask1It.Get()>0) || (mask2It.Get()>0) );
        maskUnionIt.Set( PixelValue );

        if ( PixelValue > 0 )
        {
            counter++;
        }

    }

    std::cout<<"Mask Union Points: "<<counter<<std::endl;
    return maskUnion;
}


void union_masks( unsigned char* mask1, unsigned char* mask2, unsigned char* maskUnion, int imgSize )
{
    for (int index=0; index<imgSize; index++)
    {
        if (mask1[index]>0 || mask2[index]>0)
            maskUnion[index]=1;
        else
            maskUnion[index]=0;
    }

}


template < typename TImage >
bool DiffImages( const typename TImage::Pointer image1, const typename TImage::Pointer image2 )
{
    if ( image1->GetRequestedRegion().GetSize() != image2->GetRequestedRegion().GetSize() )
    {
        std::cout<<"image1 and image2 have different sizes!"<<std::endl;
        return NULL;
    }

    typedef itk::ImageRegionConstIterator< TImage > ConstIteratorType;
    ConstIteratorType image1It( image1, image1->GetRequestedRegion() );
    ConstIteratorType image2It( image2, image2->GetRequestedRegion() );

    std::cout<<"Comparing Images..."<<std::endl;
    for ( image1It.GoToBegin(), image2It.GoToBegin(); !image1It.IsAtEnd(); ++image1It, ++image2It )
    {
        if ( image1It.Get() != image2It.Get() )
        {
            std::cout<<"Different!"<<std::endl;
            return false;
        }
    }

    std::cout<<"Same."<<std::endl;
    return true;
}


template < typename TImage >
bool ITKImage2Array( const typename TImage::Pointer image, float* u )
{
    typename TImage::IndexType PixelIndex;
    typename TImage::SizeType size = image->GetLargestPossibleRegion().GetSize();
    int xdim = size[0];
    int ydim = size[1];
    int zdim = size[2];

    int index=0;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                PixelIndex[0] = x;
                PixelIndex[1] = y;
                PixelIndex[2] = z;
                u[index] = static_cast<float>( image->GetPixel( PixelIndex ) );
                index++;
            }
        }
    }

    int imageSize = xdim*ydim*zdim;
    if ( index == imageSize )
    {
        return true;
    }
    else
    {
        std::cout<<"Errors when converting ITKImage to Array!"<<std::endl;
        std::cout<<"ITKImage size = "<<imageSize<<"; Write array space = "<<index<<std::endl;
        return false;
    }

}


template < typename TImage >
bool ITKImage2ArrayUCHAR( const typename TImage::Pointer image, unsigned char* u )
{
    typename TImage::IndexType PixelIndex;
    typename TImage::SizeType size = image->GetLargestPossibleRegion().GetSize();
    int xdim = size[0];
    int ydim = size[1];
    int zdim = size[2];

    int index=0;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                PixelIndex[0] = x;
                PixelIndex[1] = y;
                PixelIndex[2] = z;
                u[index] = static_cast<unsigned char>( image->GetPixel( PixelIndex ) );
                index++;
            }
        }
    }

    int imageSize = xdim*ydim*zdim;
    if ( index == imageSize )
    {
        return true;
    }
    else
    {
        std::cout<<"Errors when converting ITKImage to ArrayUCHAR!"<<std::endl;
        std::cout<<"ITKImage size = "<<imageSize<<"; Write array space = "<<index<<std::endl;
        return false;
    }

}


template < typename TImage >
bool Array2ITKImage( float* u, typename TImage::Pointer& image )
{
    typedef typename TImage::PixelType PixelType;
    typename TImage::IndexType PixelIndex;
    typename TImage::SizeType size = image->GetLargestPossibleRegion().GetSize();
    int xdim = size[0];
    int ydim = size[1];
    int zdim = size[2];

    int index=0;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                PixelIndex[0] = x;
                PixelIndex[1] = y;
                PixelIndex[2] = z;
                image->SetPixel( PixelIndex, static_cast<PixelType>( u[index] ) ); 
                index++;
            }
        }
    }

    int imageSize = xdim*ydim*zdim;
    if ( index == imageSize )
    {
        return true;
    }
    else
    {
        std::cout<<"Errors when converting Array to ITKImage!"<<std::endl;
        std::cout<<"ITKImage size = "<<imageSize<<"; Read array space = "<<index<<std::endl;
        return false;
    }

}

template < typename TImage >
bool ArrayUCHAR2ITKImage( unsigned char* u, typename TImage::Pointer& image )
{
    typedef typename TImage::PixelType PixelType;
    typename TImage::IndexType PixelIndex;
    typename TImage::SizeType size = image->GetLargestPossibleRegion().GetSize();
    int xdim = size[0];
    int ydim = size[1];
    int zdim = size[2];

    int index=0;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                PixelIndex[0] = x;
                PixelIndex[1] = y;
                PixelIndex[2] = z;
                image->SetPixel( PixelIndex, static_cast<PixelType>( u[index] ) ); 
                index++;
            }
        }
    }

    int imageSize = xdim*ydim*zdim;
    if ( index == imageSize )
    {
        return true;
    }
    else
    {
        std::cout<<"Errors when converting ArrayUCHAR to ITKImage!"<<std::endl;
        std::cout<<"ITKImage size = "<<imageSize<<"; Read array space = "<<index<<std::endl;
        return false;
    }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename TImage >
bool DispITKImage2Array( const typename TImage::Pointer image1, const typename TImage::Pointer image2, const typename TImage::Pointer image3, float* u1, float* u2, float* u3 )
{
    typename TImage::IndexType PixelIndex;
    typename TImage::SizeType size = image1->GetLargestPossibleRegion().GetSize();
    int xdim = size[0];
    int ydim = size[1];
    int zdim = size[2];

    int index=0;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                PixelIndex[0] = x;
                PixelIndex[1] = y;
                PixelIndex[2] = z;
                u1[index] = static_cast<float>( image1->GetPixel( PixelIndex ) );
                u2[index] = static_cast<float>( image2->GetPixel( PixelIndex ) );
                u3[index] = static_cast<float>( image3->GetPixel( PixelIndex ) );
                index++;
            }
        }
    }

    int imageSize = xdim*ydim*zdim;
    if ( index == imageSize )
    {
        return true;
    }
    else
    {
        std::cout<<"Errors when converting DispITKImages to Arrays!"<<std::endl;
        std::cout<<"ITKImage size = "<<imageSize<<"; Write array space = "<<index<<std::endl;
        return false;
    }

}

template < typename TImage >
bool DispArray2ITKImage( float* u1, float* u2, float* u3, typename TImage::Pointer& image1, typename TImage::Pointer& image2, typename TImage::Pointer& image3 )
{
    typedef typename TImage::PixelType PixelType;
    typename TImage::IndexType PixelIndex;
    typename TImage::SizeType size = image1->GetLargestPossibleRegion().GetSize();
    int xdim = size[0];
    int ydim = size[1];
    int zdim = size[2];

    int index=0;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                PixelIndex[0] = x;
                PixelIndex[1] = y;
                PixelIndex[2] = z;
                image1->SetPixel( PixelIndex, static_cast<PixelType>( u1[index] ) ); 
                image2->SetPixel( PixelIndex, static_cast<PixelType>( u2[index] ) ); 
                image3->SetPixel( PixelIndex, static_cast<PixelType>( u3[index] ) ); 
                index++;
            }
        }
    }

    int imageSize = xdim*ydim*zdim;
    if ( index == imageSize )
    {
        return true;
    }
    else
    {
        std::cout<<"Errors when converting DispArrays to ITKImages!"<<std::endl;
        std::cout<<"ITKImage size = "<<imageSize<<"; Read array space = "<<index<<std::endl;
        return false;
    }

}


template < typename TImage >
typename TImage::Pointer DilateImage( const typename TImage::Pointer inputImage, const unsigned char dilateRadius, const typename TImage::PixelType dilateValue )
{
    typedef itk::BinaryBallStructuringElement< typename TImage::PixelType, 3 >  StructuringElementType;
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( dilateRadius );  // (2n+1)x(2n+1) structuring element
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter< TImage, TImage, StructuringElementType >  DilateFilterType;
    typename DilateFilterType::Pointer dilater = DilateFilterType::New();
    dilater->SetKernel( structuringElement );
    dilater->SetInput( inputImage );
    dilater->SetDilateValue( dilateValue );

    try
    {
        dilater->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "DilateImage exception caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return dilater->GetOutput();
}






//////////////////////////////////////////////////////////////////////////////////////////////////////

template < typename TInputImage, typename TOutputImage, typename TMask >
typename TOutputImage::Pointer CalcVesselnessMeasure( const typename TInputImage::Pointer HUImage, const typename TMask::Pointer mask, float sigmaMin, float sigmaMax, int nScales )
{
    typedef double     OutputVesselnessPixelType;
    typedef itk::Image< OutputVesselnessPixelType, 3 > VesselnessOutputImageType;
    
	typedef itk::HessianToObjectnessMeasureImageFilter< double,3 > VesselnessFilterType;

    // Declare the type of multiscale vesselness filter
    typedef itk::MultiScaleHessianBasedMeasureImageFilter<
                                            TInputImage, VesselnessFilterType,
                                            VesselnessOutputImageType>  
                                            MultiScaleVesselnessFilterType;

    // Create a vesselness Filter
    typename MultiScaleVesselnessFilterType::Pointer MultiScaleVesselnessFilter = 
                                      MultiScaleVesselnessFilterType::New();
    MultiScaleVesselnessFilter->SetInput( HUImage );
    MultiScaleVesselnessFilter->SetSigmaMin( sigmaMin ); 
    MultiScaleVesselnessFilter->SetSigmaMax( sigmaMax ); 
    MultiScaleVesselnessFilter->SetNumberOfSigmaSteps( nScales ); 

	//MultiScaleVesselnessFilter->PrintSelf(std::cout,0);
    try
    {
    MultiScaleVesselnessFilter->Update();
    }
    catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return NULL;
    }

    VesselnessOutputImageType::Pointer vmOut = MultiScaleVesselnessFilter->GetOutput();
    VesselnessOutputImageType::Pointer vmMasked = MaskImage< VesselnessOutputImageType, TMask >( vmOut, mask, 0, 0 );

    typedef itk::RescaleIntensityImageFilter< VesselnessOutputImageType, TOutputImage > RescaleFilterType;
    typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetInput( vmMasked );
    //rescaler->SetInput( vmOut );
    rescaler->SetOutputMinimum( 0 );
    rescaler->SetOutputMaximum( 1 );
    rescaler->Update();

    return rescaler->GetOutput();
}


void rescaleImageIntensity( float* image, float outMin, float outMax, int xdim, int ydim, int zdim )
{
    int index;
    float inMin=1000;
    float inMax=-1000;
    float tmp;
    int planeSize=xdim*ydim;

    int x, y, z;
    for (z=0; z<zdim; z++)
        for (y=0; y<ydim; y++)
            for(x=0; x<xdim; x++)
            {
                index=z*planeSize+y*xdim+x;
                tmp = image[index];
                if ( tmp > inMax )
                    inMax = tmp;
                if ( tmp < inMin )
                    inMin = tmp;
            }


    float scale = (outMax-outMin)/(inMax-inMin);
    for (z=0; z<zdim; z++)
        for (y=0; y<ydim; y++)
            for(x=0; x<xdim; x++)
            {
                index=z*planeSize+y*xdim+x;
                tmp = image[index];
                
                image[index] = (tmp-inMin)*scale+outMin;
            }
}


template <typename TDispField, typename  TInputImage>
bool ComposeDispField( typename TDispField::Pointer& dfield, const typename TInputImage::Pointer inputimage0, const typename TInputImage::Pointer inputimage1, const typename TInputImage::Pointer inputimage2)
{
	typename TInputImage::SizeType size = inputimage0->GetLargestPossibleRegion().GetSize();
	const int xsize = size[0];
	const int ysize = size[1];
	const int zsize = size[2];

	typename TInputImage::SpacingType spacing = inputimage0->GetSpacing();
	typename TInputImage::IndexType inputindex;
	typename TInputImage::PixelType inputvalue0, inputvalue1, inputvalue2;

	typename TDispField::IndexType outputindex;
	typename TDispField::PixelType outputvalue;

	for(int k=0; k<zsize; k++) //slice
    {
		for(int j=0; j<ysize; j++)//row
        {
			for(int i=0; i<xsize; i++)//col
			{
				inputindex[0]=i;
				inputindex[1]=j;
				inputindex[2]=k;

				inputvalue0 = inputimage0->GetPixel(inputindex);
				inputvalue1 = inputimage1->GetPixel(inputindex);
				inputvalue2 = inputimage2->GetPixel(inputindex);

				//outputvalue0 = inputvalue[0]/spacing[0];
				//outputvalue1 = inputvalue[1]/spacing[1];
				//outputvalue2 = inputvalue[2]/spacing[2];

                outputvalue[0] = inputvalue0;
				outputvalue[1] = inputvalue1;
				outputvalue[2] = inputvalue2;

				outputindex[0]=i;
				outputindex[1]=j;
				outputindex[2]=k;

                dfield->SetPixel( outputindex, outputvalue );
			}
        }
    }
	
    return true;
}

template <typename TDispField, typename  TOutputImage>
bool DecomposeDispField( const typename TDispField::Pointer dfield, typename TOutputImage::Pointer &outputimage0,  typename TOutputImage::Pointer &outputimage1,  typename TOutputImage::Pointer  &outputimage2)
{
	typename TDispField::SizeType size = dfield->GetLargestPossibleRegion().GetSize();
	const int xsize = size[0];
	const int ysize = size[1];
	const int zsize = size[2];

	typename TDispField::SpacingType spacing = dfield->GetSpacing();
	typename TDispField::IndexType inputindex;
	typename TDispField::PixelType inputvalue;

	typename TOutputImage::IndexType outputindex;
	typename TOutputImage::PixelType outputvalue0, outputvalue1, outputvalue2;

	for(int k=0; k<zsize; k++) //slice
    {
		for(int j=0; j<ysize; j++)//row
        {
			for(int i=0; i<xsize; i++)//col
			{
				inputindex[0]=i;
				inputindex[1]=j;
				inputindex[2]=k;

				inputvalue = dfield->GetPixel(inputindex);

				//outputvalue0 = inputvalue[0]/spacing[0];
				//outputvalue1 = inputvalue[1]/spacing[1];
				//outputvalue2 = inputvalue[2]/spacing[2];

                outputvalue0 = inputvalue[0];
				outputvalue1 = inputvalue[1];
				outputvalue2 = inputvalue[2];

				outputindex[0]=i;
				outputindex[1]=j;
				outputindex[2]=k;

				outputimage0->SetPixel( outputindex,outputvalue0);
				outputimage1->SetPixel( outputindex,outputvalue1);
				outputimage2->SetPixel( outputindex,outputvalue2);
			}
        }
    }
	
    return true;
}

template < typename TDispField >
bool DispField2Arrays( const typename TDispField::Pointer dfield, float* u1, float* u2, float* u3 )
{
    typename TDispField::SizeType size = dfield->GetLargestPossibleRegion().GetSize();
    const int xdim = size[0];
    const int ydim = size[1];
    const int zdim = size[2];

	typename TDispField::SpacingType spacing = dfield->GetSpacing();
	typename TDispField::IndexType inputIndex;
	typename TDispField::PixelType inputValue;

    int index=0;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                inputIndex[0] = x;
                inputIndex[1] = y;
                inputIndex[2] = z;
                inputValue = dfield->GetPixel(inputIndex);

                //u1[index] = static_cast<float>( inputValue[0]/spacing[0] );
                //u2[index] = static_cast<float>( inputValue[1]/spacing[1] );
                //u3[index] = static_cast<float>( inputValue[2]/spacing[2] );

                u1[index] = static_cast<float>( inputValue[0] );
                u2[index] = static_cast<float>( inputValue[1] );
                u3[index] = static_cast<float>( inputValue[2] );
                index++;
            }
        }
    }

    int imageSize = xdim*ydim*zdim;
    if ( index == imageSize )
    {
        return true;
    }
    else
    {
        std::cout<<"Errors when converting DispITKImages to Arrays!"<<std::endl;
        std::cout<<"ITKImage size = "<<imageSize<<"; Write array space = "<<index<<std::endl;
        return false;
    }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//rescale image
void downsizeImg_byAvg_OuterZero(float* u, float* u_fullsize, int xdimo, int ydimo, int zdimo, int scale_x, int scale_y, int scale_z, int scale_u, float zero)
{
	int xdim=(xdimo-1)/scale_x+1;
	int ydim=(ydimo-1)/scale_y+1;
	int zdim=(zdimo-1)/scale_z+1;
    int planeSize=xdimo*ydimo;

	int scale_total=scale_x*scale_y*scale_z*scale_u;
    int index=0;
    int id;
	for (int z=0;z<zdim;z++)
	{
		int z_scaled=scale_z*z;
		for (int y=0;y<ydim;y++)
		{
			int y_scaled=scale_y*y;
			for (int x=0;x<xdim;x++)
			{
				int x_scaled=scale_x*x;
				u[index]=0;
				for (int dz=0;dz<scale_z;dz++)
				{
					for (int dy=0;dy<scale_y;dy++)
					{
						for (int dx=0;dx<scale_x;dx++)
						{
							if (x_scaled+dx>=xdimo || y_scaled+dy>=ydimo || z_scaled+dz>=zdimo)
								u[index]+=zero;
                            else
                            {
                                id=(z_scaled+dz)*planeSize+(y_scaled+dy)*xdimo+(x_scaled+dx);
							    u[index]+=u_fullsize[id];
                            }
						}
					}
				}
                index++;

			}
		}
	}

	for (index=0; index<xdim*ydim*zdim; index++)
	{
        u[index]/=scale_total;
    }
}

//rescale disp image
void downsizeImg_byAvg_bdry(float* u, float* u_fullsize, int xdimo, int ydimo, int zdimo, int scale_x, int scale_y, int scale_z, int scale_u)
{
    int fx,fy,fz;
	int xdim=(xdimo-1)/scale_x+1;
	int ydim=(ydimo-1)/scale_y+1;
	int zdim=(zdimo-1)/scale_z+1;
    int planeSize=xdimo*ydimo;

	int scale_total=scale_x*scale_y*scale_z*scale_u;
    int index=0;
    int id;
	for (int z=0;z<zdim;z++)
	{
		int z_scaled=scale_z*z;
		for (int y=0;y<ydim;y++)
		{
			int y_scaled=scale_y*y;
			for (int x=0;x<xdim;x++)
			{
				int x_scaled=scale_x*x;
				u[index]=0;

				for (int dz=0;dz<scale_z;dz++)
				{
                    fz=z_scaled+dz;
                    if (fz>=zdimo)
                        fz=zdimo;
					for (int dy=0;dy<scale_y;dy++)
					{
                        fy=y_scaled+dy;
                        if (fy>=ydimo)
                            fy=ydimo;
						for (int dx=0;dx<scale_x;dx++)
						{
                            fx=x_scaled+dx;
                            if (fx>=xdimo)
                                fx=xdimo;
                            
                            id=fz*planeSize+fy*xdimo+fx;
							u[index]+=u_fullsize[id];
						}
					}
				}
                index++;

			}
		}
	}

	for (index=0; index<xdim*ydim*zdim; index++)
	{
        u[index]/=scale_total;
    }
}

//rescale mask
void downsizeMask(unsigned char* u, unsigned char* u_fullsize, int xdimo, int ydimo, int zdimo, int scale_x, int scale_y, int scale_z)
{
	int xdim=(xdimo-1)/scale_x+1;
	int ydim=(ydimo-1)/scale_y+1;
	int zdim=(zdimo-1)/scale_z+1;
	int planeSize=xdimo*ydimo;

    int scale_total=scale_x*scale_y*scale_z;
    int index=0;
    int id;
	for (int z=0;z<zdim;z++)
	{
		int z_scaled=scale_z*z;
		for (int y=0;y<ydim;y++)
		{
			int y_scaled=scale_y*y;
			for (int x=0;x<xdim;x++)
			{
				int x_scaled=scale_x*x;

				float temp=0;
				for (int dz=0;dz<scale_z;dz++)
				{
					if (z_scaled+dz>=zdimo)
						continue;
					for (int dy=0;dy<scale_y;dy++)
					{
						if (y_scaled+dy>=ydimo)
							continue;
						for (int dx=0;dx<scale_x;dx++)
						{
							if (x_scaled+dx>=xdimo)
								continue;
                            
                            id=(z_scaled+dz)*planeSize+(y_scaled+dy)*xdimo+(x_scaled+dx);
							temp+=u_fullsize[id];
						}
					}
				}
                if (temp/scale_total>=0.5)
                    u[index]=1;
                else
                    u[index]=0;

                index++;

			}
		}
	}

}

void downsizeMask_wCt(unsigned char* u, unsigned char* u_fullsize, int xdimo, int ydimo, int zdimo, int scale_x, int scale_y, int scale_z, int& maskCT)
{
	int xdim=(xdimo-1)/scale_x+1;
	int ydim=(ydimo-1)/scale_y+1;
	int zdim=(zdimo-1)/scale_z+1;
    int planeSize=xdimo*ydimo;

	int scale_total=scale_x*scale_y*scale_z;
    int index=0;
    int id;
	for (int z=0;z<zdim;z++)
	{
		int z_scaled=scale_z*z;
		for (int y=0;y<ydim;y++)
		{
			int y_scaled=scale_y*y;
			for (int x=0;x<xdim;x++)
			{
				int x_scaled=scale_x*x;
				float temp=0;
				for (int dz=0;dz<scale_z;dz++)
				{
					if (z_scaled+dz>=zdimo)
						continue;
					for (int dy=0;dy<scale_y;dy++)
					{
						if (y_scaled+dy>=ydimo)
							continue;
						for (int dx=0;dx<scale_x;dx++)
						{
							if (x_scaled+dx>=xdimo)
								continue;

                            id=(z_scaled+dz)*planeSize+(y_scaled+dy)*xdimo+(x_scaled+dx);
							temp+=u_fullsize[id];
						}
					}
				}
                if (temp/scale_total>=0.5)
                {
                    u[index]=1;
                    maskCT++;
                }
                else
                    u[index]=0;

                index++;

			}
		}
	}

}

int count_MaskPnt( unsigned char* mask, unsigned char markValue, int imgSize )
{
    int maskCT=0;
    for (int i=0; i<imgSize; i++)
    {
        if (mask[i] == markValue)
            maskCT++;
    }

    return maskCT;
}

template < typename TImage >
int count_MaskPnt( const typename TImage::Pointer image, const typename TImage::PixelType bgValue )
{
    int fgCT=0;

    typedef itk::ImageRegionIterator< TImage > IteratorType;
    IteratorType imageIt( image, image->GetRequestedRegion() );
    for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt)
    {
        if (imageIt.Get() != bgValue)
            fgCT++;
    }

    return fgCT;
}


void expand_mask( unsigned char* maskExpand, unsigned char* mask, int xdim, int ydim, int zdim, int& maskCT2 )
{
    int imgSize=xdim*ydim*zdim;
    int planeSize=xdim*ydim;
    int id0, id1, id2, id3, id4, id5, id6;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                id0=z*planeSize+y*xdim+x;
                id1=z*planeSize+y*xdim+x-1;
                id2=z*planeSize+y*xdim+x+1;
                id3=z*planeSize+(y-1)*xdim+x;
                id4=z*planeSize+(y+1)*xdim+x;
                id5=(z-1)*planeSize+y*xdim+x;
                id6=(z+1)*planeSize+y*xdim+x;

                int sum=mask[id0];
                if (id1>=0)
                    sum+=mask[id1];
                if (id3>=0)
                    sum+=mask[id3];
                if (id5>=0)
                    sum+=mask[id5];
                if (id2<imgSize)
                    sum+=mask[id2];
                if (id4<imgSize)
                    sum+=mask[id4];
                if (id6<imgSize)
                    sum+=mask[id6];
                if (sum>0)
                {
                    maskExpand[id0]=1;
                    maskCT2++;
                }
                else
                {
                    maskExpand[id0]=0;
                }

            }
        }
    }

}

void expand_image( float* imageExpand, float* image, int xdim, int ydim, int zdim, float zero)
{
    //int imageSize=xdim*ydim*zdim;
    int planeSize=(xdim+1)*(ydim+1);
    //int index, index0;
    int id_diff=planeSize+xdim+2;
    //for (index=0; index<imageSize; index++)
      //  imageExpand[index]=zero;
    int index=0;
    for (int z=0; z<zdim+2; z++)
    {
        for (int y=0; y<ydim+2; y++)
        {
            for (int x=0; x<xdim+2; x++)
            {
                //index=z*planeSize+y*xdim+x;
                if (x==0 || x==xdim+1 || y==0 || y==ydim+1 || z==0 || z==zdim+1)
                    imageExpand[index]=zero;
                else
                {
                    //index0=(z-1)*planeSize+(y-1)*(xdim+1)+(x-1);
                    //index0=index-planeSize-(xdim+1)-1;
                    //index0=index-id_diff;
                    imageExpand[index]=image[index-id_diff];
                }

                index++;

            }
        }
    }

}

void expand_disp( float* dispExpand, float* disp, int xdim, int ydim, int zdim )
{
    int planeSize=(xdim+1)*(ydim+1);
    int id_diff=planeSize+xdim+2;
    int index=0;
    for (int z=0; z<zdim+2; z++)
    {
        if (z==0)
            z=1;
        else if (z==zdim+1)
            z=zdim;

        for (int y=0; y<ydim+2; y++)
        {
            if (y==0)
                y=1;
            else if (y==ydim+1)
                y=ydim;

            for (int x=0; x<xdim+2; x++)
            {
                if (x==0)
                    x=1;
                else if (x==xdim+1)
                    x=xdim;

                dispExpand[index]=disp[index-id_diff];

                index++;

            }
        }
    }

}

void recover_img( float* img, float* imgExpand, int xdim, int ydim, int zdim)
{
    int imgSize=xdim*ydim*zdim;
    int planeSize=(xdim+1)*(ydim+1);
    int id_diff=planeSize+xdim+2;
    for (int index=0; index<imgSize; index++)
        img[index]=imgExpand[index+id_diff];

    /*int index=0;
    for (int z=0; z<zdim; z++)
    {
        for (int y=0; y<ydim; y++)
        {
            for (int x=0; x<xdim; x++)
            {
                img[index]=imgExpand[index+id_diff];
                index++;
            }
        }
    }*/

}

void generate_array( float* arrayBasis0, float* arrayBasis1, float* arrayBasis2, float* arrayBasis3, int* arrayLocation, float* u_x_existing, float* u_y_existing, float* u_z_existing, int xdim, int ydim, int zdim, int grid_space_x, int grid_space_y, int grid_space_z, unsigned char* mask )
{
    int x, y, z;
    float xe, ye, ze;
    float xBasis, yBasis, zBasis;
    int i0, j0, k0;
    int i, j, k;
    int index, id;
    int startBasis, startLocation;
    int maskID=0;

    int planeSize=xdim*ydim;
    index=0;
    for (z=0; z<zdim; z++)
    {
        for (y=0; y<ydim; y++)
        {
            for(x=0; x<xdim; x++)
            {
                index=z*planeSize+y*xdim+x;
                if (mask[index]>0)
                {
                    startBasis=maskID*64;
                    startLocation=maskID*4;

                    xe=x+u_x_existing[index];
                    ye=y+u_y_existing[index];
                    ze=z+u_z_existing[index];

                    xBasis=xe/grid_space_x;
                    yBasis=ye/grid_space_y;
                    zBasis=ze/grid_space_z;

                    i0=static_cast<int>(xBasis)-1;
                    if (xBasis<0)
                    {
                        i0-=1;
                    }

                    j0=static_cast<int>(yBasis)-1;
                    if (yBasis<0)
                    {
                        j0-=1;
                    }

                    k0=static_cast<int>(zBasis)-1;
                    if (zBasis<0)
                    {
                        k0-=1;
                    }

                    id=0;
                    for (k=k0; k<k0+4; k++)
                    {
                        for (j=j0; j<j0+4; j++)
                        {
                            for (i=i0; i<i0+4; i++)
                            {
                                arrayBasis0[startBasis+id]=cubic_bsp(xBasis-i)*cubic_bsp(yBasis-j)*cubic_bsp(zBasis-k);
                                arrayBasis1[startBasis+id]=cubic_dbsp(xBasis-i)*cubic_bsp(yBasis-j)*cubic_bsp(zBasis-k)/grid_space_x;
                                arrayBasis2[startBasis+id]=cubic_bsp(xBasis-i)*cubic_dbsp(yBasis-j)*cubic_bsp(zBasis-k)/grid_space_y;
                                arrayBasis3[startBasis+id]=cubic_bsp(xBasis-i)*cubic_bsp(yBasis-j)*cubic_dbsp(zBasis-k)/grid_space_z;

                                id++;
                            }
                        }
                    }

                    arrayLocation[startLocation]=i0;
                    arrayLocation[startLocation+1]=j0;
                    arrayLocation[startLocation+2]=k0;
                    arrayLocation[startLocation+3]=index;
                    //arrayLocation[start+4]=x;
                    //arrayLocation[start+5]=y;
                    //arrayLocation[start+6]=z;
                
                    maskID++;
                }

                //index++;
            }
        }
    }

}

void generate_arrayd( float* arrayBasis0, float* arrayBasis1, float* arrayBasis2, float* arrayBasis3, int* arrayLocation, float* u_x_existing, float* u_y_existing, float* u_z_existing, int xdim, int ydim, int zdim, int grid_space_x, int grid_space_y, int grid_space_z, unsigned char* mask )
{
    int x, y, z;
    double xe, ye, ze;
    double xBasis, yBasis, zBasis;
    int i0, j0, k0;
    int i, j, k;
    int index, id;
    int startBasis, startLocation;
    int maskID=0;

    int planeSize=xdim*ydim;
    index=0;
    for (z=0; z<zdim; z++)
    {
        for (y=0; y<ydim; y++)
        {
            for(x=0; x<xdim; x++)
            {
                index=z*planeSize+y*xdim+x;
                if (mask[index]>0)
                {
                    startBasis=maskID*64;
                    startLocation=maskID*4;

                    xe=static_cast<double>(x+u_x_existing[index]);
                    ye=static_cast<double>(y+u_y_existing[index]);
                    ze=static_cast<double>(z+u_z_existing[index]);

                    xBasis=xe/grid_space_x;
                    yBasis=ye/grid_space_y;
                    zBasis=ze/grid_space_z;

                    i0=static_cast<int>(xBasis)-1;
                    if (xBasis<0)
                    {
                        i0-=1;
                    }

                    j0=static_cast<int>(yBasis)-1;
                    if (yBasis<0)
                    {
                        j0-=1;
                    }

                    k0=static_cast<int>(zBasis)-1;
                    if (zBasis<0)
                    {
                        k0-=1;
                    }

                    id=0;
                    for (k=k0; k<k0+4; k++)
                    {
                        for (j=j0; j<j0+4; j++)
                        {
                            for (i=i0; i<i0+4; i++)
                            {
                                arrayBasis0[startBasis+id]=static_cast<float>(dcubic_bsp(xBasis-i)*dcubic_bsp(yBasis-j)*dcubic_bsp(zBasis-k));
                                arrayBasis1[startBasis+id]=static_cast<float>(dcubic_dbsp(xBasis-i)*dcubic_bsp(yBasis-j)*dcubic_bsp(zBasis-k)/grid_space_x);
                                arrayBasis2[startBasis+id]=static_cast<float>(dcubic_bsp(xBasis-i)*dcubic_dbsp(yBasis-j)*dcubic_bsp(zBasis-k)/grid_space_y);
                                arrayBasis3[startBasis+id]=static_cast<float>(dcubic_bsp(xBasis-i)*dcubic_bsp(yBasis-j)*dcubic_dbsp(zBasis-k)/grid_space_z);

                                id++;
                            }
                        }
                    }

                    arrayLocation[startLocation]=i0;
                    arrayLocation[startLocation+1]=j0;
                    arrayLocation[startLocation+2]=k0;
                    arrayLocation[startLocation+3]=index;
                    //arrayLocation[start+4]=x;
                    //arrayLocation[start+5]=y;
                    //arrayLocation[start+6]=z;
                
                    maskID++;
                }

                //index++;
            }
        }
    }

}

void generate_arrayd2( float* arrayBasis0, float* arrayBasis1, float* arrayBasis2, float* arrayBasis3,  float* arrayBasis4, float* arrayBasis5, float* arrayBasis6, int* arrayLocation, float* u_x_existing, float* u_y_existing, float* u_z_existing, int xdim, int ydim, int zdim, int grid_space_x, int grid_space_y, int grid_space_z, unsigned char* mask )
{
    int x, y, z;
    double xe, ye, ze;
    double xBasis, yBasis, zBasis;
    int i0, j0, k0;
    int i, j, k;
    int index, id;
    int startBasis, startLocation;
    int maskID=0;

    int planeSize=xdim*ydim;
    index=0;
    for (z=0; z<zdim; z++)
    {
        for (y=0; y<ydim; y++)
        {
            for(x=0; x<xdim; x++)
            {
                index=z*planeSize+y*xdim+x;
                if (mask[index]>0)
                {
                    startBasis=maskID*64;
                    startLocation=maskID*4;

                    xe=static_cast<double>(x+u_x_existing[index]);
                    ye=static_cast<double>(y+u_y_existing[index]);
                    ze=static_cast<double>(z+u_z_existing[index]);

                    xBasis=xe/grid_space_x;
                    yBasis=ye/grid_space_y;
                    zBasis=ze/grid_space_z;

                    i0=static_cast<int>(xBasis)-1;
                    if (xBasis<0)
                    {
                        i0-=1;
                    }

                    j0=static_cast<int>(yBasis)-1;
                    if (yBasis<0)
                    {
                        j0-=1;
                    }

                    k0=static_cast<int>(zBasis)-1;
                    if (zBasis<0)
                    {
                        k0-=1;
                    }

                    id=0;
                    for (k=k0; k<k0+4; k++)
                    {
                        for (j=j0; j<j0+4; j++)
                        {
                            for (i=i0; i<i0+4; i++)
                            {
                                arrayBasis0[startBasis+id]=static_cast<float>(dcubic_bsp(xBasis-i)*dcubic_bsp(yBasis-j)*dcubic_bsp(zBasis-k));

                                arrayBasis1[startBasis+id]=static_cast<float>(dcubic_dbsp(xBasis-i)*dcubic_bsp(yBasis-j)*dcubic_bsp(zBasis-k)/grid_space_x);
                                arrayBasis2[startBasis+id]=static_cast<float>(dcubic_bsp(xBasis-i)*dcubic_dbsp(yBasis-j)*dcubic_bsp(zBasis-k)/grid_space_y);
                                arrayBasis3[startBasis+id]=static_cast<float>(dcubic_bsp(xBasis-i)*dcubic_bsp(yBasis-j)*dcubic_dbsp(zBasis-k)/grid_space_z);

                                arrayBasis4[startBasis+id]=static_cast<float>(dcubic_d2bsp(xBasis-i)*dcubic_bsp(yBasis-j)*dcubic_bsp(zBasis-k)/grid_space_x);
                                arrayBasis5[startBasis+id]=static_cast<float>(dcubic_bsp(xBasis-i)*dcubic_d2bsp(yBasis-j)*dcubic_bsp(zBasis-k)/grid_space_y);
                                arrayBasis6[startBasis+id]=static_cast<float>(dcubic_bsp(xBasis-i)*dcubic_bsp(yBasis-j)*dcubic_d2bsp(zBasis-k)/grid_space_z);

                                id++;
                            }
                        }
                    }

                    arrayLocation[startLocation]=i0;
                    arrayLocation[startLocation+1]=j0;
                    arrayLocation[startLocation+2]=k0;
                    arrayLocation[startLocation+3]=index;
                    //arrayLocation[start+4]=x;
                    //arrayLocation[start+5]=y;
                    //arrayLocation[start+6]=z;
                
                    maskID++;
                }

                //index++;
            }
        }
    }

}




void compose_u( float* u_x_total, float* u_y_total, float* u_z_total, float* u_x_existing, float* u_y_existing, float* u_z_existing, float* arrayBasis0, int* arrayLocation, float* c_x, float* c_y, float* c_z, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int maskCT2 )
{
    int coeffPlaneSize=num_element_x*num_element_y;
    int imgSize=xdim*ydim*zdim;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int xg,yg,zg;
    int idbx, idby, idbz;
    float basisProduct;

    for (index=0; index<imgSize; index++)
    {
        u_x_total[index]=u_x_existing[index];
        u_y_total[index]=u_y_existing[index];
        u_z_total[index]=u_z_existing[index];
    }
    
    for (maskID=0; maskID<maskCT2; maskID++)
    {
        //startBasis=maskID*64;
        startLocation=maskID*7;

        i0=arrayLocation[startLocation];
        j0=arrayLocation[startLocation+1];
        k0=arrayLocation[startLocation+2];
        index=arrayLocation[startLocation+3];
        xg=arrayLocation[startLocation+4];
        yg=arrayLocation[startLocation+5];
        zg=arrayLocation[startLocation+6];

        idbz=zg*4;
        for (k=k0; k<k0+4; k++)
        {
            idby=yg*4;
            for (j=j0; j<j0+4; j++)
            {
                idbx=xg*4;
                for (i=i0; i<i0+4; i++)
                {
                    if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                    {
                        basisProduct=arrayBasis0[idbx]*arrayBasis0[idby]*arrayBasis0[idbz];
                        /*if (index==618)
                        {
                            std::cout<<"("<<idbx+i-i0<<", "<<idby+j-j0<<", "<<idbz+k-k0<<"): "<<basisProduct<<" ";
                        }*/
                        //coeff global index (storage form)
                        coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);          
                        //bspIndex=idg*64+id;
                        u_x_total[index]+=c_x[coeffIndex]*basisProduct;
                        u_y_total[index]+=c_y[coeffIndex]*basisProduct;
                        u_z_total[index]+=c_z[coeffIndex]*basisProduct;
                    }
                    idbx++;
                }
                idby++;
            }
            idbz++;
        }

    //if (index==618)
    //std::cout<<index<<" total_disp_x: "<<u_x_total[index]<<std::endl;
    }
    

}

void compose_u10( float* u_x_total, float* u_y_total, float* u_z_total, float* u_x_existing, float* u_y_existing, float* u_z_existing, float* arrayBasis0, int* arrayLocation, float* c_x, float* c_y, float* c_z, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int maskCT2 )
{
    int coeffPlaneSize=num_element_x*num_element_y;
    int imgSize=xdim*ydim*zdim;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int xg,yg,zg;
    int idbx, idby, idbz;
    float basisProduct;

    for (index=0; index<imgSize; index++)
    {
        u_x_total[index]=u_x_existing[index];
        u_y_total[index]=u_y_existing[index];
        u_z_total[index]=u_z_existing[index];
    }
    
    for (maskID=0; maskID<maskCT2; maskID++)
    {
        //startBasis=maskID*64;
        startLocation=maskID*10;

        i0=arrayLocation[startLocation];
        j0=arrayLocation[startLocation+1];
        k0=arrayLocation[startLocation+2];
        index=arrayLocation[startLocation+3];
        xg=arrayLocation[startLocation+4];
        yg=arrayLocation[startLocation+5];
        zg=arrayLocation[startLocation+6];

        idbz=zg*4;
        for (k=k0; k<k0+4; k++)
        {
            idby=yg*4;
            for (j=j0; j<j0+4; j++)
            {
                idbx=xg*4;
                for (i=i0; i<i0+4; i++)
                {
                    if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                    {
                        basisProduct=arrayBasis0[idbx]*arrayBasis0[idby]*arrayBasis0[idbz];
                        /*if (index==618)
                        {
                            std::cout<<"("<<idbx+i-i0<<", "<<idby+j-j0<<", "<<idbz+k-k0<<"): "<<basisProduct<<" ";
                        }*/
                        //coeff global index (storage form)
                        coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);          
                        //bspIndex=idg*64+id;
                        u_x_total[index]+=c_x[coeffIndex]*basisProduct;
                        u_y_total[index]+=c_y[coeffIndex]*basisProduct;
                        u_z_total[index]+=c_z[coeffIndex]*basisProduct;
                    }
                    idbx++;
                }
                idby++;
            }
            idbz++;
        }

    //if (index==618)
    //std::cout<<index<<" total_disp_x: "<<u_x_total[index]<<std::endl;
    }
    

}


void compose_u9( float* u_x_total, float* u_y_total, float* u_z_total, float* u_x_existing, float* u_y_existing, float* u_z_existing, float* arrayBasis0, short int* arrayLocation, float* c_x, float* c_y, float* c_z, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int maskCT2 )
{
    int coeffPlaneSize=num_element_x*num_element_y;
    int planeSize=xdim*ydim;
    int imgSize=xdim*ydim*zdim;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int xg,yg,zg;
    int idbx, idby, idbz;
    float basisProduct;
    int x, y, z;

    for (index=0; index<imgSize; index++)
    {
        u_x_total[index]=u_x_existing[index];
        u_y_total[index]=u_y_existing[index];
        u_z_total[index]=u_z_existing[index];
    }
    
    for (maskID=0; maskID<maskCT2; maskID++)
    {
        //startBasis=maskID*64;
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        i0=arrayLocation[startLocation+3];
        j0=arrayLocation[startLocation+4];
        k0=arrayLocation[startLocation+5];

        xg=arrayLocation[startLocation+6];
        yg=arrayLocation[startLocation+7];
        zg=arrayLocation[startLocation+8];

        idbz=zg*4;
        for (k=k0; k<k0+4; k++)
        {
            idby=yg*4;
            for (j=j0; j<j0+4; j++)
            {
                idbx=xg*4;
                for (i=i0; i<i0+4; i++)
                {
                    if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                    {
                        basisProduct=arrayBasis0[idbx]*arrayBasis0[idby]*arrayBasis0[idbz];
                        /*if (index==618)
                        {
                            std::cout<<"("<<idbx+i-i0<<", "<<idby+j-j0<<", "<<idbz+k-k0<<"): "<<basisProduct<<" ";
                        }*/
                        //coeff global index (storage form)
                        coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);          
                        //bspIndex=idg*64+id;
                        u_x_total[index]+=c_x[coeffIndex]*basisProduct;
                        u_y_total[index]+=c_y[coeffIndex]*basisProduct;
                        u_z_total[index]+=c_z[coeffIndex]*basisProduct;
                    }
                    idbx++;
                }
                idby++;
            }
            idbz++;
        }

    //if (index==618)
    //std::cout<<index<<" total_disp_x: "<<u_x_total[index]<<std::endl;
    }
    

}






//compose disp
void compose_disp(float* u12_x, float* u12_y, float* u12_z, float* u23_x, float* u23_y, float* u23_z, float* u13_x, float* u13_y, float* u13_z, int xdim, int ydim, int zdim)
{
    int planeSize=xdim*ydim;
    float u111,u112,u121,u122,u211,u212,u221,u222;
    float u12_x_def, u12_y_def, u12_z_def;
    int index=0;
    for(int z=0;z<zdim;z++)
    {
        for(int y=0;y<ydim;y++)
        {
            for(int x=0;x<xdim;x++)
            {
                float hx=x+u23_x[index];
                float hy=y+u23_y[index];
                float hz=z+u23_z[index];
                
                //x1,y1,z1---starting point; 
                //x2,y2,z2---dist from starting point
				int x1=static_cast<int>(hx);
                if (hx<0)
                    x1-=1;
				float x2=hx-x1;
				int y1=static_cast<int>(hy);
                if (hy<0)
                    y1-=1;
				float y2=hy-y1;
				int z1=static_cast<int>(hz);
                if (hz<0)
                    z1-=1;
				float z2=hz-z1;

                //circular indexing
                if (x1<0)
                    x1=0;
                if (x1>=xdim)
                    x1=xdim-1;
                if (y1<0)
                    y1=0;
                if (y1>=ydim)
                    y1=ydim-1;
                if (z1<0)
                    z1=0;
                if (z1>=zdim)
                    z1=zdim-1;
                
                //x3,y3,z3---last point, and circular indexing
                int x3=x1+1;
                int y3=y1+1;
                int z3=z1+1;
                if (x3<0)
                    x3=0;
                if (x3>=xdim)
                    x3=xdim-1;
                if (y3<0)
                    y3=0;
                if (y3>=ydim)
                    y3=ydim-1;
                if (z3<0)
                    z3=0;
                if (z3>=zdim)
                    z3=zdim-1;

                int id1=z1*planeSize+y1*xdim+x1;
                int id2=z1*planeSize+y1*xdim+x3;
                int id3=z1*planeSize+y3*xdim+x1;
                int id4=z1*planeSize+y3*xdim+x3;
                int id5=z3*planeSize+y1*xdim+x1;
                int id6=z3*planeSize+y1*xdim+x3;
                int id7=z3*planeSize+y3*xdim+x1;
                int id8=z3*planeSize+y3*xdim+x3;
                
                u111=u12_x[id1];
				u211=u12_x[id2];
				u121=u12_x[id3];
				u221=u12_x[id4];
				u112=u12_x[id5];
				u212=u12_x[id6];
				u122=u12_x[id7];
				u222=u12_x[id8];
				u12_x_def=(1-x2)*(1-y2)*(1-z2)*u111+x2*(1-y2)*(1-z2)*u211+(1-x2)*y2*(1-z2)*u121+x2*y2*(1-z2)*u221+(1-x2)*(1-y2)*z2*u112+x2*(1-y2)*z2*u212+(1-x2)*y2*z2*u122+x2*y2*z2*u222;
                u13_x[index]=u23_x[index]+u12_x_def;

                u111=u12_y[id1];
				u211=u12_y[id2];
				u121=u12_y[id3];
				u221=u12_y[id4];
				u112=u12_y[id5];
				u212=u12_y[id6];
				u122=u12_y[id7];
				u222=u12_y[id8];
				u12_y_def=(1-x2)*(1-y2)*(1-z2)*u111+x2*(1-y2)*(1-z2)*u211+(1-x2)*y2*(1-z2)*u121+x2*y2*(1-z2)*u221+(1-x2)*(1-y2)*z2*u112+x2*(1-y2)*z2*u212+(1-x2)*y2*z2*u122+x2*y2*z2*u222;
                u13_y[index]=u23_y[index]+u12_y_def;

                u111=u12_z[id1];
				u211=u12_z[id2];
				u121=u12_z[id3];
				u221=u12_z[id4];
				u112=u12_z[id5];
				u212=u12_z[id6];
				u122=u12_z[id7];
				u222=u12_z[id8];
				u12_z_def=(1-x2)*(1-y2)*(1-z2)*u111+x2*(1-y2)*(1-z2)*u211+(1-x2)*y2*(1-z2)*u121+x2*y2*(1-z2)*u221+(1-x2)*(1-y2)*z2*u112+x2*(1-y2)*z2*u212+(1-x2)*y2*z2*u122+x2*y2*z2*u222;
                u13_z[index]=u23_z[index]+u12_z_def;
                
                index++;
            }
        }
    }

}

//deform image
void deformed1_inner(float* img, float* tmp, int xdim, int ydim, int zdim, float* u_x_for, float* u_y_for, float* u_z_for, float zero)
{
    int imgSize=xdim*ydim*zdim;
    int planeSize=xdim*ydim;
	float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;

    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                tmp[index]=zero;

				hx=x+u_x_for[index];
				hy=y+u_y_for[index];
				hz=z+u_z_for[index];

                //x1,y1,z1---starting point; 
                //x2,y2,z2---dist from starting point
				x1=static_cast<int>(hx);
                if (hx<0)
                    x1-=1;
				x2=hx-x1;

				y1=static_cast<int>(hy);
                if (hy<0)
                    y1-=1;
				y2=hy-y1;

				z1=static_cast<int>(hz);
                if (hz<0)
                    z1-=1;
				z2=hz-z1;

                //x3,y3,z3---last point
                x3=x1+1;
                y3=y1+1;
                z3=z1+1;

                //if (id0>=0 && id7<imgSize)
                if (x1>=0 && x3<xdim && y1>=0 && y3<ydim && z1>=0 && z3<zdim)
                {
                    id0=z1*planeSize+y1*xdim+x1;
                    id1=id0+1;
                    id2=id0+xdim;
                    id3=id2+1;
                    id4=id0+planeSize;
                    id5=id4+1;
                    id6=id4+xdim;
                    id7=id6+1;

                    dist0=(1-x2)*(1-y2)*(1-z2);
                    dist1=x2*(1-y2)*(1-z2);
                    dist2=(1-x2)*y2*(1-z2);
                    dist3=x2*y2*(1-z2);
                    dist4=(1-x2)*(1-y2)*z2;
                    dist5=x2*(1-y2)*z2;
                    dist6=(1-x2)*y2*z2;
                    dist7=x2*y2*z2;

			    	tmp[index]=dist0*img[id0]+dist1*img[id1]+dist2*img[id2]+dist3*img[id3]+dist4*img[id4]+dist5*img[id5]+dist6*img[id6]+dist7*img[id7];

                }

                index++;
			}

    for (int i=0; i<imgSize; i++)
    {
        img[i]=tmp[i];
    }
}

//deform image
void deformed2_inner(float* vdeformed, float* imagein, int xdim, int ydim, int zdim, float* u_x_for, float* u_y_for, float* u_z_for, float zero)
{
    int planeSize=xdim*ydim;
	float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;
    //std::cout<<"xdim="<<xdim<<", ydim="<<ydim<<", zdim="<<zdim<<std::endl;
    //std::cout<<"imgSize="<<imgSize<<std::endl;

    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                vdeformed[index]=zero;

				hx=x+u_x_for[index];
				hy=y+u_y_for[index];
				hz=z+u_z_for[index];

                //x1,y1,z1---starting point; 
                //x2,y2,z2---dist from starting point
				x1=static_cast<int>(hx);
                if (hx<0)
                    x1-=1;
				x2=hx-x1;

				y1=static_cast<int>(hy);
                if (hy<0)
                    y1-=1;
				y2=hy-y1;

				z1=static_cast<int>(hz);
                if (hz<0)
                    z1-=1;
				z2=hz-z1;

                //x3,y3,z3---last point
                x3=x1+1;
                y3=y1+1;
                z3=z1+1;

                /*if (x==53 && y==220 && z==116)
                {
                    std::cout<<"look at me"<<std::endl;
                    std::cout<<"("<<x<<", "<<y<<", "<<z<<") ---> "<<"("<<hx<<", "<<hy<<", "<<hz<<")"<<std::endl;
                    std::cout<<"id0 = "<<id0<<", id7 = "<<id7<<std::endl;
                }*/

                //if (id0>=0 && id7<imgSize)
                if (x1>=0 && x3<xdim && y1>=0 && y3<ydim && z1>=0 && z3<zdim)
                {
                    id0=z1*planeSize+y1*xdim+x1;
                    id1=id0+1;
                    id2=id0+xdim;
                    id3=id2+1;
                    id4=id0+planeSize;
                    id5=id4+1;
                    id6=id4+xdim;
                    id7=id6+1;

                    dist0=(1-x2)*(1-y2)*(1-z2);
                    dist1=x2*(1-y2)*(1-z2);
                    dist2=(1-x2)*y2*(1-z2);
                    dist3=x2*y2*(1-z2);
                    dist4=(1-x2)*(1-y2)*z2;
                    dist5=x2*(1-y2)*z2;
                    dist6=(1-x2)*y2*z2;
                    dist7=x2*y2*z2;

			    	vdeformed[index]=dist0*imagein[id0]+dist1*imagein[id1]+dist2*imagein[id2]+dist3*imagein[id3]+dist4*imagein[id4]+dist5*imagein[id5]+dist6*imagein[id6]+dist7*imagein[id7];

                    /*if (x==108 && y==161 && z==151)
                    {
                        std::cout<<x2<<", "<<y2<<", "<<z2<<";"<<std::endl;
                        std::cout<<dist0<<", "<<dist1<<", "<<dist2<<", "<<dist3<<";"<<std::endl;
                        std::cout<<dist4<<", "<<dist5<<", "<<dist6<<", "<<dist7<<";"<<std::endl;
                        std::cout<<imagein[id0]<<", "<<imagein[id1]<<", "<<imagein[id2]<<", "<<imagein[id3]<<";"<<std::endl;
                        std::cout<<imagein[id4]<<", "<<imagein[id5]<<", "<<imagein[id6]<<", "<<imagein[id7]<<"."<<std::endl;
                        std::cout<<"def value: "<<vdeformed[index]<<std::endl;
                        std::cout<<id0<<", "<<id1<<", "<<id2<<", "<<id3<<";"<<std::endl;
                        std::cout<<id4<<", "<<id5<<", "<<id6<<", "<<id7<<";"<<std::endl;
                    }*/

                }

                index++;
			}

}

//deform image
void deformed1_all(float* img, float* tmp, int xdim, int ydim, int zdim, float* u_x_for, float* u_y_for, float* u_z_for, float zero)
{
    int imgSize=xdim*ydim*zdim;
    int planeSize=xdim*ydim;
	float image111,image112,image121,image122,image211,image212,image221,image222;
    int id;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;
    //std::cout<<"xdim="<<xdim<<", ydim="<<ydim<<", zdim="<<zdim<<std::endl;

    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
				hx=x+u_x_for[index];
				hy=y+u_y_for[index];
				hz=z+u_z_for[index];

                //x1,y1,z1---starting point; 
                //x2,y2,z2---dist from starting point
				x1=static_cast<int>(hx);
                if (hx<0)
                    x1-=1;
				x2=hx-x1;

				y1=static_cast<int>(hy);
                if (hy<0)
                    y1-=1;
				y2=hy-y1;

				z1=static_cast<int>(hz);
                if (hz<0)
                    z1-=1;
				z2=hz-z1;

                //x3,y3,z3---last point
                x3=x1+1;
                y3=y1+1;
                z3=z1+1;

                image111=zero;
				image211=zero;
				image121=zero;
				image221=zero;
				image112=zero;
				image212=zero;
				image122=zero;
				image222=zero;

                //linear interpolating
                id=z1*planeSize+y1*xdim+x1;
                if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
				    image111=img[id];//(x1,y1,z1);
                if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
				    image211=img[id+1];//(x3,y1,z1);
                if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
				    image121=img[id+xdim];//(x1,y3,z1);
                if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
				    image221=img[id+xdim+1];//(x3,y3,z1);
                if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
				    image112=img[id+planeSize];//(x1,y1,z3);
                if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
				    image212=img[id+planeSize+1];//(x3,y1,z3);
                if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
				    image122=img[id+planeSize+xdim];//(x1,y3,z3);
                if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
				    image222=img[id+planeSize+xdim+1];//(x3,y3,z3);

				tmp[index]=(1-x2)*(1-y2)*(1-z2)*image111+x2*(1-y2)*(1-z2)*image211+(1-x2)*y2*(1-z2)*image121+x2*y2*(1-z2)*image221+(1-x2)*(1-y2)*z2*image112+x2*(1-y2)*z2*image212+(1-x2)*y2*z2*image122+x2*y2*z2*image222;

                index++;
			}

    for (int i=0; i<imgSize; i++)
    {
        img[i]=tmp[i];
    }
}


//deform image
void deformed2_all(float* vdeformed, float* imagein, int xdim, int ydim, int zdim, float* u_x_for, float* u_y_for, float* u_z_for, float zero)
{
    int planeSize=xdim*ydim;
	float image111,image112,image121,image122,image211,image212,image221,image222;
    int id;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;
    //std::cout<<"xdim="<<xdim<<", ydim="<<ydim<<", zdim="<<zdim<<std::endl;

    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
				hx=x+u_x_for[index];
				hy=y+u_y_for[index];
				hz=z+u_z_for[index];

                //x1,y1,z1---starting point; 
                //x2,y2,z2---dist from starting point
				x1=static_cast<int>(hx);
                if (hx<0)
                    x1-=1;
				x2=hx-x1;

				y1=static_cast<int>(hy);
                if (hy<0)
                    y1-=1;
				y2=hy-y1;

				z1=static_cast<int>(hz);
                if (hz<0)
                    z1-=1;
				z2=hz-z1;

                //x3,y3,z3---last point
                x3=x1+1;
                y3=y1+1;
                z3=z1+1;

                image111=zero;
				image211=zero;
				image121=zero;
				image221=zero;
				image112=zero;
				image212=zero;
				image122=zero;
				image222=zero;

                //linear interpolating
                id=z1*planeSize+y1*xdim+x1;
                if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
				    image111=imagein[id];//(x1,y1,z1);
                if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
				    image211=imagein[id+1];//(x3,y1,z1);
                if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
				    image121=imagein[id+xdim];//(x1,y3,z1);
                if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
				    image221=imagein[id+xdim+1];//(x3,y3,z1);
                if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
				    image112=imagein[id+planeSize];//(x1,y1,z3);
                if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
				    image212=imagein[id+planeSize+1];//(x3,y1,z3);
                if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
				    image122=imagein[id+planeSize+xdim];//(x1,y3,z3);
                if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
				    image222=imagein[id+planeSize+xdim+1];//(x3,y3,z3);

				vdeformed[index]=(1-x2)*(1-y2)*(1-z2)*image111+x2*(1-y2)*(1-z2)*image211+(1-x2)*y2*(1-z2)*image121+x2*y2*(1-z2)*image221+(1-x2)*(1-y2)*z2*image112+x2*(1-y2)*z2*image212+(1-x2)*y2*z2*image122+x2*y2*z2*image222;

                /*if (x==108 && y==161 && z==151)
                {
                    std::cout<<"look at me"<<std::endl;
                    std::cout<<"("<<x<<", "<<y<<", "<<z<<") ---> "<<"("<<hx<<", "<<hy<<", "<<hz<<")"<<std::endl;
                    std::cout<<image111<<", "<<image211<<", "<<image121<<", "<<image221<<";"<<std::endl;
                    std::cout<<image112<<", "<<image212<<", "<<image122<<", "<<image222<<"."<<std::endl;
                    std::cout<<"def value: "<<vdeformed[index]<<std::endl;
                }*/

                index++;
			}
}


//deform image
void mask_deformed2_all(unsigned char* vdeformed, unsigned char* imagein, int xdim, int ydim, int zdim, float* u_x_for, float* u_y_for, float* u_z_for, unsigned char zero)
{
    int planeSize=xdim*ydim;
    int id;
    float hx, hy, hz;
    int x1, y1, z1;
    //std::cout<<"xdim="<<xdim<<", ydim="<<ydim<<", zdim="<<zdim<<std::endl;

    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
				hx=x+u_x_for[index];
				hy=y+u_y_for[index];
				hz=z+u_z_for[index];

                //x1,y1,z1---nearest point; 
				x1=round(hx);
				y1=round(hy);
				z1=round(hz);

                //nearest interpolating
                id=z1*planeSize+y1*xdim+x1;
                if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
				    vdeformed[index]=imagein[id];
                else
                    vdeformed[index]=zero;
                /*if (x==108 && y==161 && z==151)
                {
                    std::cout<<"look at me"<<std::endl;
                    std::cout<<"("<<x<<", "<<y<<", "<<z<<") ---> "<<"("<<hx<<", "<<hy<<", "<<hz<<")"<<std::endl;
                    std::cout<<image111<<", "<<image211<<", "<<image121<<", "<<image221<<";"<<std::endl;
                    std::cout<<image112<<", "<<image212<<", "<<image122<<", "<<image222<<"."<<std::endl;
                    std::cout<<"def value: "<<vdeformed[index]<<std::endl;
                }*/

                index++;
			}
}



void deformedM(float* vdeformed, float* imagein, int xdim, int ydim, int zdim, float* u_x_for, float* u_y_for, float* u_z_for, float zero, unsigned char* mask)
{
    int planeSize=xdim*ydim;
	float image111,image112,image121,image122,image211,image212,image221,image222;
    int id;
    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (mask[index]>0)
                {
			    	float hx=x+u_x_for[index];
			    	float hy=y+u_y_for[index];
			    	float hz=z+u_z_for[index];

                    //x1,y1,z1---starting point; 
                    //x2,y2,z2---dist from starting point
			    	int x1=static_cast<int>(hx);
                    if (hx<0)
                        x1-=1;
			    	float x2=hx-x1;
			    	int y1=static_cast<int>(hy);
                    if (hy<0)
                        y1-=1;
			    	float y2=hy-y1;
			    	int z1=static_cast<int>(hz);
                    if (hz<0)
                        z1-=1;
			    	float z2=hz-z1;

                    //x3,y3,z3---last point
                    int x3=x1+1;
                    int y3=y1+1;
                    int z3=z1+1;

                    image111=zero;
			    	image211=zero;
			    	image121=zero;
			    	image221=zero;
			    	image112=zero;
			    	image212=zero;
			    	image122=zero;
			    	image222=zero;
                    //linear interpolating
                    /*if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
			    	    image111=imagein[z1*planeSize+y1*xdim+x1];//(x1,y1,z1);
                    if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
			    	    image211=imagein[z1*planeSize+y1*xdim+x3];//(x3,y1,z1);
                    if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
			    	    image121=imagein[z1*planeSize+y3*xdim+x1];//(x1,y3,z1);
                    if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
			    	    image221=imagein[z1*planeSize+y3*xdim+x3];//(x3,y3,z1);
                    if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
			    	    image112=imagein[z3*planeSize+y1*xdim+x1];//(x1,y1,z3);
                    if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
			    	    image212=imagein[z3*planeSize+y1*xdim+x3];//(x3,y1,z3);
                    if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
			    	    image122=imagein[z3*planeSize+y3*xdim+x1];//(x1,y3,z3);
                    if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
			    	    image222=imagein[z3*planeSize+y3*xdim+x3];//(x3,y3,z3);*/

                    id=z1*planeSize+y1*xdim+x1;
                    if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
				        image111=imagein[id];//(x1,y1,z1);
                    if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
				        image211=imagein[id+1];//(x3,y1,z1);
                    if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
				        image121=imagein[id+xdim];//(x1,y3,z1);
                    if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
				        image221=imagein[id+xdim+1];//(x3,y3,z1);
                    if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
				        image112=imagein[id+planeSize];//(x1,y1,z3);
                    if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
				        image212=imagein[id+planeSize+1];//(x3,y1,z3);
                    if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
				        image122=imagein[id+planeSize+xdim];//(x1,y3,z3);
                    if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
				        image222=imagein[id+planeSize+xdim+1];//(x3,y3,z3);

			    	vdeformed[index]=(1-x2)*(1-y2)*(1-z2)*image111+x2*(1-y2)*(1-z2)*image211+(1-x2)*y2*(1-z2)*image121+x2*y2*(1-z2)*image221+(1-x2)*(1-y2)*z2*image112+x2*(1-y2)*z2*image212+(1-x2)*y2*z2*image122+x2*y2*z2*image222;
                }
                else
                {
                    vdeformed[index]=zero;
                }
                index++;
			}
}


void index2coord(int index, int& x, int& y, int& z, int xdim, int ydim , int zdim)
{
    z=int(index/(xdim*ydim));
    int remainder=index%(xdim*ydim);
    if (remainder==0)
    {
        z--;
        y=ydim-1;
        x=xdim-1;
    }
    else
    {
        y=int(remainder/xdim);
        x=remainder%xdim;
        if (x==0)
        {
            y--;
            x=xdim-1;
        }
    }

}


//deform image
float deform_pnt(float* imagein, float* u_x_for, float* u_y_for, float* u_z_for, int index, int xdim, int ydim, int zdim, float zero)
{
    float defVal=0.0;
	float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    int imgSize=xdim*ydim*zdim;
    int planeSize=xdim*ydim;
	
    int x,y,z;
    index2coord(index, x, y, z, xdim, ydim, zdim);

				float hx=x+u_x_for[index];
				float hy=y+u_y_for[index];
				float hz=z+u_z_for[index];

                //x1,y1,z1---starting point; 
                //x2,y2,z2---dist from starting point
				int x1=static_cast<int>(hx);
                if (hx<0)
                    x1-=1;
				float x2=hx-x1;
				int y1=static_cast<int>(hy);
                if (hy<0)
                    y1-=1;
				float y2=hy-y1;
				int z1=static_cast<int>(hz);
                if (hz<0)
                    z1-=1;
				float z2=hz-z1;

                id0=z1*planeSize+y1*xdim+x1;
                id7=id0+planeSize+xdim+1;
                if (id0>=0 && id7<imgSize)
                {
                    id1=id0+1;
                    id2=id0+xdim;
                    id3=id2+1;
                    id4=id0+planeSize;
                    id5=id4+1;
                    id6=id4+xdim;

                    dist0=(1-x2)*(1-y2)*(1-z2);
                    dist1=x2*(1-y2)*(1-z2);
                    dist2=(1-x2)*y2*(1-z2);
                    dist3=x2*y2*(1-z2);
                    dist4=(1-x2)*(1-y2)*z2;
                    dist5=x2*(1-y2)*z2;
                    dist6=(1-x2)*y2*z2;
                    dist7=x2*y2*z2;

			    	defVal=dist0*imagein[id0]+dist1*imagein[id1]+dist2*imagein[id2]+dist3*imagein[id3]+dist4*imagein[id4]+dist5*imagein[id5]+dist6*imagein[id6]+dist7*imagein[id7];

                }


    return defVal;

}




void compute_u( float* u_x_total, float* u_y_total, float* u_z_total, float* arrayBasis, int* arrayLocation, float* c_x, float* c_y, float* c_z, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int maskCT2 )
{
    int coeffPlaneSize=num_element_x*num_element_y;
    int imgSize=xdim*ydim*zdim;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int idbx,idby,idbz;
    float basisProduct;
    int xg, yg, zg;

    for (index=0; index<imgSize; index++)
    {
        u_x_total[index]=0;
        u_y_total[index]=0;
        u_z_total[index]=0;
    }
    
    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*7;

        i0=arrayLocation[startLocation];
        j0=arrayLocation[startLocation+1];
        k0=arrayLocation[startLocation+2];
        index=arrayLocation[startLocation+3];
        xg=arrayLocation[startLocation+4];
        yg=arrayLocation[startLocation+5];
        zg=arrayLocation[startLocation+6];

        idbz=zg*4;
        for (k=k0; k<k0+4; k++)
        {
            idby=yg*4;
            for (j=j0; j<j0+4; j++)
            {
                idbx=xg*4;
                for (i=i0; i<i0+4; i++)
                {
                    if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                    {
                        basisProduct=arrayBasis[idbx]*arrayBasis[idby]*arrayBasis[idbz];
                        //coeff global index (storage form)
                        coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);              
                        u_x_total[index]+=c_x[coeffIndex]*basisProduct;
                        u_y_total[index]+=c_y[coeffIndex]*basisProduct;
                        u_z_total[index]+=c_z[coeffIndex]*basisProduct;
                    }
                    idbx++;
                }
                idby++;
            }
            idbz++;
        }

    }

}

void compute_u10( float* u_x_total, float* u_y_total, float* u_z_total, float* arrayBasis, int* arrayLocation, float* c_x, float* c_y, float* c_z, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int maskCT2 )
{
    int coeffPlaneSize=num_element_x*num_element_y;
    int imgSize=xdim*ydim*zdim;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int idbx,idby,idbz;
    float basisProduct;
    int xg, yg, zg;

    for (index=0; index<imgSize; index++)
    {
        u_x_total[index]=0;
        u_y_total[index]=0;
        u_z_total[index]=0;
    }
    
    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*10;

        i0=arrayLocation[startLocation];
        j0=arrayLocation[startLocation+1];
        k0=arrayLocation[startLocation+2];
        index=arrayLocation[startLocation+3];
        xg=arrayLocation[startLocation+4];
        yg=arrayLocation[startLocation+5];
        zg=arrayLocation[startLocation+6];

        idbz=zg*4;
        for (k=k0; k<k0+4; k++)
        {
            idby=yg*4;
            for (j=j0; j<j0+4; j++)
            {
                idbx=xg*4;
                for (i=i0; i<i0+4; i++)
                {
                    if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                    {
                        basisProduct=arrayBasis[idbx]*arrayBasis[idby]*arrayBasis[idbz];
                        //coeff global index (storage form)
                        coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);              
                        u_x_total[index]+=c_x[coeffIndex]*basisProduct;
                        u_y_total[index]+=c_y[coeffIndex]*basisProduct;
                        u_z_total[index]+=c_z[coeffIndex]*basisProduct;
                    }
                    idbx++;
                }
                idby++;
            }
            idbz++;
        }

    }

}

void compute_u9( float* u_x_total, float* u_y_total, float* u_z_total, float* arrayBasis, short int* arrayLocation, float* c_x, float* c_y, float* c_z, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int maskCT2 )
{
    int coeffPlaneSize=num_element_x*num_element_y;
    int planeSize=xdim*ydim;
    int imgSize=xdim*ydim*zdim;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int idbx,idby,idbz;
    float basisProduct;
    int xg, yg, zg;
    int x, y, z;

    for (index=0; index<imgSize; index++)
    {
        u_x_total[index]=0;
        u_y_total[index]=0;
        u_z_total[index]=0;
    }
    
    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        i0=arrayLocation[startLocation+3];
        j0=arrayLocation[startLocation+4];
        k0=arrayLocation[startLocation+5];

        xg=arrayLocation[startLocation+6];
        yg=arrayLocation[startLocation+7];
        zg=arrayLocation[startLocation+8];
        


        idbz=zg*4;
        for (k=k0; k<k0+4; k++)
        {
            idby=yg*4;
            for (j=j0; j<j0+4; j++)
            {
                idbx=xg*4;
                for (i=i0; i<i0+4; i++)
                {
                    if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                    {
                        basisProduct=arrayBasis[idbx]*arrayBasis[idby]*arrayBasis[idbz];
                        //coeff global index (storage form)
                        coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);              
                        u_x_total[index]+=c_x[coeffIndex]*basisProduct;
                        u_y_total[index]+=c_y[coeffIndex]*basisProduct;
                        u_z_total[index]+=c_z[coeffIndex]*basisProduct;
                    }
                    idbx++;
                }
                idby++;
            }
            idbz++;
        }

    }

}



//deform image
void deformed0_derivatives(float* uxx_def, float* uxy_def, float* uxz_def, float* uyx_def, float* uyy_def, float* uyz_def, float* uzx_def, float* uzy_def, float* uzz_def, float* uxx, float* uxy, float* uxz, float* uyx, float* uyy, float* uyz, float* uzx, float* uzy, float* uzz, int xdim, int ydim, int zdim, float* u_x_for, float* u_y_for, float* u_z_for, float zero)
{
    int imgSize=xdim*ydim*zdim;
    int planeSize=xdim*ydim;
	float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
				float hx=x+u_x_for[index];
				float hy=y+u_y_for[index];
				float hz=z+u_z_for[index];

                //x1,y1,z1---starting point; 
                //x2,y2,z2---dist from starting point
				int x1=static_cast<int>(hx);
                if (hx<0)
                    x1-=1;
				float x2=hx-x1;
				int y1=static_cast<int>(hy);
                if (hy<0)
                    y1-=1;
				float y2=hy-y1;
				int z1=static_cast<int>(hz);
                if (hz<0)
                    z1-=1;
				float z2=hz-z1;

                //x3,y3,z3---last point
                /*int x3=x1+1;
                int y3=y1+1;
                int z3=z1+1;*/

                id0=z1*planeSize+y1*xdim+x1;
                id7=id0+planeSize+xdim+1;
                if (id0>=0 && id7<imgSize)
                {
                    id1=id0+1;
                    id2=id0+xdim;
                    id3=id2+1;
                    id4=id0+planeSize;
                    id5=id4+1;
                    id6=id4+xdim;

                    /*if (id0>=0 && id0<imgSize)
			    	    image111=uxx[id0];//(x1,y1,z1);
                    if (id1>=0 && id1<imgSize)
			    	    image211=uxx[id1];//(x3,y1,z1);
                    if (id2>=0 && id2<imgSize)
			    	    image121=uxx[id2];//(x1,y3,z1);
                    if (id3>=0 && id3<imgSize)
			    	    image221=uxx[id3];//(x3,y3,z1);
                    if (id4>=0 && id4<imgSize)
			    	    image112=uxx[id4];//(x1,y1,z3);
                    if (id5>=0 && id5<imgSize)
			    	    image212=uxx[id5];//(x3,y1,z3);
                    if (id6>=0 && id6<imgSize)
			    	    image122=uxx[id6];//(x1,y3,z3);
                    if (id7>=0 && id7<imgSize)
			    	    image222=uxx[id7];//(x3,y3,z3);*/

                    dist0=(1-x2)*(1-y2)*(1-z2);
                    dist1=x2*(1-y2)*(1-z2);
                    dist2=(1-x2)*y2*(1-z2);
                    dist3=x2*y2*(1-z2);
                    dist4=(1-x2)*(1-y2)*z2;
                    dist5=x2*(1-y2)*z2;
                    dist6=(1-x2)*y2*z2;
                    dist7=x2*y2*z2;

			    	uxx_def[index]=dist0*uxx[id0]+dist1*uxx[id1]+dist2*uxx[id2]+dist3*uxx[id3]+dist4*uxx[id4]+dist5*uxx[id5]+dist6*uxx[id6]+dist7*uxx[id7];
			    	uxy_def[index]=dist0*uxy[id0]+dist1*uxy[id1]+dist2*uxy[id2]+dist3*uxy[id3]+dist4*uxy[id4]+dist5*uxy[id5]+dist6*uxy[id6]+dist7*uxy[id7];
			    	uxz_def[index]=dist0*uxz[id0]+dist1*uxz[id1]+dist2*uxz[id2]+dist3*uxz[id3]+dist4*uxz[id4]+dist5*uxz[id5]+dist6*uxz[id6]+dist7*uxz[id7];

			    	uyx_def[index]=dist0*uyx[id0]+dist1*uyx[id1]+dist2*uyx[id2]+dist3*uyx[id3]+dist4*uyx[id4]+dist5*uyx[id5]+dist6*uyx[id6]+dist7*uyx[id7];
			    	uyy_def[index]=dist0*uyy[id0]+dist1*uyy[id1]+dist2*uyy[id2]+dist3*uyy[id3]+dist4*uyy[id4]+dist5*uyy[id5]+dist6*uyy[id6]+dist7*uyy[id7];
			    	uyz_def[index]=dist0*uyz[id0]+dist1*uyz[id1]+dist2*uyz[id2]+dist3*uyz[id3]+dist4*uyz[id4]+dist5*uyz[id5]+dist6*uyz[id6]+dist7*uyz[id7];

			    	uzx_def[index]=dist0*uzx[id0]+dist1*uzx[id1]+dist2*uzx[id2]+dist3*uzx[id3]+dist4*uzx[id4]+dist5*uzx[id5]+dist6*uzx[id6]+dist7*uzx[id7];
			    	uzy_def[index]=dist0*uzy[id0]+dist1*uzy[id1]+dist2*uzy[id2]+dist3*uzy[id3]+dist4*uzy[id4]+dist5*uzy[id5]+dist6*uzy[id6]+dist7*uzy[id7];
			    	uzz_def[index]=dist0*uzz[id0]+dist1*uzz[id1]+dist2*uzz[id2]+dist3*uzz[id3]+dist4*uzz[id4]+dist5*uzz[id5]+dist6*uzz[id6]+dist7*uzz[id7];

                }
                else
                {
                    uxx_def[index]=zero;
                    uxy_def[index]=zero;
                    uxz_def[index]=zero;

                    uyx_def[index]=zero;
                    uyy_def[index]=zero;
                    uyz_def[index]=zero;
                    
                    uzx_def[index]=zero;
                    uzy_def[index]=zero;
                    uzz_def[index]=zero;
                }

                index++;
			}
}

template< typename TImage >
typename TImage::Pointer 
DownsampleImage0(
                    typename TImage::Pointer inputimage,
                    int scale_x,
                    int scale_y,
                    int scale_z,
                    const std::string interpolationType,
                    const typename TImage::PixelType bkValue = itk::NumericTraits<typename TImage::PixelType>::Zero
                )
{
    typedef itk::ResampleImageFilter< TImage, TImage > ResampleImageType;
    typename ResampleImageType::Pointer ResampleFilter = ResampleImageType::New();

    //set space
    typename TImage::SpacingType oldspacing;
    typename TImage::SpacingType newspacing;
    oldspacing = inputimage->GetSpacing();
    newspacing = oldspacing;
    newspacing[0] *= scale_x;
    newspacing[1] *= scale_y;
    newspacing[2] *= scale_z;
    ResampleFilter->SetOutputSpacing( newspacing );

    //set size
    typename TImage::SizeType   oldsize, newsize;
	oldsize = inputimage->GetLargestPossibleRegion().GetSize();
	newsize = oldsize;
	newsize[0] = (oldsize[0]-1)/scale_x+1; // number of pixels along X
	newsize[1] = (oldsize[1]-1)/scale_y+1; // number of pixels along Y
	newsize[2] = (oldsize[2]-1)/scale_z+1; // number of pixels along Z
    ResampleFilter->SetSize( newsize );

    //set origin
    typename TImage::PointType origin;
	origin = inputimage->GetOrigin();
    ResampleFilter->SetOutputOrigin( origin );

    //set direction
    typename TImage::DirectionType direction;
	direction = inputimage->GetDirection();
    ResampleFilter->SetOutputDirection( direction );

    //set interpolater
    if ( interpolationType == "l" )
    {
        typename itk::LinearInterpolateImageFunction< TImage, double >::Pointer interpolator
            = itk::LinearInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }
    else if ( interpolationType == "n" )
    {
        typename itk::NearestNeighborInterpolateImageFunction< TImage, double >::Pointer interpolator
            = itk::NearestNeighborInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }
    else if ( interpolationType == "b" )
    {
        typename itk::BSplineInterpolateImageFunction< TImage, double >::Pointer interpolator
            =itk::BSplineInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }

    //set transform
	typedef itk::AffineTransform< double, 3 >  TransformType;
	TransformType::Pointer transform = TransformType::New();
    //TransformType::MatrixType matrix = transform->GetMatrix();
    //TransformType::OffsetType offset = transform->GetOffset();
    //std::cout<<" Matrix = "<<std::endl<<matrix<<std::endl;
    //std::cout<<" Offset = "<<std::endl<<offset<<std::endl;

    TransformType::OutputVectorType newtranslation;
	newtranslation[0] = 0;
	newtranslation[1] = 0;
	newtranslation[2] = 0;
    transform->Translate( newtranslation );

    ResampleFilter->SetTransform( transform );
    ResampleFilter->SetDefaultPixelValue( bkValue );
    ResampleFilter->SetInput( inputimage );

    try
    {
        ResampleFilter->UpdateLargestPossibleRegion();
        ResampleFilter->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
        std::cerr<<"Error reading the series "<<std::endl;
        std::cerr<<excp<<std::endl;
        exit(-1);
    }

    typename TImage::Pointer outputimage = ResampleFilter->GetOutput();

    return outputimage;

}



template< typename TImage >
typename TImage::Pointer 
DownsampleDisp0(
                    typename TImage::Pointer inputimage,
                    int scale_x,
                    int scale_y,
                    int scale_z,
                    const std::string interpolationType,
                    const typename TImage::PixelType bkValue = itk::NumericTraits<typename TImage::PixelType>::Zero
                )
{
    typedef itk::ResampleImageFilter< TImage, TImage > ResampleImageType;
    typename ResampleImageType::Pointer ResampleFilter = ResampleImageType::New();

    //set space
    typename TImage::SpacingType oldspacing;
    typename TImage::SpacingType newspacing;
    oldspacing = inputimage->GetSpacing();
    newspacing = oldspacing;
    newspacing[0] *= scale_x;
    newspacing[1] *= scale_y;
    newspacing[2] *= scale_z;
    ResampleFilter->SetOutputSpacing( newspacing );

    //set size
    typename TImage::SizeType   oldsize, newsize;
	oldsize = inputimage->GetLargestPossibleRegion().GetSize();
	newsize = oldsize;
	newsize[0] = (oldsize[0]-1)/scale_x+1; // number of pixels along X
	newsize[1] = (oldsize[1]-1)/scale_y+1; // number of pixels along Y
	newsize[2] = (oldsize[2]-1)/scale_z+1; // number of pixels along Z
    ResampleFilter->SetSize( newsize );

    //set origin
    typename TImage::PointType origin;
	origin = inputimage->GetOrigin();
    ResampleFilter->SetOutputOrigin( origin );

    //set direction
    typename TImage::DirectionType direction;
	direction = inputimage->GetDirection();
    ResampleFilter->SetOutputDirection( direction );

    //set interpolater
    if ( interpolationType == "l" )
    {
        typename itk::LinearInterpolateImageFunction< TImage, double >::Pointer interpolator
            = itk::LinearInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }
    else if ( interpolationType == "n" )
    {
        typename itk::NearestNeighborInterpolateImageFunction< TImage, double >::Pointer interpolator
            = itk::NearestNeighborInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }
    else if ( interpolationType == "b" )
    {
        typename itk::BSplineInterpolateImageFunction< TImage, double >::Pointer interpolator
            =itk::BSplineInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }

    //set transform
	typedef itk::AffineTransform< double, 3 >  TransformType;
	TransformType::Pointer transform = TransformType::New();
    //TransformType::MatrixType matrix = transform->GetMatrix();
    //TransformType::OffsetType offset = transform->GetOffset();
    //std::cout<<" Matrix = "<<std::endl<<matrix<<std::endl;
    //std::cout<<" Offset = "<<std::endl<<offset<<std::endl;

    TransformType::OutputVectorType newtranslation;
	newtranslation[0] = 0;
	newtranslation[1] = 0;
	newtranslation[2] = 0;
    transform->Translate( newtranslation );

    ResampleFilter->SetTransform( transform );
    ResampleFilter->SetDefaultPixelValue( bkValue );
    ResampleFilter->SetInput( inputimage );

    try
    {
        ResampleFilter->UpdateLargestPossibleRegion();
        ResampleFilter->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
        std::cerr<<"Error reading the series "<<std::endl;
        std::cerr<<excp<<std::endl;
        exit(-1);
    }

    typename TImage::Pointer outputimage = ResampleFilter->GetOutput();

    typedef itk::ImageRegionIterator< TImage > IteratorType;
    IteratorType it( outputimage, outputimage->GetRequestedRegion() );

    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
        it.Set( it.Get()/static_cast<float>(scale_x) );
    }

    return outputimage;

}

////////////////////////////////////////////////////////////////////////////////

template< typename TImage >
typename TImage::Pointer 
DownsampleImage(
                    typename TImage::Pointer inputimage,
                    const typename TImage::SizeType size,
                    const typename TImage::SpacingType spacing,
                    const typename TImage::PointType origin,
                    const typename TImage::DirectionType direction,
                    const typename itk::AffineTransform< double, TImage::ImageDimension >::Pointer transform,
                    const std::string interpolationType,
                    const typename TImage::PixelType bkValue = itk::NumericTraits<typename TImage::PixelType>::Zero
                )
{
    typedef itk::ResampleImageFilter< TImage, TImage > ResampleImageType;
    typename ResampleImageType::Pointer ResampleFilter = ResampleImageType::New();

    //set space
    ResampleFilter->SetOutputSpacing( spacing );

    //set size
    ResampleFilter->SetSize( size );

    //set origin
    ResampleFilter->SetOutputOrigin( origin );

    //set direction
    ResampleFilter->SetOutputDirection( direction );

    //set interpolater
    if ( interpolationType == "l" )
    {
        typename itk::LinearInterpolateImageFunction< TImage, double >::Pointer interpolator
            = itk::LinearInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }
    else if ( interpolationType == "n" )
    {
        typename itk::NearestNeighborInterpolateImageFunction< TImage, double >::Pointer interpolator
            = itk::NearestNeighborInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }
    else if ( interpolationType == "b" )
    {
        typename itk::BSplineInterpolateImageFunction< TImage, double >::Pointer interpolator
            =itk::BSplineInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }

    //set transform
    //TransformType::MatrixType matrix = transform->GetMatrix();
    //TransformType::OffsetType offset = transform->GetOffset();
    //std::cout<<" Matrix = "<<std::endl<<matrix<<std::endl;
    //std::cout<<" Offset = "<<std::endl<<offset<<std::endl;

    ResampleFilter->SetTransform( transform );
    ResampleFilter->SetDefaultPixelValue( bkValue );
    ResampleFilter->SetInput( inputimage );

    try
    {
        ResampleFilter->UpdateLargestPossibleRegion();
        ResampleFilter->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
        std::cerr<<"Error reading the series "<<std::endl;
        std::cerr<<excp<<std::endl;
        exit(-1);
    }

    typename TImage::Pointer outputimage = ResampleFilter->GetOutput();

    return outputimage;

}

template< typename TImage >
typename TImage::Pointer 
DownsampleDisp(
                    typename TImage::Pointer inputimage,
                    const typename TImage::SizeType size,
                    const typename TImage::SpacingType spacing,
                    const typename TImage::PointType origin,
                    const typename TImage::DirectionType direction,
                    const typename itk::AffineTransform< double, TImage::ImageDimension >::Pointer transform,
                    const std::string interpolationType,
                    int scale_u,
                    const typename TImage::PixelType bkValue = itk::NumericTraits<typename TImage::PixelType>::Zero
                )
{
    typedef itk::ResampleImageFilter< TImage, TImage > ResampleImageType;
    typename ResampleImageType::Pointer ResampleFilter = ResampleImageType::New();

    //set space
    ResampleFilter->SetOutputSpacing( spacing );

    //set size
    ResampleFilter->SetSize( size );

    //set origin
    ResampleFilter->SetOutputOrigin( origin );

    //set direction
    ResampleFilter->SetOutputDirection( direction );

    //set interpolater
    if ( interpolationType == "l" )
    {
        typename itk::LinearInterpolateImageFunction< TImage, double >::Pointer interpolator
            = itk::LinearInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }
    else if ( interpolationType == "n" )
    {
        typename itk::NearestNeighborInterpolateImageFunction< TImage, double >::Pointer interpolator
            = itk::NearestNeighborInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }
    else if ( interpolationType == "b" )
    {
        typename itk::BSplineInterpolateImageFunction< TImage, double >::Pointer interpolator
            =itk::BSplineInterpolateImageFunction< TImage, double >::New();
        ResampleFilter->SetInterpolator( interpolator );
    }

    //set transform
    //TransformType::MatrixType matrix = transform->GetMatrix();
    //TransformType::OffsetType offset = transform->GetOffset();
    //std::cout<<" Matrix = "<<std::endl<<matrix<<std::endl;
    //std::cout<<" Offset = "<<std::endl<<offset<<std::endl;

    ResampleFilter->SetTransform( transform );
    ResampleFilter->SetDefaultPixelValue( bkValue );
    ResampleFilter->SetInput( inputimage );

    try
    {
        ResampleFilter->UpdateLargestPossibleRegion();
        ResampleFilter->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
        std::cerr<<"Error reading the series "<<std::endl;
        std::cerr<<excp<<std::endl;
        exit(-1);
    }

    typename TImage::Pointer outputimage = ResampleFilter->GetOutput();

    typedef itk::ImageRegionIterator< TImage > IteratorType;
    IteratorType it( outputimage, outputimage->GetRequestedRegion() );

    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
        it.Set( it.Get()/static_cast<float>(scale_u) );
    }


    return outputimage;

}

//deform image
void deformed2_inner_mask(unsigned char* vdeformed, unsigned char* imagein, int xdim, int ydim, int zdim, float* u_x_for, float* u_y_for, float* u_z_for, unsigned char zero)
{
    int planeSize=xdim*ydim;
    float hx, hy, hz;
    int x1, y1, z1;
    int id0;

    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                vdeformed[index]=zero;

				hx=x+u_x_for[index];
				hy=y+u_y_for[index];
				hz=z+u_z_for[index];

                //x1,y1,z1---nearest point; 
				x1=round(hx);
				y1=round(hy);
				z1=round(hz);

                if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
                {
                    id0=z1*planeSize+y1*xdim+x1;
			    	vdeformed[index]=imagein[id0];
                }

                index++;
			}

}
