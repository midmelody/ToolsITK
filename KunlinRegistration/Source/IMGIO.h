/*
 * =====================================================================================
 *
 *       Filename:  ImgIO.h
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  10/25/2009 03:06:10 PM CDT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */

#ifndef IMGIO_H
#define IMGIO_H

#include "stdlib.h"
#include "stdio.h"
#include <string>
#include "ImageType.h"

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"


using namespace std;
class IMGIO
{
public:
	inline IMGIO(){};
	inline ~IMGIO(){};


    template < typename TImage > 
        static typename TImage::Pointer ReadImage( std::string filename );

	template < typename TImage >
		static bool  WriteImage( const typename TImage::Pointer image, std::string filename ); 
	
    template < typename TImageOriginal, typename TImageWrite >
		static bool  CastWriteImage( const typename TImageOriginal::Pointer image, std::string filename );

    template < typename TImageOriginal, typename TImageWrite >
		static typename TImageWrite::Pointer  CastImage( const typename TImageOriginal::Pointer image );

    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /*template < typename TImage > 
		static bool ReadImage( typename TImage::Pointer& image, const char* filename );

	template < typename TImage >
		static bool  WriteImage( const typename TImage::Pointer image, const char* filename ); 
	
    template < typename TImageOriginal, typename TImageWrite >
		static bool  CastWriteImage( const typename TImageOriginal::Pointer image, const char* filename );*/

    ///////////////////////////////////////////////////////////////////////////////////////////////////

    template < typename TImage > 
        static typename TImage::Pointer initImage( typename TImage::IndexType iStart, typename TImage::SizeType size, typename TImage::DirectionType direction, typename TImage::PointType origin, typename TImage::SpacingType spacing );

    template < typename TImage >
        static typename TImage::Pointer initmhd( typename TImage::IndexType iStart, typename TImage::SizeType size, typename TImage::DirectionType direction, typename TImage::PointType origin, typename TImage::SpacingType spacing );

    template < typename TImage > 
        static typename TImage::Pointer initImageWithValue( typename TImage::IndexType iStart, typename TImage::SizeType size, typename TImage::DirectionType direction, typename TImage::PointType origin, typename TImage::SpacingType spacing, typename TImage::PixelType pixel );

};



//////////////////////////////////////////////////////////////////////

template < typename TImage > 
typename TImage::Pointer IMGIO::ReadImage( std::string filename )
{
	typedef itk::ImageFileReader< TImage > ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( filename );

    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "ReadImage exception caught ! --- " << filename << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return reader->GetOutput();
}

template < typename TImage >
bool IMGIO::WriteImage( const typename TImage::Pointer image, std::string filename )
{
    typedef itk::ImageFileWriter< TImage > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( filename );
	writer->SetInput( image );
    
    try
    {
        writer->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "WriteImage exception caught ! --- " << filename << std::endl;
        std::cerr << err << std::endl;
        return false;
    }
    
    return true;
}

template < typename TImageOriginal, typename TImageWrite >
bool IMGIO::CastWriteImage( const typename TImageOriginal::Pointer image, std::string filename )
{
    typedef itk::CastImageFilter< TImageOriginal, TImageWrite > CastFilterType;
    typename CastFilterType::Pointer caster = CastFilterType::New();
    caster->SetInput( image );
    caster->Update();

    typedef itk::ImageFileWriter< TImageWrite > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( filename );
	writer->SetInput( caster->GetOutput() );
    
    try
    {
        writer->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "CastWriteImage exception caught ! --- " << filename << std::endl;
        std::cerr << err << std::endl;
        return false;
    }
    
    return true;
}

template < typename TImageOriginal, typename TImageWrite >
typename TImageWrite::Pointer IMGIO::CastImage( const typename TImageOriginal::Pointer image )
{
    typedef itk::CastImageFilter< TImageOriginal, TImageWrite > CastFilterType;
    typename CastFilterType::Pointer caster = CastFilterType::New();
    caster->SetInput( image );

    try
    {
        caster->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "CastImage exception caught ! --- "  << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }
    
    return caster->GetOutput();
}


//////////////////////////////////////////////////////////////////////////////////////////

/*template < typename TImage > 
bool IMGIO::ReadImage( typename TImage::Pointer& image, const char* filename )
{
	typedef itk::ImageFileReader< TImage > ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( filename );

    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "ReadImage exception caught ! --- " << filename << std::endl;
        std::cerr << err << std::endl;
        return false;
    }

    image = reader->GetOutput();
	return true;
}

template < typename TImage >
bool IMGIO::WriteImage( const typename TImage::Pointer image, const char* filename )
{
    typedef itk::ImageFileWriter< TImage > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( filename );
	writer->SetInput( image );
    
    try
    {
        writer->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "WriteImage exception caught ! --- " << filename << std::endl;
        std::cerr << err << std::endl;
        return false;
    }
    
    return true;
}

template < typename TImageOriginal, typename TImageWrite >
bool IMGIO::CastWriteImage( const typename TImageOriginal::Pointer image, const char* filename )
{
    typedef itk::CastImageFilter< TImageOriginal, TImageWrite > CastFilterType;
    typename CastFilterType::Pointer caster = CastFilterType::New();
    caster->SetInput( image );
    caster->Update();

    typedef itk::ImageFileWriter< TImageWrite > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( filename );
	writer->SetInput( caster->GetOutput() );
    
    try
    {
        writer->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "CastWriteImage exception caught ! --- " << filename << std::endl;
        std::cerr << err << std::endl;
        return false;
    }
    
    return true;
}*/

//////////////////////////////////////////////////////////////////////////////////////////

template < typename TImage > 
typename TImage::Pointer IMGIO::initImage( typename TImage::IndexType iStart, typename TImage::SizeType size, typename TImage::DirectionType direction, typename TImage::PointType origin, typename TImage::SpacingType spacing )
{
	typename TImage::RegionType region;
	region.SetIndex(iStart);
	region.SetSize(size);
	
	typename TImage::Pointer image = TImage::New();
    //InpuTImage->SetLargestPossibleRegion( Region );
	//InpuTImage->SetBufferedRegion( Region );
	//InpuTImage->SetRequestedRegion( Region );
	image->SetRegions( region );
	image->Allocate();
	
    image->SetDirection( direction );
    image->SetOrigin( origin );
	image->SetSpacing( spacing );
	
    /*typename itk::ImageRegionIterator< TImage>  imageIt( image, image->GetRequestedRegion() );
	for( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
	{
		imageIt.Set(0);
	}*/
    image->FillBuffer(0);
	
	return image;
}

template < typename TImage >
typename TImage::Pointer IMGIO::initmhd( typename TImage::IndexType iStart, typename TImage::SizeType size, typename TImage::DirectionType direction, typename TImage::PointType origin, typename TImage::SpacingType spacing )
{
	typename TImage::RegionType region;
	region.SetIndex(iStart);
	region.SetSize(size);
	
	typename TImage::Pointer image = TImage::New();
    //InpuTImage->SetLargestPossibleRegion( Region );
	//InpuTImage->SetBufferedRegion( Region );
	//InpuTImage->SetRequestedRegion( Region );
	image->SetRegions( region );
	image->Allocate();
	
    image->SetDirection( direction );
    image->SetOrigin( origin );
	image->SetSpacing( spacing );
	
    typename itk::ImageRegionIterator< TImage>  imageIt( image, image->GetRequestedRegion() );
    typename TImage::PixelType pixel;
	for( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
	{
        pixel[0]=0;
        pixel[1]=0;
        pixel[2]=0;
		imageIt.Set(pixel);
	}
    //image->FillBuffer(0);
	
	return image;
}

template < typename TImage > 
typename TImage::Pointer IMGIO::initImageWithValue( typename TImage::IndexType iStart, typename TImage::SizeType size, typename TImage::DirectionType direction, typename TImage::PointType origin, typename TImage::SpacingType spacing, typename TImage::PixelType initPixelValue )
{
	typename TImage::RegionType region;
	region.SetIndex(iStart);
	region.SetSize(size);
	
	typename TImage::Pointer image = TImage::New();
	image->SetRegions( region );
	image->Allocate();
	
    image->SetDirection( direction );
    image->SetOrigin( origin );
	image->SetSpacing( spacing );
	
    image->FillBuffer( initPixelValue );
	
	return image;
}


extern IMGIO ImageIO;



#endif
