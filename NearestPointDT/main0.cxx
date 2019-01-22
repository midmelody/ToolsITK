#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <stdlib.h>
#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " distanceImage binaryImage outputLinkImage outputRoIImage";
      std::cerr << " " << std::endl;
      return EXIT_FAILURE;
    }

  typedef float PixelType;
  typedef unsigned char OutPixelType;
  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::Image< OutPixelType, 3 > OutImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< OutImageType >  WriterType;
  ReaderType::Pointer readerDT = ReaderType::New();
  ReaderType::Pointer readerSeg = ReaderType::New();
  WriterType::Pointer writerLK = WriterType::New();
  WriterType::Pointer writerROI = WriterType::New();
  readerDT->SetFileName( argv[1] );
  readerSeg->SetFileName( argv[2] );
  writerLK->SetFileName( argv[3] );
  writerROI->SetFileName( argv[4] );

  try
    {
      readerDT->Update();
      readerSeg->Update();
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem reading image file : " << argv[1] << std::endl;
      std::cerr << excp << std::endl;
      return -1;
    }

  int i,j,k; //counters

  //Image specs
  int pRow,pCol,pSli; //Size of Image: row, column, slice
  int spaRow, spaCol, spaSli; //Spacing: row, column, slice
  OutImageType::SpacingType spacing = readerDT->GetOutput()->GetSpacing();
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];
  OutImageType::PointType origin = readerDT->GetOutput()->GetOrigin(); 
  OutImageType::DirectionType direction = readerDT->GetOutput()->GetDirection();
  OutImageType::SizeType  size = readerDT->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  OutImageType::RegionType region;
  OutImageType::IndexType start = readerDT->GetOutput()->GetRequestedRegion().GetIndex();
  region.SetSize( size );
  region.SetIndex( start );
  OutImageType::IndexType pixelIndex;

  //Find nearest in airway
  float min = 1000000;
  int label;
  float dist;
  int curI, curJ, curK;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  dist = readerDT->GetOutput()->GetPixel(pixelIndex);
	  label = readerSeg->GetOutput()->GetPixel(pixelIndex);
	  if(label==2)
	    if(dist<min)
	      {
		min = dist;
		curI = i;
		curJ = j;
		curK = k;
	      }
	}
  std::cout<<curI<<" "<<curJ<<" "<<curK<<": "<<min<<std::endl;

  //Find direction at (curI, curJ, curK) - Sobel operator
  float dx, dy, dz;
  float posX,negX,posY,negY,posZ,negZ; 
  float weight;
  float norm;
  float DTValue; 
  int a,b;
  posX = 0;
  negX = 0;
  posY = 0;
  negY = 0;
  posZ = 0;
  negZ = 0;
  for(a=-1;a<=1;a++)
    for(b=-1;b<=1;b++)
      {
	if((abs(a)+abs(b))==0) weight = 4.0;
	if((abs(a)+abs(b))==1) weight = 2.0;
	if((abs(a)+abs(b))==2) weight = 1.0;
	pixelIndex[0]=curI+1;
	pixelIndex[1]=curJ+a;
	pixelIndex[2]=curK+b;
	DTValue = readerDT->GetOutput()->GetPixel(pixelIndex);
	posX = posX + DTValue*weight;

	pixelIndex[0]=curI-1;
	pixelIndex[1]=curJ+a;
	pixelIndex[2]=curK+b;
	DTValue = readerDT->GetOutput()->GetPixel(pixelIndex);
	negX = negX + DTValue*weight;

	pixelIndex[0]=curI+a;
	pixelIndex[1]=curJ+1;
	pixelIndex[2]=curK+b;
	DTValue = readerDT->GetOutput()->GetPixel(pixelIndex);
	posY = posY + DTValue*weight;

	pixelIndex[0]=curI+a;
	pixelIndex[1]=curJ-1;
	pixelIndex[2]=curK+b;
	DTValue = readerDT->GetOutput()->GetPixel(pixelIndex);
	negY = negY + DTValue*weight;

	pixelIndex[0]=curI+a;
	pixelIndex[1]=curJ+b;
	pixelIndex[2]=curK+1;
	DTValue = readerDT->GetOutput()->GetPixel(pixelIndex);
	posZ = posZ + DTValue*weight;

	pixelIndex[0]=curI+a;
	pixelIndex[1]=curJ+b;
	pixelIndex[2]=curK-1;
	DTValue = readerDT->GetOutput()->GetPixel(pixelIndex);
	negZ = negZ + DTValue*weight;
      }
  dx = (posX-negX)/32.0;
  dy = (posY-negY)/32.0;
  dz = (posZ-negZ)/32.0;
  norm = dx*dx + dy*dy + dz*dz;
  norm = sqrt(norm);
  dx = dx/norm;
  dy = dy/norm;
  dz = dz/norm;

  //Output
  OutImageType::Pointer linkImage = OutImageType::New();
  linkImage->SetRegions( region );
  linkImage->SetSpacing( spacing );
  linkImage->SetOrigin( origin );
  linkImage->SetDirection( direction );
  linkImage->Allocate();

  OutImageType::Pointer ROIImage = OutImageType::New();
  ROIImage->SetRegions( region );
  ROIImage->SetSpacing( spacing );
  ROIImage->SetOrigin( origin );
  ROIImage->SetDirection( direction );
  ROIImage->Allocate();

  int l;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  linkImage->SetPixel(pixelIndex,readerSeg->GetOutput()->GetPixel(pixelIndex));
	  ROIImage->SetPixel(pixelIndex,0);
	  dist = readerDT->GetOutput()->GetPixel(pixelIndex);
	  if((dist>=-2)&(dist<=min+10))
	    ROIImage->SetPixel(pixelIndex,1);
	  for(l=0;l<min+1;l++)
	    {
	      pixelIndex[0]=round(curI - float(l)*dx);
	      pixelIndex[1]=round(curJ - float(l)*dy);
	      pixelIndex[2]=round(curK - float(l)*dz);
	      linkImage->SetPixel(pixelIndex,3);
	    }	    
	}
  writerLK->SetInput( linkImage );
  writerLK->Update();
  writerROI->SetInput( ROIImage );
  writerROI->Update();
  return EXIT_SUCCESS;
}

  /*
  //Get bounding box
  int Rmax, Rmin, Cmax, Cmin, Smax, Smin;
  Rmax = 0;
  Cmax = 0;
  Smax = 0;
  Rmin = 100000;
  Cmin = 100000;
  Smin = 100000;
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow = size[0];
  int pCol = size[1];
  int pSli = size[2]; 
  
  int i,j,k,value;
  ImageType::IndexType pixelIndex;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  value = readerSeg->GetOutput()->GetPixel(pixelIndex);
	  if(value)
	    {
	      if(Rmax<i)
		Rmax = i;
	      if(Rmin>i)
		Rmin = i;
	      if(Cmax<j)
		Cmax = j;
	      if(Cmin>j)
		Cmin = j;
	      if(Smax<k)
		Smax = k;
	      if(Smin>k)
		Smin = k;
	    }
	}

  //Extract VOI
  int lRow, lCol, lSli; //low
  int hRow, hCol, hSli; //high
  int wRow, wCol, wSli; //width
  lRow = 0;
  lCol = 0;
  lSli = 0;
  hRow = pRow;
  hCol = pCol;
  hSli = pSli;
  wRow = Rmax - Rmin;
  wCol = Cmax - Cmin;
  wSli = Smax - Smin;
  if(Rmin - wRow > lRow)
    lRow = Rmin - wRow;
  if(Rmax + wRow < hRow)
    hRow = Rmax + wRow;
  if(Cmin - wCol > lCol)
    lCol = Cmin - wCol;
  if(Cmax + wCol < hCol)
    hCol = Cmax + wCol;
  if(Smin - wSli > lSli)
    lSli = Smin - wSli;
  if(Smax + wSli < hSli)
    hSli = Smax + wSli;

  //Resize image
  size[0] = hRow - lRow;
  size[1] = hCol - lCol;
  size[2] = hSli - lSli;
  //Read image
  ImageType::RegionType region;
  region.SetSize( size );
  ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions( region );
  outImage->SetRegions( size );
  outImage->SetSpacing( spacing );
  outImage->SetOrigin( origin );
  outImage->SetDirection( direction );
  outImage->Allocate();
  outImage->FillBuffer( 0.0 );
  ImageType::Pointer outImageSeg = ImageType::New();
  outImageSeg->SetRegions( region );
  outImageSeg->SetRegions( size );
  outImageSeg->SetSpacing( spacing );
  outImageSeg->SetOrigin( origin );
  outImageSeg->SetDirection( direction );
  outImageSeg->Allocate();
  outImageSeg->FillBuffer( 0.0 );
  ImageType::PixelType pixel;
  ImageType::IndexType inputIndex;
  ImageType::IndexType outputIndex;
	
  for( outputIndex[2]=0, inputIndex[2]=lSli; outputIndex[2]<size[2], inputIndex[2]<hSli; outputIndex[2]++, inputIndex[2]++ )
    {
      for( outputIndex[1]=0, inputIndex[1]=lCol; outputIndex[1]<size[1], inputIndex[1]<hCol; outputIndex[1]++, inputIndex[1]++ )
	{
	  for( outputIndex[0]=0, inputIndex[0]=lRow; outputIndex[0]<size[0], inputIndex[0]<hRow; outputIndex[0]++, inputIndex[0]++ )
	    {
	      pixel = reader->GetOutput()->GetPixel( inputIndex );
	      outImage->SetPixel( outputIndex, pixel );
	      pixel = readerSeg->GetOutput()->GetPixel( inputIndex );
	      outImageSeg->SetPixel( outputIndex, pixel );
	    }
	}
    }

  writer->SetInput( outImage );
  writer->Update();
  writerSeg->SetInput( outImageSeg );
  writerSeg->Update();

  */


