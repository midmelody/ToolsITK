#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "3DShow.h" 

int main(int argc, char *argv[])
{
  int i,j,k,a,b,t; //counters

 //Read image
  ImageReaderType::Pointer Size1Reader = ImageReaderType::New();
  Size1Reader->SetFileName("ShapeTensor/Size1.img");
  Size1Reader->Update();
  ImageReaderType::Pointer Ori1XReader = ImageReaderType::New();
  Ori1XReader->SetFileName("ShapeTensor/Ori1X.img");
  Ori1XReader->Update(); 
  ImageReaderType::Pointer Ori1YReader = ImageReaderType::New();
  Ori1YReader->SetFileName("ShapeTensor/Ori1Y.img");
  Ori1YReader->Update();
  ImageReaderType::Pointer Ori1ZReader = ImageReaderType::New();
  Ori1ZReader->SetFileName("ShapeTensor/Ori1Z.img");
  Ori1ZReader->Update();
  ImageReaderType::Pointer Size2Reader = ImageReaderType::New();
  Size2Reader->SetFileName("ShapeTensor/Size2.img");
  Size2Reader->Update();
  ImageReaderType::Pointer Ori2XReader = ImageReaderType::New();
  Ori2XReader->SetFileName("ShapeTensor/Ori2X.img");
  Ori2XReader->Update(); 
  ImageReaderType::Pointer Ori2YReader = ImageReaderType::New();
  Ori2YReader->SetFileName("ShapeTensor/Ori2Y.img");
  Ori2YReader->Update();
  ImageReaderType::Pointer Ori2ZReader = ImageReaderType::New();
  Ori2ZReader->SetFileName("ShapeTensor/Ori2Z.img");
  Ori2ZReader->Update();
  ImageReaderType::Pointer Size3Reader = ImageReaderType::New();
  Size3Reader->SetFileName("ShapeTensor/Size3.img");
  Size3Reader->Update();
  ImageReaderType::Pointer Ori3XReader = ImageReaderType::New();
  Ori3XReader->SetFileName("ShapeTensor/Ori3X.img");
  Ori3XReader->Update(); 
  ImageReaderType::Pointer Ori3YReader = ImageReaderType::New();
  Ori3YReader->SetFileName("ShapeTensor/Ori3Y.img");
  Ori3YReader->Update();
  ImageReaderType::Pointer Ori3ZReader = ImageReaderType::New();
  Ori3ZReader->SetFileName("ShapeTensor/Ori3Z.img");
  Ori3ZReader->Update();
  //Get image specs
  ImageType::SpacingType spacing = Size1Reader->GetOutput()->GetSpacing();
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];
  ImageType::SizeType  size = Size1Reader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];
  ImageType::IndexType start = Size1Reader->GetOutput()->GetRequestedRegion().GetIndex();
  ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );


  //Get eigen system for every voxel 
  //Compute projection on XY, YZ and XZ plane
  //Decomposite the projection result  
  TENSOR * tensorXY  = (TENSOR *)malloc(sizeof(TENSOR) * pRow * pCol * pSli); 
  TENSOR * tensorYZ  = (TENSOR *)malloc(sizeof(TENSOR) * pCol * pSli * pRow); 
  TENSOR * tensorXZ  = (TENSOR *)malloc(sizeof(TENSOR) * pRow * pSli * pCol); 

  vnl_matrix<float> xy(3,2);
  xy(0,0)=1; xy(1,0)=0; xy(2,0)=0;
  xy(0,1)=0; xy(1,1)=1; xy(2,1)=0;
  vnl_matrix<float> yz(3,2);
  yz(0,0)=0; yz(1,0)=1; yz(2,0)=0;
  yz(0,1)=0; yz(1,1)=0; yz(2,1)=1;
  vnl_matrix<float> xz(3,2);
  xz(0,0)=1; xz(1,0)=0; xz(2,0)=0;
  xz(0,1)=0; xz(1,1)=0; xz(2,1)=1;
  vnl_matrix<float> eigenVec(3,3);
  eigenVec.fill(0.0);
  vnl_matrix<float> eigenVal(3,3);
  eigenVal.fill(0.0);
  vnl_matrix<float> coVar(3,3);
  coVar.fill(0.0);
  vnl_matrix<float> coVarProXY(2,2);
  coVarProXY.fill(0.0);
  vnl_matrix<float> coVarProYZ(2,2);
  coVarProYZ.fill(0.0);
  vnl_matrix<float> coVarProXZ(2,2);
  coVarProXZ.fill(0.0);

  float alength, blength;
  ImageType::IndexType pixelIndex;


  int specX=78, specY=58, specZ=48;

  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;	  

	  eigenVec(0,0) = Ori1XReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVec(1,0) = Ori1YReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVec(2,0) = Ori1ZReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVal(0,0) = Size1Reader->GetOutput()->GetPixel(pixelIndex);
	  eigenVal(0,0) = eigenVal(0,0) * eigenVal(0,0);

	  eigenVec(0,1) = Ori2XReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVec(1,1) = Ori2YReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVec(2,1) = Ori2ZReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVal(1,1) = Size2Reader->GetOutput()->GetPixel(pixelIndex);

	  eigenVec(0,2) = Ori3XReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVec(1,2) = Ori3YReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVec(2,2) = Ori3ZReader->GetOutput()->GetPixel(pixelIndex);
	  eigenVal(2,2) = Size3Reader->GetOutput()->GetPixel(pixelIndex);

	  /////////////////////////////////////////////////////////////////////////////
	  for(a=-1;a<=0;a++)
	    for(b=-1;b<=0;b++)
	      if((i+a==specX)&&(j+b==specY)&&(k==specZ))
		{
		  std::cout<<a<<" "<<b<<": "<<std::endl;
		  std::cout<<std::endl<<"Shape Tensor Vec:"<<std::endl<<eigenVec<<std::endl;
		  std::cout<<std::endl<<"Shape Tensor Val:"<<std::endl<<eigenVal<<std::endl;
		}
	  /////////////////////////////////////////////////////////////////////////////

	  coVar = eigenVec * eigenVal * vnl_matrix_inverse<float>(eigenVec);
	  coVarProXY = xy.transpose() * coVar * xy;
	  coVarProYZ = yz.transpose() * coVar * yz;
	  coVarProXZ = xz.transpose() * coVar * xz;

	  vnl_symmetric_eigensystem<float> eigXY(coVarProXY);
	  vnl_vector<float> a = eigXY.get_eigenvector(0);
 	  vnl_vector<float> b = eigXY.get_eigenvector(1);
	  alength = eigXY.get_eigenvalue(0);
	  blength = eigXY.get_eigenvalue(1);
	  if(alength>blength)
	    {
	      tensorXY(i,j,k).x = a(0);
	      tensorXY(i,j,k).y = a(1);
	      tensorXY(i,j,k).lengthMajor = alength;
	      tensorXY(i,j,k).lengthMinor = blength;
	    }
	  else
	    {
	      tensorXY(i,j,k).x = b(0);
	      tensorXY(i,j,k).y = b(1);
	      tensorXY(i,j,k).lengthMajor = blength;
	      tensorXY(i,j,k).lengthMinor = alength;
	    }
	  
	  vnl_symmetric_eigensystem<float> eigYZ(coVarProYZ);
	  a = eigYZ.get_eigenvector(0);
 	  b = eigYZ.get_eigenvector(1);
	  alength = eigYZ.get_eigenvalue(0);
	  blength = eigYZ.get_eigenvalue(1);
	  if(alength>blength)
	    {
	      tensorYZ(i,j,k).x = a(0);
	      tensorYZ(i,j,k).y = a(1);
	      tensorYZ(i,j,k).lengthMajor = alength;  
	      tensorYZ(i,j,k).lengthMinor = blength;  
	    }
	  else
	    { 
	      tensorYZ(i,j,k).x = b(0);
	      tensorYZ(i,j,k).y = b(1);
	      tensorYZ(i,j,k).lengthMajor = blength;
	      tensorYZ(i,j,k).lengthMinor = alength;
	    }

	  vnl_symmetric_eigensystem<float> eigXZ(coVarProXZ);
	  a = eigXZ.get_eigenvector(0);
 	  b = eigXZ.get_eigenvector(1);
	  alength = eigXZ.get_eigenvalue(0);
	  blength = eigXZ.get_eigenvalue(1);
	  if(alength>blength)
	    {
	      tensorXZ(i,j,k).x = a(0);
	      tensorXZ(i,j,k).y = a(1);
	      tensorXZ(i,j,k).lengthMajor = alength;  
	      tensorXZ(i,j,k).lengthMinor = blength;  
	    }
	  else
	    {
	      tensorXZ(i,j,k).x = b(0);
	      tensorXZ(i,j,k).y = b(1);
	      tensorXZ(i,j,k).lengthMajor = blength;
	      tensorXZ(i,j,k).lengthMinor = alength;
	    }

	  /*
	  if((i==31)&&(j==66)&&(k==74))
	    {
	      std::cout<<eigenVec<<std::endl;
	      std::cout<<eigenVal<<std::endl;
	      std::cout<<coVar<<std::endl;
	      std::cout<<coVarProXY<<std::endl;
	      std::cout<<coVarProYZ<<std::endl;
	      std::cout<<coVarProXZ<<std::endl;
	      std::cout<<tensorXY(i,j,k).x<<" "<<tensorXY(i,j,k).y<<" "<<tensorXY(i,j,k).lengthMajor<<" "<<tensorXY(i,j,k).lengthMinor<<std::endl;
	      std::cout<<std::endl;
	      std::cout<<tensorYZ(i,j,k).x<<" "<<tensorYZ(i,j,k).y<<" "<<tensorYZ(i,j,k).lengthMajor<<" "<<tensorYZ(i,j,k).lengthMinor<<std::endl;
	      std::cout<<std::endl;
	      std::cout<<tensorXZ(i,j,k).x<<" "<<tensorXZ(i,j,k).y<<" "<<tensorXZ(i,j,k).lengthMajor<<" "<<tensorXZ(i,j,k).lengthMinor<<std::endl;
	      std::cout<<std::endl;
	    }
	  */
	}


  //Assign color coding
  RGBPixelType oneRGBPixel; 
  RGB3DImageType::IndexType RGBIndex;
  RGB3DImageType::IndexType RGBStart;
  RGBStart[0] =   0; 
  RGBStart[1] =   0; 
  RGBStart[2] =   0;
  RGB3DImageType::SizeType  RGBSize;
 
  RGB3DImageType::RegionType regionXY;
  RGBSize[0]  = pRow;  
  RGBSize[1]  = pCol;  
  RGBSize[2]  = pSli;
  regionXY.SetSize( RGBSize );
  regionXY.SetIndex( RGBStart );
  RGB3DImageType::RegionType regionYZ;
  RGBSize[0] = pCol;
  RGBSize[1] = pSli;
  RGBSize[2] = pRow;
  regionYZ.SetSize( RGBSize );
  regionYZ.SetIndex( RGBStart ); 
  RGB3DImageType::RegionType regionXZ;
  RGBSize[0] = pRow;
  RGBSize[1] = pSli;
  RGBSize[2] = pCol;
  regionXZ.SetSize( RGBSize );
  regionXZ.SetIndex( RGBStart );
  RGB3DImageType::Pointer XYImage = RGB3DImageType::New();
  XYImage->SetRegions( regionXY );
  XYImage->Allocate();
  RGB3DImageType::Pointer YZImage = RGB3DImageType::New();
  YZImage->SetRegions( regionYZ );
  YZImage->Allocate();
  RGB3DImageType::Pointer XZImage = RGB3DImageType::New();
  XZImage->SetRegions( regionXZ );
  XZImage->Allocate();
  float H,S,I,RGB[3];
  float major, minor, tanX, tanY;
  float thre = 0.001;
  float maxDT = 20;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  //---------------------XY----------------------------
	  major = tensorXY(i,j,k).lengthMajor;
	  minor = tensorXY(i,j,k).lengthMinor;
	  tanX = tensorXY(i,j,k).x;
	  tanY = tensorXY(i,j,k).y;
	  //Get angle
	  if((tanX>-thre)&&(tanX<thre)) tanX = 0;
	  if((tanY>-thre)&&(tanY<thre)) tanY = 0;
	  if(tanX*tanY>0) 
	    H = atan(tanY/tanX);
	  else if(tanX*tanY<0)
	    H = atan(tanY/tanX)+PI;
	  else if(tanX==0)
	    H = PI/2;
	  else
	    H = 0;
	  //Assign HSI
	  H = H/PI*360;
	  S = sqrt(1.0-minor/major);
	  if(major>maxDT) major = maxDT;
	  if(minor>maxDT) minor = maxDT;
	  I = minor/maxDT;
	  I = sqrt(I);
	  //S = 1;
	  //H = 0;
	  //I = 1;
	  //Convert to RGB and output
	  RGBIndex[0]=i;
	  RGBIndex[1]=j;
 	  RGBIndex[2]=k;
	  HSVtoRGB(RGB, H, S, I);
	  oneRGBPixel[0] = int(RGB[0]*255);  
	  oneRGBPixel[1] = int(RGB[1]*255);
	  oneRGBPixel[2] = int(RGB[2]*255);
	  XYImage->SetPixel( RGBIndex,oneRGBPixel);

	  //---------------------YZ----------------------------
	  major = tensorYZ(i,j,k).lengthMajor;
	  minor = tensorYZ(i,j,k).lengthMinor;
	  tanX = tensorYZ(i,j,k).x;
	  tanY = tensorYZ(i,j,k).y;
	  //Get angle
	  if((tanX>-thre)&&(tanX<thre)) tanX = 0;
	  if((tanY>-thre)&&(tanY<thre)) tanY = 0;
	  if(tanX*tanY>0) 
	    H = atan(tanY/tanX);
	  else if(tanX*tanY<0)
	    H = atan(tanY/tanX)+PI;
	  else if(tanX==0)
	    H = PI/2;
	  else
	    H = 0;
	  //Assign HSI
	  H = H/PI*360;
	  S = sqrt(1.0-minor/major);
	  if(major>maxDT) major = maxDT;
	  if(minor>maxDT) minor = maxDT;
	  I = minor/maxDT;
	  I = sqrt(I);
	  //S = 1;
	  //H = 0;
	  //I = 1;
	  //Convert to RGB and output
	  RGBIndex[0]=j;
	  RGBIndex[1]=k;
 	  RGBIndex[2]=i;
	  HSVtoRGB(RGB, H, S, I);
	  oneRGBPixel[0] = int(RGB[0]*255);  
	  oneRGBPixel[1] = int(RGB[1]*255);
	  oneRGBPixel[2] = int(RGB[2]*255);
	  YZImage->SetPixel( RGBIndex,oneRGBPixel);

	  //---------------------XZ----------------------------
	  major = tensorXZ(i,j,k).lengthMajor;
	  minor = tensorXZ(i,j,k).lengthMinor;
	  tanX = tensorXZ(i,j,k).x;
	  tanY = tensorXZ(i,j,k).y;
	  //Get angle
	  if((tanX>-thre)&&(tanX<thre)) tanX = 0;
	  if((tanY>-thre)&&(tanY<thre)) tanY = 0;
	  if(tanX*tanY>0) 
	    H = atan(tanY/tanX);
	  else if(tanX*tanY<0)
	    H = atan(tanY/tanX)+PI;
	  else if(tanX==0)
	    H = PI/2;
	  else
	    H = 0;
	  //Assign HSI
	  H = H/PI*360;
	  S = sqrt(1.0-minor/major);
	  if(major>maxDT) major = maxDT;
	  if(minor>maxDT) minor = maxDT;
	  I = minor/maxDT;
	  I = sqrt(I);
	  //S = 1;
	  //H = 0;
	  //I = 1;
	  //Convert to RGB and output
	  RGBIndex[0]=i;
	  RGBIndex[1]=k;
 	  RGBIndex[2]=j;
	  HSVtoRGB(RGB, H, S, I);
	  oneRGBPixel[0] = int(RGB[0]*255);  
	  oneRGBPixel[1] = int(RGB[1]*255);
	  oneRGBPixel[2] = int(RGB[2]*255);
	  XZImage->SetPixel( RGBIndex,oneRGBPixel);
	}

  typedef itk::NumericSeriesFileNames    NameGeneratorType;
  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
  nameGenerator->SetIncrementIndex( 1 );

  SeriesWriterType::Pointer XYWriter = SeriesWriterType::New();
  XYWriter->SetInput( XYImage );
  nameGenerator->SetStartIndex( 0 );
  nameGenerator->SetEndIndex( pSli-1 );
  nameGenerator->SetSeriesFormat( "XY/output%03d.png" );
  XYWriter->SetFileNames( nameGenerator->GetFileNames() );
  XYWriter->Update();

  SeriesWriterType::Pointer YZWriter = SeriesWriterType::New();
  YZWriter->SetInput( YZImage );
  nameGenerator->SetStartIndex( 0 );
  nameGenerator->SetEndIndex( pRow-1 );
  nameGenerator->SetSeriesFormat( "YZ/output%03d.png" );
  YZWriter->SetFileNames( nameGenerator->GetFileNames() );
  YZWriter->Update();

  SeriesWriterType::Pointer XZWriter = SeriesWriterType::New();
  XZWriter->SetInput( XZImage );
  nameGenerator->SetStartIndex( 0 );
  nameGenerator->SetEndIndex( pCol-1 );
  nameGenerator->SetSeriesFormat( "XZ/output%03d.png" );
  XZWriter->SetFileNames( nameGenerator->GetFileNames() );
  XZWriter->Update();

  free(tensorXY);
  free(tensorYZ);
  free(tensorXZ);
  return EXIT_SUCCESS;
}







	  /*	  
	  if((i==54)&&(j==22)&&(k==74))
	    {
	      std::cout<<eigenVec<<std::endl;
	      std::cout<<eigenVal<<std::endl;
	      std::cout<<coVar<<std::endl;
	      std::cout<<coVarProXY<<std::endl;
	      std::cout<<coVarProYZ<<std::endl;
	      std::cout<<coVarProXZ<<std::endl;
	      std::cout<<tensorXY(i,j,k).x<<" "<<tensorXY(i,j,k).y<<" "<<tensorXY(i,j,k).lengthMajor<<" "<<tensorXY(i,j,k).lengthMinor<<std::endl;
	      std::cout<<std::endl;
	      std::cout<<tensorYZ(i,j,k).x<<" "<<tensorYZ(i,j,k).y<<" "<<tensorYZ(i,j,k).lengthMajor<<" "<<tensorYZ(i,j,k).lengthMinor<<std::endl;
	      std::cout<<std::endl;
	      std::cout<<tensorXZ(i,j,k).x<<" "<<tensorXZ(i,j,k).y<<" "<<tensorXZ(i,j,k).lengthMajor<<" "<<tensorXZ(i,j,k).lengthMinor<<std::endl;
	      std::cout<<std::endl;
	    }
	  */

