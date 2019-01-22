#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 

struct VoxelType
{
  float x;
  float y;
  float z;
};

//Branch struct 
struct BranchType
{
  int voxelCt;
  int label;
  int parLabel;
  VoxelType * voxelIndex;
  VoxelType * voxelGradient;
};

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputImageFile outputTreeFile outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  //Read tree
  VoxelType index;
  int label, parLabel;
  //find total branch number
  FILE *treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  int branchCt = 0;
  int curLabel = 0;
  while(fscanf(treeFile, "%f %f %f %d %d", &index.x, &index.y, &index.z, &label, &parLabel) == 5)
    {
      if(label > curLabel)
	{	  
	  curLabel = label;
	}
    }
  branchCt = curLabel;
  fclose(treeFile);

  //get information
  //index in branchList: branchLabel-1
  BranchType * branchList = new BranchType[branchCt];
  int a;
  for(a=0;a<branchCt;a++)
    branchList[a].voxelCt = 0;
  treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(treeFile, "%f %f %f %d %d", &index.x, &index.y, &index.z, &label, &parLabel) == 5)
    {
      branchList[label-1].voxelCt++;
      branchList[label-1].label = label;
      branchList[label-1].parLabel = parLabel;      
    }
  fclose(treeFile);

  //Record point index
  int curVoxel[branchCt];
  for(a=0;a<branchCt;a++)
    {
      if(branchList[a].voxelCt>0)
	{
	  branchList[a].voxelIndex = new VoxelType[branchList[a].voxelCt];
	  branchList[a].voxelGradient = new VoxelType[branchList[a].voxelCt];
	  curVoxel[a] = 0;
	}
    }

  treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(treeFile, "%f %f %f %d %d", &index.x, &index.y, &index.z, &label, &parLabel) == 5)
    {
      branchList[label-1].voxelIndex[curVoxel[label-1]] = index;
      branchList[label-1].voxelGradient[curVoxel[label-1]].x = 0;
      branchList[label-1].voxelGradient[curVoxel[label-1]].y = 0;
      branchList[label-1].voxelGradient[curVoxel[label-1]].z = 0;
      curVoxel[label-1]++;
    }
  fclose(treeFile);

  int b;
  //Compute tangent with least square over whole branch using SVD
  for(a=0;a<branchCt;a++)
    {
      if(branchList[a].voxelCt>0)
	{
	  float meanX = 0, meanY = 0, meanZ = 0;
	  for(b=0;b<branchList[a].voxelCt;b++)
	    {
	      meanX = meanX + branchList[a].voxelIndex[b].x;
	      meanY = meanY + branchList[a].voxelIndex[b].y;
	      meanZ = meanZ + branchList[a].voxelIndex[b].z;
	    }
	  meanX = meanX/float(branchList[a].voxelCt);
	  meanY = meanY/float(branchList[a].voxelCt);
	  meanZ = meanZ/float(branchList[a].voxelCt);
  
	  vnl_matrix<double> points(branchList[a].voxelCt,3,0.0);
	  for(b=0;b<branchList[a].voxelCt;b++)
	    {
	      points(b,0) = branchList[a].voxelIndex[b].x - meanX;
	      points(b,1) = branchList[a].voxelIndex[b].y - meanY;
	      points(b,2) = branchList[a].voxelIndex[b].z - meanZ;
	    }
	  vnl_matrix<double> SVDV = vnl_svd<double> (points).V(); 
	  for(b=0;b<branchList[a].voxelCt;b++)
	    {     
	      branchList[a].voxelGradient[b].x = SVDV(0,0);
	      branchList[a].voxelGradient[b].y = SVDV(1,0);
	      branchList[a].voxelGradient[b].z = SVDV(2,0);
	    }
	}
    }

  //Output file with tangent vector
  FILE *outFile = fopen(argv[3], "w+");
  for(a=0;a<branchCt;a++)
    {
      for(b = 0; b<branchList[a].voxelCt; b++)
	fprintf(outFile, "%f %f %f %d %d %f %f %f\n", branchList[a].voxelIndex[b].x, branchList[a].voxelIndex[b].y, branchList[a].voxelIndex[b].z, branchList[a].label, branchList[a].parLabel, branchList[a].voxelGradient[b].x, branchList[a].voxelGradient[b].y, branchList[a].voxelGradient[b].z);
      //fprintf(outFile, "%f %f %f %d %d\n", branchList[a].voxelIndex[b].x, branchList[a].voxelIndex[b].y, branchList[a].voxelIndex[b].z, branchList[a].label, branchList[a].parLabel);
    }
  fclose(outFile);


  //ITK settings
  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //Parameters
  reader->SetFileName( argv[2] );
  writer->SetFileName( argv[4] );
 
  //Pipeline
  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate image for output
  ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions( region );
  outImage->SetSpacing( spacing );
  outImage->SetOrigin( origin );
  outImage->SetDirection( direction );
  outImage->Allocate();  
  //Initialize
  ImageType::IndexType pixelIndex;
  int i, j, k;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  outImage->SetPixel(pixelIndex, reader->GetOutput()->GetPixel(pixelIndex));
	}

  //Write tree and corresponding gradient vector
  
  int c,d;
  /*
  for(a=0;a<branchCt;a++)
    {
      for(b=0;b<branchList[a].voxelCt;b++)
	{
	  pixelIndex[0]=int(branchList[a].voxelIndex[b].x);
	  pixelIndex[1]=int(branchList[a].voxelIndex[b].y);
	  pixelIndex[2]=int(branchList[a].voxelIndex[b].z);    
	  outImage->SetPixel(pixelIndex,2);
	}
    }
  */
  a=2;
  for(c=70;c<branchList[a].voxelCt-2;c=c+branchList[a].voxelCt)
    {
      /*
      //Show tangent line
      for(d=1;d<=10;d=d++)
	{
	  pixelIndex[0]=int(branchList[a].voxelIndex[c].x + d*branchList[a].voxelGradient[c].x);
	  pixelIndex[1]=int(branchList[a].voxelIndex[c].y + d*branchList[a].voxelGradient[c].y);
	  pixelIndex[2]=int(branchList[a].voxelIndex[c].z + d*branchList[a].voxelGradient[c].z);  
	  if((pixelIndex[0]>=0)&&(pixelIndex[0]<pRow)&&(pixelIndex[1]>=0)&&(pixelIndex[1]<pCol)&&(pixelIndex[2]>=0)&&(pixelIndex[2]<pSli))  
	    outImage->SetPixel(pixelIndex, 3);
	}
      */
      float xT = branchList[a].voxelGradient[c].x;
      float yT = branchList[a].voxelGradient[c].y;
      float zT = branchList[a].voxelGradient[c].z;
      //Compute Rotation system
      float sin1, sin2, cos1, cos2;
      if((yT!=0)||(zT!=0))
	{
	  sin1 = yT/sqrt(yT*yT+zT*zT);
	  cos1 = zT/sqrt(yT*yT+zT*zT);
	  sin2 = xT/sqrt(xT*xT+yT*yT+zT*zT);
	  cos2 = sqrt(yT*yT+zT*zT)/sqrt(xT*xT+yT*yT+zT*zT);
	}
      else
	{
	  sin1 = 1;
	  cos1 = 0;
	  sin2 = 1;
	  cos2 = 0;
	}
      //Show plane
      for(int x=-20; x<=20; x++)
	for(int y=-20; y<=20; y++)
	  {
	    //Get coordinate after rotation
	    float xG = cos2*float(x);
	    float yG = -sin1*sin2*float(x)+cos1*float(y);
	    float zG = -cos1*sin2*float(x)-sin1*float(y);
	    //std::cout<<xG<<" "<<yG<<" "<<zG<<std::endl;
	    xG = xG + branchList[a].voxelIndex[c].x;
	    yG = yG + branchList[a].voxelIndex[c].y;
	    zG = zG + branchList[a].voxelIndex[c].z;

	    pixelIndex[0] = int(xG);
	    pixelIndex[1] = int(yG);
	    pixelIndex[2] = int(zG);
	    if((pixelIndex[0]>=0)&&(pixelIndex[0]<pRow)&&(pixelIndex[1]>=0)&&(pixelIndex[1]<pCol)&&(pixelIndex[2]>=0)&&(pixelIndex[2]<pSli))  
	      outImage->SetPixel(pixelIndex, 2);
	  }
    }
  
  //Write output  
  writer->SetInput( outImage );
  try
    {
      writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
