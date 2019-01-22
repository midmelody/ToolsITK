#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMedianImageFilter.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "stdio.h" 
#define PI 3.141592654
#define sampleSize 10
#define angleSize 360

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
  VoxelType * voxelTangent;
};

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeFile inputImageFile outputImageFile airwayMaskFile" << std::endl;
      return EXIT_FAILURE;
    }
  int i, j;
  //Read tree
  VoxelType index;
  VoxelType tangent;
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
  while(fscanf(treeFile, "%f %f %f %d %d %f %f %f", &index.x, &index.y, &index.z, &label, &parLabel, &tangent.x, &tangent.y, &tangent.z) == 8)
    {
      if(label > curLabel)
	{
	  branchCt++;
	  curLabel = label;
	}
    }
  fclose(treeFile);

  //get information
  //index in branchList: branchLabel-1
  //branchCt+1: for adding missed part
  BranchType * branchList = new BranchType[branchCt+1];
  int a, b;
  for(a=0;a<branchCt;a++)
    branchList[a].voxelCt = 0;
  treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }

  while(fscanf(treeFile, "%f %f %f %d %d %f %f %f", &index.x, &index.y, &index.z, &label, &parLabel, &tangent.x, &tangent.y, &tangent.z) == 8)
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
      branchList[a].voxelIndex = new VoxelType[branchList[a].voxelCt];
      branchList[a].voxelTangent = new VoxelType[branchList[a].voxelCt];
      curVoxel[a] = 0;
    }
  treeFile = fopen(argv[1], "r");
  if(treeFile == NULL) 
    {
      std::cerr << argv[1] << " Tree file doesn't exists"<< std::endl;
      return EXIT_FAILURE;
    }
  while(fscanf(treeFile, "%f %f %f %d %d %f %f %f", &index.x, &index.y, &index.z, &label, &parLabel, &tangent.x, &tangent.y, &tangent.z) == 8)
    {
      branchList[label-1].voxelIndex[curVoxel[label-1]] = index;
      branchList[label-1].voxelTangent[curVoxel[label-1]] = tangent;
      curVoxel[label-1]++;
    }
  fclose(treeFile);

  //Resample original image at a point along centerline on orthogonal plane to tangent vector
  //ITK settings
  //Read: 3D original image
  //Write: 3D result
  const unsigned int Dimension = 3;
  typedef int PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writerInner = WriterType::New();
  WriterType::Pointer writerOuter = WriterType::New();
  ReaderType::Pointer readerMask = ReaderType::New();
  //Parameters
  reader->SetFileName( argv[2] );
  readerMask->SetFileName( argv[4] );
  //Pipeline
  try
    {
      reader->Update();
      readerMask->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  //Get in image specs
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  ImageType::Pointer imageOutInner = ImageType::New();
  imageOutInner->SetRegions( region );
  imageOutInner->SetSpacing( spacing );
  imageOutInner->SetOrigin( origin );
  imageOutInner->SetDirection( direction );
  imageOutInner->Allocate();  
  imageOutInner->FillBuffer(0);
  ImageType::Pointer imageOutOuter = ImageType::New();
  imageOutOuter->SetRegions( region );
  imageOutOuter->SetSpacing( spacing );
  imageOutOuter->SetOrigin( origin );
  imageOutOuter->SetDirection( direction );
  imageOutOuter->Allocate();  
  imageOutOuter->FillBuffer(0);

  ImageType::IndexType pixelIndex;


  //Add missed part in trachea
  int zExtend = size[2]-1;
  int xExtend = 0;
  int yExtend = 0;
  bool statZero = true;
  while(statZero)
    {
      int addTotal = 0;
      for(int i=0; i<size[0]; i++)
	for(int j=0; j<size[1]; j++)
	  {
	    pixelIndex[0] = i;
	    pixelIndex[1] = j;
	    pixelIndex[2] = zExtend;
	    if(readerMask->GetOutput()->GetPixel(pixelIndex)!=0)
	      {
		addTotal = addTotal+readerMask->GetOutput()->GetPixel(pixelIndex);
		xExtend = xExtend + i;
		yExtend = yExtend + j;
	      }
	  }
      if(addTotal==0)
	zExtend--;
      else
	{
	  xExtend = int(float(xExtend)/float(addTotal));
	  yExtend = int(float(yExtend)/float(addTotal));
	  statZero = false;
	}
    }
  float addLength = (xExtend-branchList[0].voxelIndex[0].x)*(xExtend-branchList[0].voxelIndex[0].x);
  addLength = addLength + (yExtend-branchList[0].voxelIndex[0].y)*(yExtend-branchList[0].voxelIndex[0].y);
  addLength = addLength + (zExtend-branchList[0].voxelIndex[0].z)*(zExtend-branchList[0].voxelIndex[0].z);
  addLength = sqrt(addLength);
  branchList[branchCt].voxelCt = int(addLength); 
  branchList[branchCt].label = 0;
  branchList[branchCt].parLabel = 0;
  branchList[branchCt].voxelIndex = new VoxelType[branchList[branchCt].voxelCt*2];
  branchList[branchCt].voxelTangent = new VoxelType[branchList[branchCt].voxelCt*2];
  for(int i=0; i<branchList[branchCt].voxelCt; i++)
    {
      //Initialize at branchList[0].voxelIndex[0] 
      branchList[branchCt].voxelIndex[i].x = int(xExtend - float(i)*(xExtend-branchList[0].voxelIndex[0].x)/addLength); 
      branchList[branchCt].voxelIndex[i].y = int(yExtend - float(i)*(yExtend-branchList[0].voxelIndex[0].y)/addLength); 
      branchList[branchCt].voxelIndex[i].z = int(zExtend - float(i)*(zExtend-branchList[0].voxelIndex[0].z)/addLength); 
      branchList[branchCt].voxelTangent[i].x = 0;  
      branchList[branchCt].voxelTangent[i].y = 0; 
      branchList[branchCt].voxelTangent[i].z = 1; 
      branchList[branchCt].voxelIndex[i+branchList[branchCt].voxelCt].x = int(xExtend - float(i)*(xExtend-branchList[0].voxelIndex[0].x)/addLength); 
      branchList[branchCt].voxelIndex[i+branchList[branchCt].voxelCt].y = int(yExtend - float(i)*(yExtend-branchList[0].voxelIndex[0].y)/addLength); 
      branchList[branchCt].voxelIndex[i+branchList[branchCt].voxelCt].z = int(zExtend - float(i)*(zExtend-branchList[0].voxelIndex[0].z)/addLength); 
      branchList[branchCt].voxelTangent[i+branchList[branchCt].voxelCt].x = branchList[0].voxelTangent[0].x;  
      branchList[branchCt].voxelTangent[i+branchList[branchCt].voxelCt].y = branchList[0].voxelTangent[0].y; 
      branchList[branchCt].voxelTangent[i+branchList[branchCt].voxelCt].z = branchList[0].voxelTangent[0].z;
   }  
  branchList[branchCt].voxelCt = branchList[branchCt].voxelCt*2; 


  //Take sample lines, max line length 10 voxels
  //Sample lines start from outermost voxel of lumen
  // 0 0 0 0 0 0 0 0 9000 10000 11000 10500 10900 8000 7000 3000 0 0 0 0 0
  //                 |<-From here---------------------------To here->| 
  int sampleLengthMax = 100;
  for(a=0;a<branchCt+1;a++)
  //for(a=4;a<5;a++)
    {
      for(b=0;b<branchList[a].voxelCt;b++)
	//for(b=6;b<7;b++)
	{ 
	  index = branchList[a].voxelIndex[b];
	  tangent = branchList[a].voxelTangent[b]; 
	  //Compute Rotation system
	  float sin1, sin2, cos1, cos2;
	  if((tangent.y!=0)||(tangent.z!=0))
	    {
	      sin1 = tangent.y/sqrt(tangent.y*tangent.y+tangent.z*tangent.z);
	      cos1 = tangent.z/sqrt(tangent.y*tangent.y+tangent.z*tangent.z);
	      sin2 = tangent.x/sqrt(tangent.x*tangent.x+tangent.y*tangent.y+tangent.z*tangent.z);
	      cos2 = sqrt(tangent.y*tangent.y+tangent.z*tangent.z)/sqrt(tangent.x*tangent.x+tangent.y*tangent.y+tangent.z*tangent.z);
	    }
	  else
	    {
	      sin1 = 1;
	      cos1 = 0;
	      sin2 = 1;
	      cos2 = 0;
	    }
	  
	  //Plane orthogonal to tangent and through index
	  //In theta-length manner
	  int valueRec[angleSize+1][sampleSize];
	  int innerLoc[angleSize+1];
	  int widthLoc[angleSize+1];
	  int angleCt, sampleCt;
	  int theta;
	  for(angleCt=0;angleCt<angleSize;angleCt++)
	    {
	      innerLoc[angleCt]=0;
	      widthLoc[angleCt]=0;
	      for(sampleCt=0;sampleCt<sampleSize;sampleCt++)
		valueRec[angleCt][sampleCt]=0;
	    }


	  for(angleCt=0;angleCt<angleSize;angleCt++)
	  //for(angleCt=71;angleCt<angleSize;angleCt++)
	    {
	      theta = (float(angleCt))*360.0/angleSize;
	      bool flagInner = true;
	      bool flagTrack = true;
	      int innerCt = 1;
	      int sampleLineCt = 0;

	      for(sampleCt=1; (sampleCt<sampleLengthMax)&&flagTrack ; sampleCt++)
		{
		  //Get coordinate after rotation
		  float x = float(sampleCt)*cos(theta*PI/180);
		  float y = float(sampleCt)*sin(theta*PI/180);
		  float xG = cos2*x;
		  float yG = -sin1*sin2*x+cos1*y;
		  float zG = -cos1*sin2*x-sin1*y;
		  xG = xG + index.x;
		  yG = yG + index.y;
		  zG = zG + index.z;
		  if(xG<0) xG = 0;
		  if(xG>size[0]-1) xG = size[0]-1; 
		  if(yG<0) yG = 0;
		  if(yG>size[1]-1) yG = size[1]-1; 
		  if(zG<0) zG = 0;
		  if(zG>size[2]-1) zG = size[2]-1; 
		  
		  pixelIndex[0] = int(xG);
		  pixelIndex[1] = int(yG);
		  pixelIndex[2] = int(zG);

		  int voxelValue;
		  voxelValue = reader->GetOutput()->GetPixel(pixelIndex); 
		  //imageOut->SetPixel(pixelIndex,1);
		  //std::cout<<innerCt<<","<<sampleLineCt<<" ";		  

		  if(voxelValue!=0)
		    {
		      flagInner = false;
		    }
		  if(flagInner)
		    innerCt++;
		  if(!flagInner)
		    {
		      valueRec[angleCt][sampleLineCt] = voxelValue;
		      sampleLineCt++;		      
		      if((sampleLineCt>sampleSize)||(sampleLineCt>innerCt+3))
			flagTrack = false;
		    }
		}
	      innerLoc[angleCt] = innerCt;

	      //Find local max
	      i = 0;
	      bool flagLoc = true;
	      while(flagLoc)
		{
		  i++;
		  if(((valueRec[angleCt][i]>=valueRec[angleCt][i-1])&&(valueRec[angleCt][i]>=valueRec[angleCt][i+1]))||(i>8)||((valueRec[angleCt][i]-valueRec[angleCt][i+1])>2000))
		    flagLoc = false;
		}
	      if(i>8) i=8;
	      if(valueRec[angleCt][i+1]!=0)
		widthLoc[angleCt] = i;
	      else
		widthLoc[angleCt] = 1;

	      if(widthLoc[i]<1)
		widthLoc[i]=1;
	    }

	  //Angular sampling finished
	  //Available:
	  //innerLoc
	  //widthLoc
	  //valueRec
	  //Compute histogram and select the first three largest group
	  int hist[10];
	  for(i=0;i<10;i++)
	    hist[i]=0;
	  for(i=0;i<angleSize;i++)
	    hist[widthLoc[i]]++;
	  std::vector<int> histVec;
	  for(i=0;i<10;i++)
	    {
	      histVec.push_back(hist[i]);
	      //std::cout<<hist[i]<<" ";
	    }
	  std::sort(histVec.begin(),histVec.end());
	  int first = histVec[histVec.size()-1];
	  int second = histVec[histVec.size()-2];
	  int third = histVec[histVec.size()-3];
	  if((second==0)||((first-second)>15)||(second<10))
	    second = first;	  
	  if((third==0)||((second-third)>10)||(third<10))
	    third = second;


	  int locThre;
	  for(i=0;i<angleSize;i++)
	    if((hist[i]==first)||(hist[i]==second)||(hist[i]==third))
	      locThre = i;
	  //std::cout<<": "<<locThre<<"; "<<b<<std::endl;
	  if(locThre>3) locThre = 3;

	  float xG, yG, zG;
	  for(i=0;i<angleSize;i++)
	    //for(i=71;i<angleSize;i++)
	    {
	      theta = i*360/angleSize;

	      if(widthLoc[i]>locThre)
		widthLoc[i]=locThre;

	      if((innerLoc[i]>10)&&(widthLoc[i]<2))
		widthLoc[i]=2;
	      /*
	      if(innerLoc[i]>10)
		{
		  bool hitflag=true;
		  while(hitflag)
		    {
		      float x = float(innerLoc[i]+widthLoc[i])*cos(theta*PI/180);
		      float y = float(innerLoc[i]+widthLoc[i])*sin(theta*PI/180);
		      xG = cos2*x;
		      yG = -sin1*sin2*x+cos1*y;
		      zG = -cos1*sin2*x-sin1*y;
		      xG = xG + index.x;
		      yG = yG + index.y;
		      zG = zG + index.z;
		      if(xG<0) xG = 0;
		      if(xG>size[0]-1) xG = size[0]-1; 
		      if(yG<0) yG = 0;
		      if(yG>size[1]-1) yG = size[1]-1; 
		      if(zG<0) zG = 0;
		      if(zG>size[2]-1) zG = size[2]-1; 
		  
		      pixelIndex[0] = int(xG);
		      pixelIndex[1] = int(yG);
		      pixelIndex[2] = int(zG);

		      int voxelValue = reader->GetOutput()->GetPixel(pixelIndex); 
		      if(voxelValue<10000)
			widthLoc[i]++;
		      else
			hitflag = false;
		    }
		}
	      */
	      for(j=0;j<=widthLoc[i];j++)
		{
		  float x = float(innerLoc[i]+j)*cos(theta*PI/180);
		  float y = float(innerLoc[i]+j)*sin(theta*PI/180);
		  xG = cos2*x;
		  yG = -sin1*sin2*x+cos1*y;
		  zG = -cos1*sin2*x-sin1*y;
		  xG = xG + index.x;
		  yG = yG + index.y;
		  zG = zG + index.z;
		  if(xG<0) xG = 0;
		  if(xG>size[0]-1) xG = size[0]-1; 
		  if(yG<0) yG = 0;
		  if(yG>size[1]-1) yG = size[1]-1; 
		  if(zG<0) zG = 0;
		  if(zG>size[2]-1) zG = size[2]-1; 
		  
		  pixelIndex[0] = int(xG);
		  pixelIndex[1] = int(yG);
		  pixelIndex[2] = int(zG);

		  int voxelValue = reader->GetOutput()->GetPixel(pixelIndex); 
		  if(voxelValue!=0)
		    imageOutInner->SetPixel(pixelIndex,1);		    
		} 
	      /*
	      std::cout<<i<<": ";
	      for(j=0;j<sampleSize;j++)
		std::cout<<valueRec[i][j]<<" ";
	      std::cout<<"; "<<innerLoc[i]<<" "<<widthLoc[i]<<" "<<reader->GetOutput()->GetPixel(pixelIndex)<<std::endl;
	      */

	      float x = float(innerLoc[i]+widthLoc[i] + 3)*cos(theta*PI/180);
	      float y = float(innerLoc[i]+widthLoc[i] + 3)*sin(theta*PI/180);
	      xG = cos2*x;
	      yG = -sin1*sin2*x+cos1*y;
	      zG = -cos1*sin2*x-sin1*y;
	      xG = xG + index.x;
	      yG = yG + index.y;
	      zG = zG + index.z;
	      if(xG<0) xG = 0;
	      if(xG>size[0]-1) xG = size[0]-1; 
	      if(yG<0) yG = 0;
	      if(yG>size[1]-1) yG = size[1]-1; 
	      if(zG<0) zG = 0;
	      if(zG>size[2]-1) zG = size[2]-1; 
		  
	      pixelIndex[0] = int(xG);
	      pixelIndex[1] = int(yG);
	      pixelIndex[2] = int(zG);
	      int seedValue = imageOutInner->GetPixel(pixelIndex);	
	      int voxelValue = reader->GetOutput()->GetPixel(pixelIndex); 
	      if((seedValue!=1)&&(voxelValue!=0))     
		imageOutOuter->SetPixel(pixelIndex,1);
	    }
	}
    }

  typedef itk::MedianImageFilter< ImageType, ImageType > MedianFilterType;
  MedianFilterType::Pointer medianFilterInner = MedianFilterType::New();
  MedianFilterType::Pointer medianFilterOuter = MedianFilterType::New();
  ImageType::SizeType indexRadius;
  indexRadius[0] = 0; 
  indexRadius[1] = 0; 
  indexRadius[2] = 0; 
  medianFilterInner->SetRadius( indexRadius );
  medianFilterOuter->SetRadius( indexRadius );
  medianFilterInner->SetInput( imageOutInner );
  medianFilterOuter->SetInput( imageOutOuter );
  medianFilterInner->Update();
  medianFilterOuter->Update();




  //Write output  
  writerInner->SetInput( medianFilterInner->GetOutput() );
  writerInner->SetFileName( argv[3] );

  try
    {
      writerInner->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }
  
  return EXIT_SUCCESS;
}
