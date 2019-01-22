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
#include <queue>

#define input(i,j,k) input[(k)*pRow*pCol+(j)*pRow+(i)]
#define label(i,j,k) label[(k)*pRow*pCol+(j)*pRow+(i)]
#define parentLabel(i,j,k) parentLabel[(k)*pRow*pCol+(j)*pRow+(i)]

//ITK settings
const unsigned int Dimension = 3;
typedef unsigned char PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

//struct for point index 
struct TreePointType
{
  int x;
  int y;
  int z;
  int label;
  int parentLabel;
};
typedef std::queue< TreePointType > rootListType;

//Trace one branch
void traceBranch(TreePointType root, rootListType &rootList, int *input, int pRow, int pCol, int pSli, int & curLabel)
{
  FILE *treeFile = fopen("tree.txt", "a+"); 
  //std::cout<<"*********************************************"<<std::endl;
  //std::cout<<root.x+1<<" "<<root.y+1<<" "<<root.z+1<<" "<<curLabel<<" "<<root.parentLabel<<std::endl;
  fprintf(treeFile, "%d %d %d %d %d\n", root.x, root.y, root.z, curLabel, root.parentLabel) ;
  int parLabel = root.parentLabel;
  int x = root.x;
  int y = root.y;
  int z = root.z;

  input(x,y,z) = 0;
  
  //trace branch
  //Mark = 1: end point;
  //Mark = 2: line point;
  //Mark > 2: branching point;
  //Visited point set "input" to 0
  int a,b,c;
  int xNei, yNei, zNei;
  int xNext, yNext, zNext;
  int xNN, yNN, zNN;
  int curMark = 2;
  while(curMark == 2)
    {
      //std::cout<<"-----------------------------"<<std::endl;
      //Examine neighborhood
      int nonZeroCt = 0;
      for(a=-1;a<=1;a++)
	for(b=-1;b<=1;b++)
	  for(c=-1;c<=1;c++)
	    {
	      xNei = x+a;
	      yNei = y+b;
	      zNei = z+c;
	      if((xNei>=1)&&(xNei<pRow)&&(yNei>=0)&&(yNei<pCol)&&(zNei>=0)&&(zNei<pSli))
		{	     
		  if(input(xNei,yNei,zNei))
		    {
		      //std::cout<<"Current "<<xNei+1<<" "<<yNei+1<<" "<<zNei+1<<" "<<curLabel<<" "<<parLabel<<": "<<input(xNei,yNei,zNei)<<std::endl;
		      fprintf(treeFile, "%d %d %d %d %d\n", xNei, yNei, zNei, curLabel, parLabel) ;
		      //end point, mark and stop
		      if(input(xNei,yNei,zNei) == 1)
			{
			  input(xNei,yNei,zNei) = 0;
			  curMark = 1;
			  nonZeroCt++;
			}
		      //line point, mark and continue tracking
		      if(input(xNei,yNei,zNei) == 2)
			{
			  input(xNei,yNei,zNei) = 0;
			  nonZeroCt++;
			  xNext = xNei;
			  yNext = yNei;
			  zNext = zNei;
			}
		      //branching point, mark, do a local search further until hit all 2s then put 2s to root list and stop
		      if(input(xNei,yNei,zNei) > 2)
			{
			  input(xNei,yNei,zNei) = 0;	      
			  curMark = input(xNei,yNei,zNei);
			  nonZeroCt++;
			  //put to temp root list
			  int curParLabel = parLabel;
			  parLabel = curLabel;
			  int curCand = curMark;
			  rootListType tempRootList;	  
			  TreePointType tempRootPoint;
			  tempRootPoint.x = xNei;
			  tempRootPoint.y = yNei;
			  tempRootPoint.z = zNei;
			  tempRootPoint.label = 0;
			  tempRootPoint.parentLabel = 0;
			  tempRootList.push(tempRootPoint);
			  //search neighborhood, until find all 2s in the local neighborhood
			  while(!tempRootList.empty())
			    {
			      TreePointType candRootPoint;
			      int xTemp, yTemp, zTemp;
			      candRootPoint = tempRootList.front();
			      tempRootList.pop();
			      xTemp = candRootPoint.x;
			      yTemp = candRootPoint.y;
			      zTemp = candRootPoint.z;
			      for(int i=-1;i<=1;i++)
				for(int j=-1;j<=1;j++)
				  for(int k=-1;k<=1;k++)
				    {    
				      xNN = xTemp+i;
				      yNN = yTemp+j;
				      zNN = zTemp+k;
				      if((xNN>=1)&&(xNN<pRow)&&(yNN>=0)&&(yNN<pCol)&&(zNN>=0)&&(zNN<pSli))
					{
					  if(input(xNN,yNN,zNN)==1)
					    {
					      input(xNN,yNN,zNN) = 0;			      
					    }
					  if(input(xNN,yNN,zNN)==2)
					    {
					      //record root
					      TreePointType rootPoint;
					      rootPoint.x = xNN;
					      rootPoint.y = yNN;
					      rootPoint.z = zNN;
					      rootPoint.label = curLabel;
					      rootPoint.parentLabel = parLabel;
					      rootList.push(rootPoint);  
					      input(xNN,yNN,zNN) = 0;
					      //std::cout<<rootPoint.x+1<<" "<<rootPoint.y+1<<" "<<rootPoint.z+1<<std::endl;
					    }				
					  else if(input(xNN,yNN,zNN)!=0)
					    {
					      //continue search
					      fprintf(treeFile, "%d %d %d %d %d\n", xNN, yNN, zNN, curLabel, curParLabel) ;
					      tempRootPoint.x = xNN;
					      tempRootPoint.y = yNN;
					      tempRootPoint.z = zNN;
					      tempRootPoint.label = 0;
					      tempRootPoint.parentLabel = 0;
					      tempRootList.push(tempRootPoint);  
					      input(xNN,yNN,zNN) = 0;
					    }	    
					}
				    }			      
			    }
			}
		    }
		}
	    }
      //if isolated root, stop
      if(nonZeroCt == 0)
	curMark = 1;
      //renew current candidate point
      x = xNext;
      y = yNext;
      z = zNext;

    }

  fclose(treeFile);
}

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile treacheaImageFile sliceOrientation (0-low2high; 1-high2low)" << std::endl;
      return EXIT_FAILURE;
    }
 
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerT = ReaderType::New();
  //Parameters
  reader->SetFileName( argv[1] );
  readerT->SetFileName( argv[2] );
 
  //Pipeline
  try
    {
      reader->Update();
      readerT->Update();
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
 
  //temp space for computation
  ImageType::IndexType pixelIndex;
  int i, j, k;
  int *input;
  input = (int *)malloc(sizeof( int ) * pRow * pCol * pSli);
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  input(i,j,k) = reader->GetOutput()->GetPixel(pixelIndex);
	}
  

  //Record and find the uppermost end point;
  int uppX,uppY,uppZ;
  bool foundFlag = false;
  bool dirAir;
  bool traLoc;
  
  if(atoi(argv[3])==0)
    dirAir = true;
  else
    dirAir = false;

  if(dirAir)
    k = 0;
  else
    k = pSli-1;

  //std::cout<<dirAir<<" "<<k<<std::endl;
  
  while(!foundFlag)
    {
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {
	    pixelIndex[0]=i;
	    pixelIndex[1]=j;
	    pixelIndex[2]=k;
	    traLoc  = readerT->GetOutput()->GetPixel(pixelIndex);
	    if((input(i,j,k)==1)&&(traLoc==1))
	      {
		uppX = i;
		uppY = j;
		uppZ = k;
		foundFlag = true;
	      }
	  }
      if(dirAir)
	k = k+1;
      else
	k = k-1;
    }
	
  std::cout<<"Start: "<<uppX<<" "<<uppY<<" "<<uppZ<<std::endl;

  TreePointType rootPoint;
  rootListType rootList;
  //Initialize
  rootPoint.x = uppX;
  rootPoint.y = uppY;
  rootPoint.z = uppZ;
  rootPoint.label = 1;
  rootPoint.parentLabel = 0;
  rootList.push(rootPoint);
  int curLabel = 1;

  //tree construction

  while(!rootList.empty())
    {
      rootPoint = rootList.front();
      rootList.pop();
      traceBranch(rootPoint, rootList, input, pRow, pCol, pSli, curLabel);
      //increse current label
      curLabel++;
    }
 
  return EXIT_SUCCESS;
}
