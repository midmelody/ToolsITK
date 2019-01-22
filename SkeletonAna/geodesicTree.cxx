#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkPoint.h"
#include <iostream>

//ITK settings
const unsigned int Dimension = 3;
typedef unsigned char PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;
const unsigned int Dimension2D = 2;
typedef itk::Image< PixelType, Dimension2D > ImageType2D;
typedef itk::ConnectedThresholdImageFilter< ImageType2D, ImageType2D > ConnectedThresholdFilterType;
typedef itk::BinaryContourImageFilter < ImageType2D, ImageType2D > BinaryContourImageFilterType;
typedef itk::ImageFileWriter< ImageType2D > WriterType2D;
typedef unsigned short LabelType;
typedef itk::ShapeLabelObject< LabelType, Dimension2D > ShapeLabelObjectType;
typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;
typedef itk::LabelImageToShapeLabelMapFilter< ImageType2D, LabelMapType> ImageToLableFilterType;
//Voxel struct
struct VoxelType
{
  float x;
  float y;
  float z;
  float tanX;
  float tanY;
  float tanZ;
  int label;
  int parLabel; 
  int genLabel;
  float distToRoot;
  float lumenRadius;
  float wallThickness;
  float roundnessLumen;
  float roundnessWall;
};

//Voxel Simple struct
struct VoxelSimpleType
{
  int x;
  int y;
  int z;
};

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputTreeTangentFile inputTreeGenFile inputAirwayWallImageFile outputFullTreeInfoFile [outputGeodesicTreeImageFile]" << std::endl;
      return EXIT_FAILURE;
    }

  //Counter and other voxel values
  int i, j, k;
  float x, y, z;
  float tanX, tanY, tanZ;
  int label, parLabel, genLabel;
  float distToRoot;
 
  //To record: accumulated length for each branch segment, parent branch label, last voxel number for each branch 
  float lenRec[1000];
  int parRec[1000];
  int lastVoxNum[1000];
  bool lenRecUpdateFlag[1000];
  for(i=0; i<1000; i++)
    {
      lenRec[i] = 0.0;
      parRec[i] = 0;  
      lenRecUpdateFlag[i] = false;
      lastVoxNum[i] = 0;
    }
  lenRecUpdateFlag[0] = true;
  lenRecUpdateFlag[1] = true;
  //Record total voxel ct
  VoxelType voxel;
  int totalVoxelCt = 0;
  FILE *treeTanFile = fopen(argv[1], "r");
  while(fscanf(treeTanFile, "%f %f %f %d %d %f %f %f", &x, &y, &z, &label, &parLabel, &tanX, &tanY, &tanZ) == 8)
    totalVoxelCt++;
  fclose(treeTanFile);
 
  //Import all data
  VoxelType* voxelList = new VoxelType[totalVoxelCt]; 
  treeTanFile = fopen(argv[1], "r");
  FILE *treeGenFile = fopen(argv[2], "r");
  i = 0;
  while(fscanf(treeTanFile, "%f %f %f %d %d %f %f %f", &x, &y, &z, &label, &parLabel, &tanX, &tanY, &tanZ) == 8)
    {
      fscanf(treeGenFile, "%f %f %f %d %d", &x, &y, &z, &label, &genLabel);  
      voxelList[i].x = x;
      voxelList[i].y = y;
      voxelList[i].z = z;
      voxelList[i].tanX = tanX;
      voxelList[i].tanY = tanY;
      voxelList[i].tanZ = tanZ;
      voxelList[i].label = label;
      voxelList[i].parLabel = parLabel;
      voxelList[i].genLabel = genLabel;
      voxelList[i].distToRoot = 0;
      voxelList[i].roundnessLumen = 0;
      voxelList[i].roundnessWall = 0;
      parRec[label] = parLabel;
      i++;
    }    
  fclose(treeTanFile);
  fclose(treeGenFile);

  //Read airway wall image 
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[3] );
  reader->Update();
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
  ImageType::Pointer geodesicImage = ImageType::New();
  geodesicImage->SetRegions( region );
  geodesicImage->SetSpacing( spacing );
  geodesicImage->SetOrigin( origin );
  geodesicImage->SetDirection( direction );
  geodesicImage->Allocate();  
  ImageType::IndexType pixelIndex;

  //Assign geodesic distance record for all voxels in the list
  ImageType::IndexType preIdx, curIdx;
  ImageType::PointType prePhy, curPhy;
  preIdx[0] = voxelList[0].x;
  preIdx[1] = voxelList[0].y;
  preIdx[2] = voxelList[0].z;
  curIdx[0] = voxelList[0].x;
  curIdx[1] = voxelList[0].y;
  curIdx[2] = voxelList[0].z;
  geodesicImage->TransformIndexToPhysicalPoint(preIdx,prePhy);
  geodesicImage->TransformIndexToPhysicalPoint(curIdx,curPhy);
  int preVoxelLab = 0;
  int preVoxelNum;
  float tempLength;
  float maxDist = 0;

  //2D local ROI analysis
  int sizeROI = 30;
  ImageType2D::SizeType size2D;
  size2D[0] = 2*sizeROI+1;
  size2D[1] = 2*sizeROI+1;
  ImageType2D::RegionType region2D;
  region2D.SetSize( size2D );
  VoxelSimpleType** recIdx = new VoxelSimpleType*[size2D[0]];
  for(j=0; j<size2D[0]; j++)
    recIdx[j] = new VoxelSimpleType[size2D[1]];
  ImageType::IndexType indexCenter, indexCurrent;
  ImageType::PointType pointCenter, pointCurrent;
  int lumenFlag, wallFlag;
  float lumenRadius, wallThickness;
  float lumenCount, wallCount;
  float roundnessLumen, roundnessWall;

  //Examine every centerline points for distToRoot, lumenRadius, and wallThickness
  for(i=0;i<totalVoxelCt;i++)
    {  
      //Set current point
      curIdx[0] = voxelList[i].x;
      curIdx[1] = voxelList[i].y;
      curIdx[2] = voxelList[i].z;
      //Use the last voxel of parent segment as previous voxel
      //at the beginning of each new branch segment (label)
      if(voxelList[i].label!=preVoxelLab)
	{
	  lastVoxNum[preVoxelLab] = i-1;
	  //Set new preIdx from the second branch segment
	  if(voxelList[i].label>1)
	    preVoxelNum = lastVoxNum[parRec[voxelList[i].label]];
	  else
	    preVoxelNum = 0;
	  //Set new preVoxelLab
	  preVoxelLab = voxelList[i].label;
	}
      else
	preVoxelNum = i-1;
      preIdx[0] = voxelList[preVoxelNum].x;
      preIdx[1] = voxelList[preVoxelNum].y;
      preIdx[2] = voxelList[preVoxelNum].z;

      //Compute current distance in physical space
      geodesicImage->TransformIndexToPhysicalPoint(preIdx,prePhy);
      geodesicImage->TransformIndexToPhysicalPoint(curIdx,curPhy);
      tempLength = curPhy.EuclideanDistanceTo(prePhy);
      lenRec[voxelList[i].label] = lenRec[voxelList[i].label] + tempLength;
      voxelList[i].distToRoot = voxelList[preVoxelNum].distToRoot + tempLength;  
      if(maxDist < voxelList[i].distToRoot)
	maxDist = voxelList[i].distToRoot;  

      //Resample 2D local plan according to tangent vector
      int tempX, tempY, tempZ;
      int cenX, cenY, cenZ;
      float graX, graY, graZ;
      cenX = voxelList[i].x;
      cenY = voxelList[i].y;
      cenZ = voxelList[i].z;
      graX = voxelList[i].tanX;
      graY = voxelList[i].tanY;
      graZ = voxelList[i].tanZ;
      //Compute Rotation system
      float sin1, sin2, cos1, cos2;
      if((graY!=0)||(graZ!=0))
	{
	  sin1 = graY/sqrt(graY*graY+graZ*graZ);
	  cos1 = graZ/sqrt(graY*graY+graZ*graZ);
	  sin2 = graX/sqrt(graX*graX+graY*graY+graZ*graZ);
	  cos2 = sqrt(graY*graY+graZ*graZ)/sqrt(graX*graX+graY*graY+graZ*graZ);
	}
      else
	{
	  sin1 = 1;
	  cos1 = 0;
	  sin2 = 1;
	  cos2 = 0;
	}
      
      //Generate 2D orthogonal image and record 3D-2D index relation
      ImageType2D::Pointer crossImage = ImageType2D::New();
      crossImage->SetRegions( region2D );
      crossImage->Allocate();  
      ImageType2D::IndexType pixelIndex2D;
      PixelType pixelValue;
      for(tempX=-sizeROI; tempX<=sizeROI; tempX++)
	for(tempY=-sizeROI; tempY<=sizeROI; tempY++)
	  {
	    //Get coordinate after rotation
	    float xG = cos2*float(tempX);
	    float yG = -sin1*sin2*float(tempX)+cos1*float(tempY);
	    float zG = -cos1*sin2*float(tempX)-sin1*float(tempY);
	    xG = xG + cenX;
	    yG = yG + cenY;
	    zG = zG + cenZ;
	    pixelIndex[0] = int(xG);
	    pixelIndex[1] = int(yG);
	    pixelIndex[2] = int(zG);
	    if((pixelIndex[0]>=0)&&(pixelIndex[0]<pRow)&&(pixelIndex[1]>=0)&&(pixelIndex[1]<pCol)&&(pixelIndex[2]>=0)&&(pixelIndex[2]<pSli))
	      {  
		pixelValue = reader->GetOutput()->GetPixel(pixelIndex);
		pixelIndex2D[0] = tempX+sizeROI;
		pixelIndex2D[1] = tempY+sizeROI;	    
		crossImage->SetPixel(pixelIndex2D,pixelValue*100);
		
		recIdx[pixelIndex2D[0]][pixelIndex2D[1]].x = pixelIndex[0];
		recIdx[pixelIndex2D[0]][pixelIndex2D[1]].y = pixelIndex[1];
		recIdx[pixelIndex2D[0]][pixelIndex2D[1]].z = pixelIndex[2];
	      }
	  }
      
      //Extract lumen and wall
      ImageType2D::IndexType indexCen;
      indexCen[0] = sizeROI;
      indexCen[1] = sizeROI;
      ConnectedThresholdFilterType::Pointer connectedThreshold = ConnectedThresholdFilterType::New();
      connectedThreshold->SetInput( crossImage );
      connectedThreshold->SetReplaceValue( 255 );
      connectedThreshold->SetSeed( indexCen );
      connectedThreshold->SetLower( 199 );
      connectedThreshold->SetUpper( 201 );
      connectedThreshold->Update();

      //Use lumen to compute roundness measurement
      ImageType2D::Pointer lumenImage = ImageType2D::New();
      lumenImage = connectedThreshold->GetOutput();
      //WriterType2D::Pointer writer2DT = WriterType2D::New();
      //writer2DT->SetFileName( "../Experiment/Temp/lumen.png" );
      //writer2DT->SetInput( lumenImage );
      //writer2DT->SetInput(crossImage);
      //writer2DT->Update();
      ImageToLableFilterType::Pointer imageToLabel = ImageToLableFilterType::New();
      imageToLabel->SetInput( lumenImage );
      LabelMapType::Pointer labelMap = LabelMapType::New();
      labelMap = imageToLabel->GetOutput();
      imageToLabel->Update();
      int numObj = labelMap->GetNumberOfLabelObjects();
      if(numObj == 0)
	{
	  lumenRadius = 0;
	  wallThickness = 0;
	  roundnessLumen = 0;
	  roundnessWall = 0;
	}
      else
	{
	  //Get local lumen roundness
	  ShapeLabelObjectType::Pointer labelObject = ShapeLabelObjectType::New();
	  labelObject = labelMap->GetNthLabelObject(0);
	  roundnessLumen = labelObject->GetRoundness();
	  lumenImage->DisconnectPipeline();
	  //Get local wall roundness	  
	  connectedThreshold->SetLower( 99 );
	  connectedThreshold->SetUpper( 201 );
	  ImageType2D::Pointer wallImage = ImageType2D::New();
	  wallImage = connectedThreshold->GetOutput();
	  wallImage->Update();
	  imageToLabel->SetInput( wallImage );
	  imageToLabel->Update();
	  labelMap = imageToLabel->GetOutput();
	  labelObject = labelMap->GetNthLabelObject(0);
	  roundnessWall = labelObject->GetRoundness();
	  wallImage->DisconnectPipeline();
	  //Get contours
	  BinaryContourImageFilterType::Pointer binaryContour = BinaryContourImageFilterType::New();
	  binaryContour->SetInput( connectedThreshold->GetOutput() );
	  //Get lumen contour
	  connectedThreshold->SetLower( 199 );
	  connectedThreshold->SetUpper( 201 );
	  ImageType2D::Pointer lumenContour = ImageType2D::New();
	  lumenContour = binaryContour->GetOutput();
	  binaryContour->Update();
	  lumenContour->DisconnectPipeline();
	  //Get wall contour
	  connectedThreshold->SetLower( 99 );
	  connectedThreshold->SetUpper( 201 );
	  ImageType2D::Pointer wallContour = ImageType2D::New();
	  wallContour = binaryContour->GetOutput();
	  binaryContour->Update();
	  wallContour->DisconnectPipeline();
      
	  //Write temp output
	  if(i==180)
	    {
	      WriterType2D::Pointer writer2D = WriterType2D::New();
	      //writer2D->SetFileName( "../Experiment/Temp/lumenContour.png" );
	      //writer2D->SetInput( lumenContour );
	      //writer2D->Update();
	      //writer2D->SetFileName( "../Experiment/Temp/wallContour.png" );
	      //writer2D->SetInput( wallContour );
	      //writer2D->Update();
	      writer2D->SetFileName( "../Experiment/Temp/crossImage.png" );
	      writer2D->SetInput( crossImage );
	      writer2D->Update();
	    }
      
	  //Examine the lumen radius and wall thickness
	  indexCenter[0] = recIdx[sizeROI][sizeROI].x;
	  indexCenter[1] = recIdx[sizeROI][sizeROI].y;	  
	  indexCenter[2] = recIdx[sizeROI][sizeROI].z;	
	  geodesicImage->TransformIndexToPhysicalPoint(indexCenter,pointCenter);
      
	  lumenRadius = 0;
	  wallThickness = 0;
	  lumenCount = 0;
	  wallCount = 0;
	  for(tempX=-sizeROI; tempX<=sizeROI; tempX++)
	    for(tempY=-sizeROI; tempY<=sizeROI; tempY++)
	      {
		pixelIndex2D[0] = tempX+sizeROI;
		pixelIndex2D[1] = tempY+sizeROI;		
		lumenFlag = lumenContour->GetPixel(pixelIndex2D);
		wallFlag = wallContour->GetPixel(pixelIndex2D);
		
		indexCurrent[0] =  recIdx[pixelIndex2D[0]][pixelIndex2D[1]].x;
		indexCurrent[1] =  recIdx[pixelIndex2D[0]][pixelIndex2D[1]].y;
		indexCurrent[2] =  recIdx[pixelIndex2D[0]][pixelIndex2D[1]].z;
	    
		if(lumenFlag == 255)
		  {
		    geodesicImage->TransformIndexToPhysicalPoint(indexCurrent,pointCurrent);
		    tempLength = pointCenter.EuclideanDistanceTo(pointCurrent);
		    lumenRadius = lumenRadius+tempLength;
		    lumenCount++;
		  }
		
		if(wallFlag == 255)
		  {
		    geodesicImage->TransformIndexToPhysicalPoint(indexCurrent,pointCurrent);
		    tempLength = pointCenter.EuclideanDistanceTo(pointCurrent);
		    wallThickness = wallThickness+tempLength;
		    wallCount++;
		  }
	      }
	  
	  lumenRadius = lumenRadius/lumenCount;
	  wallThickness = wallThickness/wallCount - lumenRadius;
	}
      voxelList[i].lumenRadius = lumenRadius;
      voxelList[i].wallThickness = wallThickness;
      voxelList[i].roundnessLumen = roundnessLumen;
      voxelList[i].roundnessWall = roundnessWall;
    }
  
  //Calculate "base length" for every branch segment
  for(i=1;i<1000;i++)
    {
      //if((lenRec[i]>0)&&(parRec[i]>=0))
      //std::cout<<i<<": "<<lenRec[i]<<" "<<parRec[i]<<" ";
      if(lenRecUpdateFlag[parRec[i]])
	{
	  if(parRec[i]>0)
	    {
	      lenRec[i] = lenRec[i] + lenRec[parRec[i]];
	      lenRecUpdateFlag[i] = true;
	    }
	  //if((lenRec[i]>0)&&(parRec[i]>=0))
	  //{
	  //  std::cout<<lenRec[i]<<" ";
	  //  std::cout<<voxelList[lastVoxNum[i]].x<<" ";
	  //  std::cout<<voxelList[lastVoxNum[i]].y<<" ";	      
	  //  std::cout<<voxelList[lastVoxNum[i]].z<<std::endl;
	  //}
	}
      else
	{
	  std::cerr << "Parent branch length has not been updated!" << std::endl;
	  return EXIT_FAILURE;
	}
    }


  //Output
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  geodesicImage->SetPixel(pixelIndex, 0);
	}

  FILE *treeTotalFile = fopen(argv[4], "w");
  fprintf(treeTotalFile, "%s %s %s %s %s %s %s %s %s %s %s\n", "x", "y", "z", "label", "parent", "generation", "distToRoot", "lumenRadius", "wallThickness", "roundnessLumen", "roundnessWall");
  for(i=0;i<totalVoxelCt;i++)
    {
      x = voxelList[i].x;
      y = voxelList[i].y;
      z = voxelList[i].z;
      label = voxelList[i].label;
      parLabel = voxelList[i].parLabel;
      genLabel = voxelList[i].genLabel;
      distToRoot = voxelList[i].distToRoot;
      lumenRadius = voxelList[i].lumenRadius;
      wallThickness = voxelList[i].wallThickness;
      roundnessLumen = voxelList[i].roundnessLumen;
      roundnessWall = voxelList[i].roundnessWall;
      fprintf(treeTotalFile, "%d %d %d %d %d %d %f %f %f %f %f\n", int(x), int(y), int(z), label, parLabel, genLabel, distToRoot, lumenRadius, wallThickness, roundnessLumen, roundnessWall);

      pixelIndex[0]=int(x);
      pixelIndex[1]=int(y);
      pixelIndex[2]=int(z);
      geodesicImage->SetPixel(pixelIndex, round(distToRoot/maxDist*255));
    }
  fclose(treeTotalFile);

  if(argc == 6)
    {
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( argv[5] );
      writer->SetInput( geodesicImage );
      writer->Update();
    }

  for(j=0; j<size2D[0]; j++) 
    {
      delete [] recIdx[j];
    }
  delete [] recIdx;
  delete [] voxelList;
  return EXIT_SUCCESS;
}
