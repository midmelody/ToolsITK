void RasterScanEDT(int numberOfZeroCrossing, VOXEL *zeroCrossing, float *DT, VOXEL *NEP, int neighborCount) 
{
  int t;
  int i,j,k,a,b,c;
  float edgeX, edgeY, edgeZ;
  float distX, distY, distZ;
  float distOne;
  //neighborCount n;
  //Total member: n*n*n-1
  //Member for each mask: (n*n*n-1)/2
  int maskSize = (pow(neighborCount,3)-1)/2;
  int oneSideSize = (neighborCount-1)/2;
  VOXEL *forMask;
  VOXEL *backMask;
  forMask = (VOXEL *)malloc(sizeof(VOXEL)*maskSize); 
  backMask = (VOXEL *)malloc(sizeof(VOXEL)*maskSize); 

  //compute mask
  t = 0;
  for(i=-oneSideSize;i<=oneSideSize;i++)
    for(j=-oneSideSize;j<=oneSideSize;j++)
      for(k=-oneSideSize;k<=oneSideSize;k++)
	if((i!=0)||(j!=0)||(k!=0))
	  {
	    if(t<maskSize)
	      {
		forMask[t].x = i;
		forMask[t].y = j;
		forMask[t].z = k;
	      }
	    else
	      {
		backMask[t-maskSize].x = i;
		backMask[t-maskSize].y = j;
		backMask[t-maskSize].z = k;
	      }
	    t++;
	  }

  //Initialize distance transform and nearest edge point
  for(t=0;t<numberOfZeroCrossing;t++)
    {
      edgeX = zeroCrossing[t].x;
      edgeY = zeroCrossing[t].y;
      edgeZ = zeroCrossing[t].z;
      i = floor(edgeX);
      j = floor(edgeY);
      k = floor(edgeZ);
      distX = edgeX - i;
      distY = edgeY - j;
      distZ = edgeZ - k;
      for(a=0;a<=1;a++)
	for(b=0;b<=1;b++)
	  for(c=0;c<=1;c++)
	    {
	      if(((i+a)>=0)&&((i+a)<pRow)&&((j+b)>=0)&&((j+b)<pCol)&&((k+c)>=0)&&((k+c)<pSli))
		{
		  if(DT(i+a,j+b,k+c)>pow((a-distX),2)+pow((b-distY),2)+pow((c-distZ),2))
		    {
		      DT(i+a,j+b,k+c)=pow((a-distX),2)+pow((b-distY),2)+pow((c-distZ),2);
		      NEP(i+a,j+b,k+c).x = edgeX;
		      NEP(i+a,j+b,k+c).y = edgeY;
		      NEP(i+a,j+b,k+c).z = edgeZ;
		    }
		}
	    }
    }
  //Output zeroCrossing
  
  ImageType::SpacingType spacing;
  spacing[0] = spaRow;
  spacing[1] = spaCol;
  spacing[2] = spaSli;
  ImageType::SizeType  size;
  size[0] = pRow;
  size[1] = pCol;
  size[2] = pSli;
  ImageType::RegionType region;
  ImageType::IndexType start = {0, 0, 0};
  region.SetSize( size );
  region.SetIndex( start );
  ImageType::IndexType pixelIndex;
  ImageType::Pointer ZeroCrossingImage = ImageType::New();
  ZeroCrossingImage->SetRegions( region );
  ZeroCrossingImage->SetSpacing( spacing );
  ZeroCrossingImage->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  ZeroCrossingImage->SetPixel(pixelIndex,DT(i,j,k));
	}
 
  
  ImageWriterType::Pointer writerZeroCrossing = ImageWriterType::New();
  writerZeroCrossing->SetInput( ZeroCrossingImage );
  writerZeroCrossing->SetFileName("ZeroCrossing.img");
  writerZeroCrossing->Update();
  
  

  //Raster scan for DT 
  bool flagChange = true;
  float dist;
  int u,v,w;
  int countTimes = 0;
  while(flagChange)
    { 
      countTimes++;
      std::cout<<"Scan: "<<countTimes<<std::endl;
      flagChange = false;     
      //Forward
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  for(k=0;k<pSli;k++)
	    if(NEP(i,j,k).x!=-1)
	      {
		for(t=0;t<maskSize;t++)
		  {
		    u = i+int(forMask[t].x);
		    v = j+int(forMask[t].y);
		    w = k+int(forMask[t].z);
		    if((u>=0)&&(u<pRow)&&(v>=0)&&(v<pCol)&&(w>=0)&&(w<pSli)) 
		      {
			distOne = (float(u)-NEP(i,j,k).x);
			dist = distOne*distOne;
			distOne = (float(v)-NEP(i,j,k).y);
			dist = dist + distOne*distOne;
			distOne = (float(w)-NEP(i,j,k).z);
			dist = dist + distOne*distOne;
			//dist = sqrt((double(u)-NEP(i,j,k).x)+pow((double(v)-NEP(i,j,k).y),2)+pow((double(w)-NEP(i,j,k).z),2));
			if(DT(u,v,w)>dist)
			  {
			    DT(u,v,w) = dist;
			    NEP(u,v,w).x = NEP(i,j,k).x;
			    NEP(u,v,w).y = NEP(i,j,k).y;
			    NEP(u,v,w).z = NEP(i,j,k).z;			    
			    flagChange = true;
			  }
		      }
		  }
	      }
      //BackWard
      for(i=pRow-1;i>=0;i--)
	for(j=pCol-1;j>=0;j--)
	  for(k=pSli-1;k>=0;k--)
	    if(NEP(i,j,k).x!=-1)
	      {
		for(t=0;t<maskSize;t++)
		  {
		    u = i+int(backMask[t].x);
		    v = j+int(backMask[t].y);
		    w = k+int(backMask[t].z);
		    if((u>=0)&&(u<pRow)&&(v>=0)&&(v<pCol)&&(w>=0)&&(w<pSli)) 
		      {
			distOne = (float(u)-NEP(i,j,k).x);
			dist = distOne*distOne;
			distOne = (float(v)-NEP(i,j,k).y);
			dist = dist + distOne*distOne;
			distOne = (float(w)-NEP(i,j,k).z);
			dist = dist + distOne*distOne;
			//dist = sqrt(pow((double(u)-NEP(i,j,k).x),2)+pow((double(v)-NEP(i,j,k).y),2)+pow((double(w)-NEP(i,j,k).z),2));
			if(DT(u,v,w)>dist)
			  {
			    DT(u,v,w) = dist;
			    NEP(u,v,w).x = NEP(i,j,k).x;
			    NEP(u,v,w).y = NEP(i,j,k).y;
			    NEP(u,v,w).z = NEP(i,j,k).z;			    
			    flagChange = true;
			  }
		      }
		  }
	      } 
    }    
  free(forMask); 
  free(backMask);
  std::cout<<std::endl<<"Number of raster scans: "<<countTimes<<std::endl; 
}
