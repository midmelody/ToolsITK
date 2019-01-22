//Compute zero crossing location, stored in zeroCrossing list and return the number of zero crossings 
int ComputeZeroCrossing(float *LoGImageHigh,float *LoGImageLow,float *gradientImage, float graThreHigh, float graThreLow, VOXEL *zeroCrossing)
{
  int numberOfZeroCrossing=0;
  VOXEL group[8]; //Find the pos/neg group with less members, allocate max space
  int memberNum; //number of members
  int dist, minDist;//city block distance between 2 members
  int plus, minus; //number of positive and negative value
  int i,j,k;
  int a,b,c;
  int u,v,w;
  int p,q;
  int t;
  bool edgeLow, edgeHigh;
  float edgeNum; //number of possible edges
  float refValue, curValue;
  float edge; //store edge position
  int   refX, refY, refZ; // Point from "smaller group"
  float edgeLowX, edgeLowY, edgeLowZ;
  float edgeHighX, edgeHighY, edgeHighZ;
  float gradient;
  float i1,i2,j1,j2,w1,w2;

  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  //std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	  edgeLow = false;
	  edgeHigh = false;
	  /////////////////////////////////////////////////
	  //find edge in case of low:
	  //First judge if there may exist an edge point 
	  plus = 0;
	  minus = 0;
	  for(a=0;a<=1;a++)
	    for(b=0;b<=1;b++)
	      for(c=0;c<=1;c++)
		{
		  if(LoGImageLow(i+a,j+b,k+c)>0) plus++;
		  if(LoGImageLow(i+a,j+b,k+c)<0) minus++;  
		  /*//////////////////////////////////////////////////
		  if((i==68)&&(j==15)&&((k==1)||(k==2)))
		    std::cout<<"k= "<<k<<" :"<<LoGImageLow(i+a,j+b,k+c)<<std::endl;
		  *///////////////////////////////////////////////////
		}
	  //if sign change exists, there is potential edge
	  if(plus*minus!=0) edgeLow = true;
	  if(plus+minus!=8) edgeLow = false;

	  //if low edge may exist, make further validation
	  if(edgeLow)
	    {
	      //initialize temp parameters
	      //Find smaller group 
	      //assign space and record position
	      if(plus>=minus) 
		{
		  memberNum = minus; 
		  p=0; 
		  for(a=0;a<=1;a++)
		    for(b=0;b<=1;b++)
		      for(c=0;c<=1;c++)
			if(LoGImageLow(i+a,j+b,k+c)<0) 
			  { 
			    group[p].x = a;
			    group[p].y = b;
			    group[p].z = c;
			    p = p+1;
			  } 
		}
	      else
		{
		  memberNum = plus; 
		  p=0; 
		  for(a=0;a<=1;a++)
		    for(b=0;b<=1;b++)
		      for(c=0;c<=1;c++)
			if(LoGImageLow(i+a,j+b,k+c)>0) 
			  { 
			    group[p].x = a;
			    group[p].y = b;
			    group[p].z = c;
			    p = p+1;
			  } 
		}

	      //verify edge possibilities
	      //in the "group" list, there should not be any isolated point.
	      //isolated point: dist to all other points are >1
	      if(memberNum>1)
		for(p=0;p<memberNum;p++)
		  {
		    minDist = 100;
		    for(q=0;q<memberNum;q++)
		      if(q!=p)
			{
			  dist=abs(group[p].x-group[q].x)+abs(group[p].y-group[q].y)+abs(group[p].z-group[q].z);
			  if(minDist>dist) 
			    minDist=dist;
			} 
		    if(minDist>1)
		      edgeLow = false;
		  }
	    }

	  //If edge do exist, find it
	  if(edgeLow)
	    {
	      //initialize temp parameters
	      edgeNum = 0;
	      gradient = 0;
	      edgeLowX = 0;
	      edgeLowY = 0;
	      edgeLowZ = 0;	    
	      //go through the list of starting points
	      for(p=0;p<memberNum;p++)
		{
		  refX = int(group[p].x);
		  refY = int(group[p].y);
		  refZ = int(group[p].z);
		  //std::cout<<refX<<" "<<refY<<" "<<refZ<<std::endl;
		  //go through 8 points in the cube that has an edge in between with Ref
		  for(a=0;a<=1;a++)
		    for(b=0;b<=1;b++)
		      for(c=0;c<=1;c++)
			if(LoGImageLow(i+a,j+b,k+c)*LoGImageLow(i+refX,j+refY,k+refZ)<0)
			  {
			    //find edge
			    edgeNum = edgeNum + 1;
			    refValue = absl(LoGImageLow(i+refX,j+refY,k+refZ));
			    curValue = absl(LoGImageLow(i+a,j+b,k+c));
			    edge = refValue/(refValue+curValue);
			    //add contributions
			    edgeLowX = edgeLowX + refX + pow(-1,refX)*abs(a-refX)*edge;
			    edgeLowY = edgeLowY + refY + pow(-1,refY)*abs(b-refY)*edge;
			    edgeLowZ = edgeLowZ + refZ + pow(-1,refZ)*abs(c-refZ)*edge;
			  }		
		}
	      edgeLowX = edgeLowX/edgeNum;
	      edgeLowY = edgeLowY/edgeNum;
	      edgeLowZ = edgeLowZ/edgeNum;
	      //3D trilinear interpolation for gradient
	      i1 = gradientImage(i,j,k)*(1-edgeLowZ) + gradientImage(i,j,k+1)*edgeLowZ;
	      i2 = gradientImage(i,j+1,k)*(1-edgeLowZ) + gradientImage(i,j+1,k+1)*edgeLowZ;
	      j1 = gradientImage(i+1,j,k)*(1-edgeLowZ) + gradientImage(i+1,j,k+1)*edgeLowZ;
	      j2 = gradientImage(i+1,j+1,k)*(1-edgeLowZ) + gradientImage(i+1,j+1,k+1)*edgeLowZ;
	      w1 = i1*(1-edgeLowY) + i2*edgeLowY;
	      w2 = j1*(1-edgeLowY) + j2*edgeLowY;
	      gradient = w1*(1-edgeLowX) + w2*edgeLowX;
	      //check gradient to discard weak edge
	      if(gradient<graThreLow)
		{
		  edgeLow = false;
		  edgeLowX = 0;
		  edgeLowY = 0;
		  edgeLowZ = 0;
		}
	    }
	  //std::cout<<"here!!!!!!!!!!!!"<<std::endl;
	  
	  /////////////////////////////////////////////////////////////////////////
	  //find edge in case of high:
	  //First judge if there may exist an edge point 
	  plus = 0;
	  minus = 0;
	  for(a=0;a<=1;a++)
	    for(b=0;b<=1;b++)
	      for(c=0;c<=1;c++)
		{
		  if(LoGImageHigh(i+a,j+b,k+c)>0) plus++;
		  if(LoGImageHigh(i+a,j+b,k+c)<0) minus++; 
		  /*//////////////////////////////////////////////////
		  if((i==68)&&(j==15)&&((k==1)||(k==2)))
		    std::cout<<"k= "<<k<<" :"<<LoGImageHigh(i+a,j+b,k+c)<<std::endl;
		  *///////////////////////////////////////////////////
		}
	  //if sign change exists, there is potential edge
	  if(plus*minus!=0) edgeHigh = true;
	  if(plus+minus!=8) edgeHigh = false;
	  //if high edge may exist, make further validation
	  if(edgeHigh)
	    {
	      //initialize temp parameters
	      //Find smaller group 
	      //assign space and record position
	      if(plus>=minus) 
		{
		  memberNum = minus; 
		  p=0; 
		  for(a=0;a<=1;a++)
		    for(b=0;b<=1;b++)
		      for(c=0;c<=1;c++)
			if(LoGImageHigh(i+a,j+b,k+c)<0) 
			  { 
			    group[p].x = a;
			    group[p].y = b;
			    group[p].z = c;
			    p = p+1;
			  } 
		}
	      else
		{		
		  memberNum = plus; 
		  p=0; 
		  for(a=0;a<=1;a++)
		    for(b=0;b<=1;b++)
		      for(c=0;c<=1;c++)
			if(LoGImageHigh(i+a,j+b,k+c)>0) 
			  { 
			    group[p].x = a;
			    group[p].y = b;
			    group[p].z = c;
			    p = p+1;
			  } 
		}
	      //verify edge possibilities
	      //in the "group" list, there should not be any isolated point.
	      //isolated point: dist to all other points are >1
	      if(memberNum>1)
		for(p=0;p<memberNum;p++)
		  {
		    minDist = 100;
		    for(q=0;q<memberNum;q++)
		      if(q!=p)
			{
			  dist=abs(group[p].x-group[q].x)+abs(group[p].y-group[q].y)+abs(group[p].z-group[q].z);
			  if(minDist>dist) 
			    minDist=dist;
			} 
		    if(minDist>1)
		      edgeHigh = false;
		  }
	    }
	  //If edge do exist, find it
	  if(edgeHigh)
	    {
	      //initialize temp parameters
	      edgeNum = 0;
	      gradient = 0;
	      edgeHighX = 0;
	      edgeHighY = 0;
	      edgeHighZ = 0;	    
	      //go through the list of starting points
	      for(p=0;p<memberNum;p++)
		{
		  refX = int(group[p].x);
		  refY = int(group[p].y);
		  refZ = int(group[p].z);
		  //go through 8 points in the cube that has an edge in between with Ref
		  for(a=0;a<=1;a++)
		    for(b=0;b<=1;b++)
		      for(c=0;c<=1;c++)
			if(LoGImageHigh(i+a,j+b,k+c)*LoGImageHigh(i+refX,j+refY,k+refZ)<0)
			  {
			    //find edge
			    edgeNum = edgeNum + 1;
			    refValue = absl(LoGImageHigh(i+refX,j+refY,k+refZ));
			    curValue = absl(LoGImageHigh(i+a,j+b,k+c));
			    edge = refValue/(refValue+curValue);
			    //add contributions
			    edgeHighX = edgeHighX + refX + pow(-1,refX)*abs(a-refX)*edge;
			    edgeHighY = edgeHighY + refY + pow(-1,refY)*abs(b-refY)*edge;
			    edgeHighZ = edgeHighZ + refZ + pow(-1,refZ)*abs(c-refZ)*edge;
			  }		
		}
	      edgeHighX = edgeHighX/edgeNum;
	      edgeHighY = edgeHighY/edgeNum;
	      edgeHighZ = edgeHighZ/edgeNum;
	      //3D trilinear interpolation for gradient
	      i1 = gradientImage(i,j,k)*(1-edgeHighZ) + gradientImage(i,j,k+1)*edgeHighZ;
	      i2 = gradientImage(i,j+1,k)*(1-edgeHighZ) + gradientImage(i,j+1,k+1)*edgeHighZ;
	      j1 = gradientImage(i+1,j,k)*(1-edgeHighZ) + gradientImage(i+1,j,k+1)*edgeHighZ;
	      j2 = gradientImage(i+1,j+1,k)*(1-edgeHighZ) + gradientImage(i+1,j+1,k+1)*edgeHighZ;
	      w1 = i1*(1-edgeHighY) + i2*edgeHighY;
	      w2 = j1*(1-edgeHighY) + j2*edgeHighY;
	      gradient = w1*(1-edgeHighX) + w2*edgeHighX;
	      //check gradient to discard weak edge
	      if(gradient<graThreHigh)
		{
		  edgeHigh = false;
		  edgeHighX = 0;
		  edgeHighY = 0;
		  edgeHighZ = 0;
		}
	    }

	  /*//////////////////////////////////////////////////
	  if((i==15)&&(j==68)&&((k==1)||(k==2)))
	    std::cout<<edgeHigh<<" "<<edgeLow<<std::endl;
	  *///////////////////////////////////////////////////
	  
	  ////////////////////////////////////////////////////////////////////////////////////
	  //Get final edge situation
	  if(edgeHigh==true)
	    {
	      //both edges exist or only high exists, keep the high
	      zeroCrossing[numberOfZeroCrossing].x = i + edgeHighX;
	      zeroCrossing[numberOfZeroCrossing].y = j + edgeHighY;
	      zeroCrossing[numberOfZeroCrossing].z = k + edgeHighZ;
	      numberOfZeroCrossing = numberOfZeroCrossing+1;	    
	    }
	  else if(edgeLow==true)
	    {
	      zeroCrossing[numberOfZeroCrossing].x = i + edgeLowX;
	      zeroCrossing[numberOfZeroCrossing].y = j + edgeLowY;
	      zeroCrossing[numberOfZeroCrossing].z = k + edgeLowZ;
	      numberOfZeroCrossing = numberOfZeroCrossing+1;	    
	    }
	}

  //Drop single point
  //record number of edges around certain neighbor
  int *edgeCount;
  edgeCount =(int *)malloc(sizeof(int) * pRow * pCol * pSli);   
  //record edge position
  VOXEL *edgePosition;
  edgePosition =(VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli); 
  //initialize invalid position and no count  
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  edgePosition(i,j,k).x = -1;
	  edgePosition(i,j,k).y = -1;
	  edgePosition(i,j,k).z = -1;
	  edgeCount(i,j,k) = 0;
	}

  //for all edge points, record position and set count to 1
  for(t=0;t<numberOfZeroCrossing;t++)
    {
      i = floor(zeroCrossing[t].x);
      j = floor(zeroCrossing[t].y);
      k = floor(zeroCrossing[t].z);
      if((i>=0)&&(i<pRow)&&(j>=0)&&(j<pCol)&&(k>=0)&&(k<pSli))
	{
	  edgePosition(i,j,k).x = zeroCrossing[t].x;
	  edgePosition(i,j,k).y = zeroCrossing[t].y;
	  edgePosition(i,j,k).z = zeroCrossing[t].z;
	  edgeCount(i,j,k) = 1;
	}
    }
 
  //compute edge number in neighborhood
  int accumEdge = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  if(edgeCount(i,j,k)>0)
	    {
	      accumEdge = 0;
	      //test 5*5*5 neigborhood
	      for(a=-2;a<=2;a++)
		for(b=-2;b<=2;b++)
		  for(c=-2;c<=2;c++)	  
		    {
		      u = i+a;
		      v = j+b;
		      w = k+c;
		      if(u<0) u=0;
		      if(u>pRow-1) u=pRow-1;
		      if(v<0) v=0; 
		      if(v>pCol-1) v=pCol-1;
		      if(w<0) w=0;
		      if(w>pSli-1) w=pSli-1;
		      if(edgeCount(u,v,w)>0)
			accumEdge++;
		    }
	      edgeCount(i,j,k) = accumEdge;
	    }
	}
 	 
  //drop edges isolated
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  if( edgeCount(i,j,k) < 10) 
	    {
	      edgePosition(i,j,k).x = -1;
	      edgePosition(i,j,k).y = -1;
	      edgePosition(i,j,k).z = -1;
	      edgeCount(i,j,k) = 0;
	    }
	}

  /*/////////////////////////////////////////////////
  std::cout<<edgePosition(15,68,1).x<<" "<<edgePosition(15,68,1).y<<" "<<edgePosition(15,68,1).z<<std::endl;
  std::cout<<edgePosition(15,68,2).x<<" "<<edgePosition(15,68,2).y<<" "<<edgePosition(15,68,2).z<<std::endl;
  *//////////////////////////////////////////////////

  //re-generate edge point list
  t = 0;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  if(edgeCount(i,j,k)>0)
	    {
	      zeroCrossing[t].x = edgePosition(i,j,k).x;
	      zeroCrossing[t].y = edgePosition(i,j,k).y;
	      zeroCrossing[t].z = edgePosition(i,j,k).z;	    
	      t = t+1;
	    }
	}

  numberOfZeroCrossing = t;
  free(edgeCount);
  free(edgePosition);
  return(numberOfZeroCrossing);
}
 
