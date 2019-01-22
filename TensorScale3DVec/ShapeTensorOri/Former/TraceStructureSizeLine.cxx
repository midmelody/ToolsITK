void TraceStructureSize(float *DT, VOXEL *vecOri, float *localSize, int maxTrace)
{
  int i,j,k;
  int a; //count for tracing
  float delta = 0.5; //sample interval
  float oriX,oriY,oriZ;
  float currentX, currentY, currentZ; //float index of pixel 
  int cornerX, cornerY, cornerZ; //corner pixel index for bilinear interpolation
  float weightX, weightY, weightZ; //weight for bilinear interpolation
  float sizeL, sizeR; //resulting structure size of left and right trace line
  //temp vector to store DT values along sample line
  float traceDTValue, preDTValue, traceGradient, localScale;
  float errorTolerance  = 0.05;
  float acceptableGradientNearBoundary = 0.1;
  float i1,i2,j1,j2,w1,w2;

  //Trace along both sides
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{ 
	  /*	   
	  ///////////////////////////////////////////
	  localSize(i,j,k) = sqrt(DT(i,j,k));
	}
  int xCoor = 39;
  int yCoor = 56;
  int zCoor = 74;
  for(i=xCoor;i<xCoor+1;i++)
    for(j=yCoor;j<yCoor+1;j++)
      for(k=zCoor;k<zCoor+1;k++)
	{
	  ///////////////////////////////////////////
	  */

	  oriX = vecOri(i,j,k).x;
	  oriY = vecOri(i,j,k).y;
	  oriZ = vecOri(i,j,k).z;

	  //////////////////////////////////////////////	  
	  //std::cout<<oriX<<" "<<oriY<<" "<<oriZ<<std::endl;
	  //////////////////////////////////////////////

	  //Get value for trace line:
	  //--------------------------------------------------------------
	  //left half line
	  a = 0;
	  traceDTValue = 10000;
	  preDTValue = DT(i,j,k);
	  while(a<maxTrace)
	    {
	      currentX = i + float(a) * delta * oriX;
	      currentY = j + float(a) * delta * oriY;
	      currentZ = k + float(a) * delta * oriZ;
	      cornerX = floor(currentX);
	      cornerY = floor(currentY);
	      cornerZ = floor(currentZ);
	      weightX = currentX - cornerX;
	      weightY = currentY - cornerY;
	      weightZ = currentZ - cornerZ;
      
	      if(cornerX>=0&&(cornerX<pRow-1)&&(cornerY>=0)&&(cornerY<pCol-1)&&(cornerZ>=0)&&(cornerZ<pSli-1))
		{
		  //3D trilinear interpolation for sample value
		  i1 = DT(cornerX,cornerY,cornerZ)*(1-weightZ) + DT(cornerX,cornerY,cornerZ+1)*weightZ;
		  i2 = DT(cornerX,cornerY+1,cornerZ)*(1-weightZ) + DT(cornerX,cornerY+1,cornerZ+1)*weightZ;
		  j1 = DT(cornerX+1,cornerY,cornerZ)*(1-weightZ) + DT(cornerX+1,cornerY,cornerZ+1)*weightZ;
		  j2 = DT(cornerX+1,cornerY+1,cornerZ)*(1-weightZ) + DT(cornerX+1,cornerY+1,cornerZ+1)*weightZ;
		  w1 = i1*(1-weightY) + i2*weightY;
		  w2 = j1*(1-weightY) + j2*weightY;
		  traceDTValue = w1*(1-weightX) + w2*weightX;
		  traceGradient = (preDTValue - traceDTValue)/delta;
		  traceGradient = absl(traceGradient);
		  if(traceGradient > 1.0-errorTolerance) 
		    {
		      localScale = traceDTValue+float(a)*delta;
		      break;
		    }
		  else if((traceDTValue<1.0)&&(traceGradient>acceptableGradientNearBoundary))
		    {
		      localScale = float(a)*delta + traceDTValue/fabs(traceGradient);
		      break;
		    }
		  //////////////////////////////////
		  //localSize(cornerX,cornerY,cornerZ) = 3;
		  //////////////////////////////////
		}
	      else 
		{
		  localScale = double(a)*delta;
		  break;
		}
	      preDTValue = traceDTValue;
	      a++;
	    }
	  //find sizeL
	  sizeL = localScale;
	  if(a==maxTrace) sizeL = maxTrace * delta;

	  //////////////////////////////////////////////////////////
	  /*
	  currentX = i + float(sizeL) * oriX;
	  currentY = j + float(sizeL) * oriY;
	  currentZ = k + float(sizeL) * oriZ;
	  cornerX = floor(currentX);
	  cornerY = floor(currentY);
	  cornerZ = floor(currentZ);
	  localSize(cornerX,cornerY,cornerZ) = 10;	
	  */
	  //////////////////////////////////////////////////////////
	  
	  //------------------------------------------------------------------------------  
	  //right half line 
	  a = 0;
	  traceDTValue = 10000;
	  preDTValue = DT(i,j,k);
	  while(a<maxTrace)
	    {
	      currentX = i - float(a) * delta * oriX;
	      currentY = j - float(a) * delta * oriY;
	      currentZ = k - float(a) * delta * oriZ;
	      cornerX = floor(currentX);
	      cornerY = floor(currentY);
	      cornerZ = floor(currentZ);
	      weightX = currentX - cornerX;
	      weightY = currentY - cornerY;
	      weightZ = currentZ - cornerZ;
     
	      if(cornerX>=0&&(cornerX<pRow-1)&&(cornerY>=0)&&(cornerY<pCol-1)&&(cornerZ>=0)&&(cornerZ<pSli-1))
		{
		  //3D trilinear interpolation for gradient
		  i1 = DT(cornerX,cornerY,cornerZ)*(1-weightZ) + DT(cornerX,cornerY,cornerZ+1)*weightZ;
		  i2 = DT(cornerX,cornerY+1,cornerZ)*(1-weightZ) + DT(cornerX,cornerY+1,cornerZ+1)*weightZ;
		  j1 = DT(cornerX+1,cornerY,cornerZ)*(1-weightZ) + DT(cornerX+1,cornerY,cornerZ+1)*weightZ;
		  j2 = DT(cornerX+1,cornerY+1,cornerZ)*(1-weightZ) + DT(cornerX+1,cornerY+1,cornerZ+1)*weightZ;
		  w1 = i1*(1-weightY) + i2*weightY;
		  w2 = j1*(1-weightY) + j2*weightY;
		  traceDTValue = w1*(1-weightX) + w2*weightX;
		  traceGradient = (preDTValue - traceDTValue)/delta;		  
		  if(traceGradient < -1.0+errorTolerance) 
		    {
		      localScale = traceDTValue+float(a)*delta;
		      break;
		    }
		  else if((traceDTValue<1.0)&&(traceGradient<-acceptableGradientNearBoundary))
		    {
		      localScale = float(a)*delta + traceDTValue/fabs(traceGradient);
		      break;
		    }
		  //////////////////////////////////
		  //localSize(cornerX,cornerY,cornerZ) = 3;
		  //////////////////////////////////
		}
	      else 
		{
		  localScale = double(a)*delta;
		  break;
		}
	      preDTValue = traceDTValue;
	      a++;
	    }
	  //find sizeR
	  sizeR = localScale;
	  if(a==maxTrace) sizeR=maxTrace * delta;
	  //////////////////////////////////////////////////////////
	  /*
	  currentX = i - float(sizeR) * oriX;
	  currentY = j - float(sizeR) * oriY;
	  currentZ = k - float(sizeR) * oriZ;
	  cornerX = floor(currentX);
	  cornerY = floor(currentY);
	  cornerZ = floor(currentZ);
	  localSize(cornerX,cornerY,cornerZ) = 10;	
	  */
	  //////////////////////////////////////////////////////////
	  //choose smaller one as size
	  if(sizeL<sizeR) localSize(i,j,k) = sizeL;
	  else  localSize(i,j,k) = sizeR;
	  //////////////////////
	  //localSize(i,j,k) = 8;
	  //////////////////////
	  //localSize(i,j,k) = sqrt(localSize(i,j,k));
	}
}
