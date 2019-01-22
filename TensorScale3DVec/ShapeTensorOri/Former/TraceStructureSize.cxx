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
  int flagPre, flagPost;
  float preA, postA, dtThreshold, preDT, postDT;

  //Trace along both sides
  for(i=0;i<pRow-1;i++)
    for(j=0;j<pCol-1;j++)
      for(k=0;k<pSli-1;k++)
	{ 
	  oriX = vecOri(i,j,k).x;
	  oriY = vecOri(i,j,k).y;
	  oriZ = vecOri(i,j,k).z;

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
	      flagPre = 0; flagPost = 0;   
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
		  else 
		    if((traceDTValue<1.0)&&(traceGradient>acceptableGradientNearBoundary))
		      {
			localScale = float(a)*delta + traceDTValue/fabs(traceGradient);
			break;
		      }
		    /*
		    {
		      if (!flagPre &&(traceGradient < -0.1 && traceDTValue < dtThreshold))
			{
			  flagPre = 1;
			  preA = a;
			  preDT = traceDTValue;
			}
		      if(flagPre &&(traceGradient > 0.1 && traceDTValue > dtThreshold))
			{
			  flagPost = 1;
			  postA = a;
			  postDT = traceDTValue;
			}
		      if(flagPre && flagPost) 
			{
			  localScale = (preA+(postA-preA)*preDT/(postDT+preDT))*delta;
			  break;
			}
		    }
		    */
		}
	      else 
		{		  
		  if((cornerX<0)||(cornerY<0)||(cornerZ<0))
		    localScale = double(a)*delta;
		  else
		    localScale = double(a-1)*delta; 
		  //localScale = double(maxTrace)*delta;
		  break;
		}
	      preDTValue = traceDTValue;
	      a++;
	    }
	  //find sizeL
	  sizeL = localScale;
	  if(a==maxTrace) sizeL = double(maxTrace) * delta;

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
	      flagPre = 0; flagPost = 0;
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
		  else 
		    if((traceDTValue<1.0)&&(traceGradient<-acceptableGradientNearBoundary))
		      {
			localScale = float(a)*delta + traceDTValue/fabs(traceGradient);
			break;
		      }
		  /*
		    {
		      if (!flagPre &&(traceGradient < -0.1 && traceDTValue < dtThreshold))
			{
			  flagPre = 1;
			  preA = a;
			  preDT = traceDTValue;

			}
		      if(flagPre &&(traceGradient > 0.1 && traceDTValue > dtThreshold))
			{
			  flagPost = 1;
			  postA = a;
			  postDT = traceDTValue;

			}
		      if(flagPre && flagPost) 
			{
			  localScale = (preA+(postA-preA)*preDT/(postDT+preDT))*delta;
			  break;
			}
		    }
		  */
		}
	      else 
		{
		  if((cornerX<0)||(cornerY<0)||(cornerZ<0))
		    localScale = double(a)*delta;
		  else
		    localScale = double(a-1)*delta;
		  //localScale = double(maxTrace)*delta;
		  break;
		}
	      preDTValue = traceDTValue;
	      a++;
	    }
	  //find sizeR
	  sizeR = localScale;
	  if(a==maxTrace) sizeR= double(maxTrace) * delta;

	  //choose smaller one as size
	  if(sizeL<sizeR) localSize(i,j,k) = sizeL;
	  else  localSize(i,j,k) = sizeR;
	}
}
