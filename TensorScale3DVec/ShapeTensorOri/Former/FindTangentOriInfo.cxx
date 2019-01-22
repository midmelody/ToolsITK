void FindTangentOri(float *DT, VOXEL *vecOri, VOXEL *tangentOri1, VOXEL *tangentOri2)
{
  int i,j,k;//index
  float x,y,z;//orientation vector (z axis after rotation)
  float xx,xy,xz;//x axis after rotation
  float kVx, kVy, kVz; //vertical component
  float kHx, kHy, kHz; //horizontal component
  float xt1,yt1,zt1,xt2,yt2,zt2;//tangent orientation vector
  float dotProduct; //get cosine orientation
  float norm;//length of the vector
  float sin1,sin2,cos1,cos2;
  int a,b;//5*5 grid index
  VOXEL grid[5][5]; //converted coordinate
  VOXEL orthoSampleVec[5][5];
  vnl_matrix<float> samplePCA(8*2,3,0.0);
  float sampleValue[25];
  const int neighborSize = 2; //neighborhood Size*2+1
  int cornerX, cornerY, cornerZ;
  float weightX, weightY, weightZ;
  int t; //1d counter
  float graX,graY;
  float temp;
  int specX=40, specY=57, specZ=74;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{ 
	  //std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	  //shape vec at p - normalized as local z unit vector
	  x = vecOri(i,j,k).x;
	  y = vecOri(i,j,k).y;
	  z = vecOri(i,j,k).z;
	  /*////////////////////////////////////////////////////////////////////////////
	    if((i==specX)&&(j==specY)&&(k==specZ))
	    std::cout<<std::endl<<"Shape Tensor Vec:"<<std::endl<<x<<" "<<y<<" "<<z<<" "<<y/x<<std::endl<<std::endl;
	  */////////////////////////////////////////////////////////////////////////////
	  //Compute Rotation system
	  if((y!=0)||(z!=0))
	    {
	      sin1 = y/sqrt(y*y+z*z);
	      cos1 = z/sqrt(y*y+z*z);
	      sin2 = x/sqrt(x*x+y*y+z*z);
	      cos2 = sqrt(y*y+z*z)/sqrt(x*x+y*y+z*z);
	    }
	  else
	    {
	      sin1 = 1;
	      cos1 = 0;
	      sin2 = 1;
	      cos2 = 0;
	    }
	  //Get the shape tensor value at grid locations after rotation
	  /*////////////////////////////////////////////////////////////////////////////
	    if((i==specX)&&(j==specY)&&(k==specZ))
	    std::cout<<"Sample Vec:"<<std::endl;
	  */////////////////////////////////////////////////////////////////////////////
	  for(a=-neighborSize;a<=neighborSize;a++)
	    for(b=-neighborSize;b<=neighborSize;b++)
	      {
		//Get coordinate after rotation
		grid[a+neighborSize][b+neighborSize].x = cos2*a;
		grid[a+neighborSize][b+neighborSize].y = -sin1*sin2*a+cos1*b;
		grid[a+neighborSize][b+neighborSize].z = -cos1*sin2*a-sin1*b;
		/*////////////////////////////////////////////////////////////////////////////
		if((i==specX)&&(j==specY)&&(k==specZ)&&(a==1)&(b==0))
		  std::cout<<"X Vec:"<<std::endl<<grid[a+neighborSize][b+neighborSize].x<<" "<<grid[a+neighborSize][b+neighborSize].y<<" "<<grid[a+neighborSize][b+neighborSize].z<<" "<<std::endl;
		*/////////////////////////////////////////////////////////////////////////////
		grid[a+neighborSize][b+neighborSize].x = grid[a+neighborSize][b+neighborSize].x+i;
		grid[a+neighborSize][b+neighborSize].y = grid[a+neighborSize][b+neighborSize].y+j;
		grid[a+neighborSize][b+neighborSize].z = grid[a+neighborSize][b+neighborSize].z+k;
		
		//PCA from 8 neighbors to get value -- distance weighted
		cornerX = floor(grid[a+neighborSize][b+neighborSize].x);
		cornerY = floor(grid[a+neighborSize][b+neighborSize].y);
		cornerZ = floor(grid[a+neighborSize][b+neighborSize].z);
		weightX = grid[a+neighborSize][b+neighborSize].x - cornerX;
		weightY = grid[a+neighborSize][b+neighborSize].y - cornerY;
		weightZ = grid[a+neighborSize][b+neighborSize].z - cornerZ;
		if(cornerX<0) cornerX=0;
		if(cornerX>pRow-2) cornerX=pRow-2;
		if(cornerY<0) cornerY=0;
		if(cornerY>pCol-2) cornerY=pCol-2;
		if(cornerZ<0) cornerZ=0;
		if(cornerZ>pSli-2) cornerZ=pSli-2;
		samplePCA(0,0) = vecOri(cornerX,cornerY,cornerZ).x*(1-weightX)*(1-weightY)*(1-weightZ);
		samplePCA(0,1) = vecOri(cornerX,cornerY,cornerZ).y*(1-weightX)*(1-weightY)*(1-weightZ);
		samplePCA(0,2) = vecOri(cornerX,cornerY,cornerZ).z*(1-weightX)*(1-weightY)*(1-weightZ);
		samplePCA(1,0) = vecOri(cornerX+1,cornerY,cornerZ).x*weightX*(1-weightY)*(1-weightZ);
		samplePCA(1,1) = vecOri(cornerX+1,cornerY,cornerZ).y*weightX*(1-weightY)*(1-weightZ);
		samplePCA(1,2) = vecOri(cornerX+1,cornerY,cornerZ).z*weightX*(1-weightY)*(1-weightZ);
		samplePCA(2,0) = vecOri(cornerX,cornerY+1,cornerZ).x*(1-weightX)*weightY*(1-weightZ);
		samplePCA(2,1) = vecOri(cornerX,cornerY+1,cornerZ).y*(1-weightX)*weightY*(1-weightZ);
		samplePCA(2,2) = vecOri(cornerX,cornerY+1,cornerZ).z*(1-weightX)*weightY*(1-weightZ);
		samplePCA(3,0) = vecOri(cornerX,cornerY,cornerZ+1).x*(1-weightX)*(1-weightY)*weightZ;
		samplePCA(3,1) = vecOri(cornerX,cornerY,cornerZ+1).y*(1-weightX)*(1-weightY)*weightZ;
		samplePCA(3,2) = vecOri(cornerX,cornerY,cornerZ+1).z*(1-weightX)*(1-weightY)*weightZ;
		samplePCA(4,0) = vecOri(cornerX+1,cornerY+1,cornerZ).x*weightX*weightY*(1-weightZ);
		samplePCA(4,1) = vecOri(cornerX+1,cornerY+1,cornerZ).y*weightX*weightY*(1-weightZ);
		samplePCA(4,2) = vecOri(cornerX+1,cornerY+1,cornerZ).z*weightX*weightY*(1-weightZ);
		samplePCA(5,0) = vecOri(cornerX,cornerY+1,cornerZ+1).x*(1-weightX)*weightY*weightZ;
		samplePCA(5,1) = vecOri(cornerX,cornerY+1,cornerZ+1).y*(1-weightX)*weightY*weightZ;
		samplePCA(5,2) = vecOri(cornerX,cornerY+1,cornerZ+1).z*(1-weightX)*weightY*weightZ;
		samplePCA(6,0) = vecOri(cornerX+1,cornerY,cornerZ+1).x*weightX*(1-weightY)*weightZ;
		samplePCA(6,1) = vecOri(cornerX+1,cornerY,cornerZ+1).y*weightX*(1-weightY)*weightZ;
		samplePCA(6,2) = vecOri(cornerX+1,cornerY,cornerZ+1).z*weightX*(1-weightY)*weightZ;	
		samplePCA(7,0) = vecOri(cornerX+1,cornerY+1,cornerZ+1).x*weightX*weightY*weightZ;
		samplePCA(7,1) = vecOri(cornerX+1,cornerY+1,cornerZ+1).y*weightX*weightY*weightZ;
		samplePCA(7,2) = vecOri(cornerX+1,cornerY+1,cornerZ+1).z*weightX*weightY*weightZ;

		//Apply SVD to get sample value
		vnl_matrix<float> SVDV = vnl_svd<float> (samplePCA).V(); 
		orthoSampleVec[a+neighborSize][b+neighborSize].x = SVDV(0,0);
 		orthoSampleVec[a+neighborSize][b+neighborSize].y = SVDV(1,0); 
		orthoSampleVec[a+neighborSize][b+neighborSize].z = SVDV(2,0);

		//Generate unit vector
		norm = orthoSampleVec[a+neighborSize][b+neighborSize].x * orthoSampleVec[a+neighborSize][b+neighborSize].x;
		norm = norm + orthoSampleVec[a+neighborSize][b+neighborSize].y * orthoSampleVec[a+neighborSize][b+neighborSize].y;
 		norm = norm + orthoSampleVec[a+neighborSize][b+neighborSize].z * orthoSampleVec[a+neighborSize][b+neighborSize].z;
		norm = sqrt(norm);
		orthoSampleVec[a+neighborSize][b+neighborSize].x = orthoSampleVec[a+neighborSize][b+neighborSize].x/norm;
		orthoSampleVec[a+neighborSize][b+neighborSize].y = orthoSampleVec[a+neighborSize][b+neighborSize].y/norm;
		orthoSampleVec[a+neighborSize][b+neighborSize].z = orthoSampleVec[a+neighborSize][b+neighborSize].z/norm;
		/*////////////////////////////////////////////////////////////////////////////
		if((i==specX)&&(j==specY)&&(k==specZ))
		  std::cout<<orthoSampleVec[a+neighborSize][b+neighborSize].x<<" "<<orthoSampleVec[a+neighborSize][b+neighborSize].y<<" "<<orthoSampleVec[a+neighborSize][b+neighborSize].z<<std::endl;
		*/////////////////////////////////////////////////////////////////////////////
	      }

	  //From projection in the sample plane, get inplane orientation as graX and graY
	  /*////////////////////////////////////////////////////////////////////////////
	  if((i==specX)&&(j==specY)&&(k==specZ))
	    std::cout<<std::endl<<"Inplane:"<<std::endl;
	  */////////////////////////////////////////////////////////////////////////////

	  //Constract matrix for PCA
	  vnl_matrix<float> planePCA(25*2,2,0.0);
	  float tempX, tempY, tempZ; 
	  //x axis after rotation (unit vector)
	  tempX = 1;
	  tempY = 0;
	  tempZ = 0;
	  xx = cos2*tempX;
	  xy = -sin1*sin2*tempX+cos1*tempY;
	  xz = -cos1*sin2*tempX-sin1*tempY;
	  norm = xx*xx + xy*xy + xz*xz;
	  norm = sqrt(norm);
	  xx = xx/norm;
	  xy = xy/norm;
	  xz = xz/norm;
	  /*////////////////////////////////////////////////////////////////////////////
	  if((i==specX)&&(j==specY)&&(k==specZ))
	    std::cout<<"X Vec:"<<std::endl<<xx<<" "<<xy<<" "<<xz<<" "<<xy/xx<<" "<<std::endl;
	  */////////////////////////////////////////////////////////////////////////////  
	  t=0;
	  for(a=-neighborSize;a<=neighborSize;a++)
	    for(b=-neighborSize;b<=neighborSize;b++)
	      {
		dotProduct = orthoSampleVec[a+neighborSize][b+neighborSize].x * x;
		dotProduct = dotProduct +  orthoSampleVec[a+neighborSize][b+neighborSize].y * y;
		dotProduct = dotProduct +  orthoSampleVec[a+neighborSize][b+neighborSize].z * z;
		kVx = dotProduct * x;
		kVy = dotProduct * y;
		kVz = dotProduct * z;
		kHx = orthoSampleVec[a+neighborSize][b+neighborSize].x - kVx;
		kHy = orthoSampleVec[a+neighborSize][b+neighborSize].y - kVy;
		kHz = orthoSampleVec[a+neighborSize][b+neighborSize].z - kVz;
		/*////////////////////////////////////////////////////////////////////////////	
		if((i==specX)&&(j==specY)&&(k==specZ))
		  std::cout<<"horizontal: "<<kHx<<" "<<kHy<<" "<<kHz<<" "<<kHy/kHx<<" "<<std::endl;
		*/////////////////////////////////////////////////////////////////////////////
		norm = kHx*kHx + kHy*kHy + kHz*kHz;
		norm = sqrt(norm);
		if(norm==0) norm=1;
		dotProduct = kHx*xx/norm + kHy*xy/norm + kHz*xz/norm;
		/*////////////////////////////////////////////////////////////////////////////
		if((i==specX)&&(j==specY)&&(k==specZ))
		  std::cout<<"horizontal: "<<dotProduct<<std::endl;
		*/////////////////////////////////////////////////////////////////////////////
		//assign sample value
		planePCA(t,0) = norm*dotProduct;
		planePCA(t,1) = norm*(1-sqrt(dotProduct*dotProduct));
		planePCA(49-t,0) = norm*dotProduct;
		planePCA(49-t,1) = norm*(1-sqrt(dotProduct*dotProduct));
		t = t+1;
	      }
	  /*////////////////////////////////////////////////////////////////////////////
	  if((i==specX)&&(j==specY)&&(k==specZ))
	    std::cout<<"plane PCA:"<<std::endl<<planePCA;
	  */////////////////////////////////////////////////////////////////////////////

	  //Analyse the direction
	  vnl_matrix<float> SVDplane = vnl_svd<float> (planePCA).V(); 
	  if(absl(SVDplane(0,0))<absl(SVDplane(1,0)))
	    {
	      graX = SVDplane(0,0);
	      graY = SVDplane(1,0); 
	    }
	  else
	    {
	      graX = SVDplane(1,0);
	      graY = SVDplane(0,0); 
	    }
	  /*////////////////////////////////////////////////////////////////////////////
	  if((i==specX)&&(j==specY)&&(k==specZ))
	    std::cout<<"Inplane direction: "<<graX<<" "<<graY<<std::endl;
	  */////////////////////////////////////////////////////////////////////////////

	  //Apply 3D rotation for the orientation vector (graX, graY, 0)
	  //save the normalized result in tangentOri1
	  tangentOri1(i,j,k).x = cos2*graX;
	  tangentOri1(i,j,k).y = -sin1*sin2*graX+cos1*graY;
	  tangentOri1(i,j,k).z = -cos1*sin2*graX-sin1*graY;
	  xt1 = tangentOri1(i,j,k).x;
	  yt1 = tangentOri1(i,j,k).y;
	  zt1 = tangentOri1(i,j,k).z;
	  if(absl(xt1)<0.0001) xt1 = 0;
	  if(absl(yt1)<0.0001) yt1 = 0;
	  if(absl(zt1)<0.0001) zt1 = 0;
	  //use upper space 
	  if(zt1<0)
	    {
	      xt1 = -xt1;
	      yt1 = -yt1;
	      zt1 = -zt1;
	    }
	  /*////////////////////////////////////////////////////////////////////////////
	  if((i==specX)&&(j==specY)&&(k==specZ))
	    std::cout<<xt1<<" "<<yt1<<" "<<zt1<<std::endl;
	  */////////////////////////////////////////////////////////////////////////////
	  norm = sqrt(xt1*xt1+yt1*yt1+zt1*zt1);
	  tangentOri1(i,j,k).x = xt1/norm;
	  tangentOri1(i,j,k).y = yt1/norm;
	  tangentOri1(i,j,k).z = zt1/norm;

	  //record the normalized cross product in tangentOri2
	  tangentOri2(i,j,k).x = y*zt1-z*yt1;
	  tangentOri2(i,j,k).y = z*xt1-x*zt1;
	  tangentOri2(i,j,k).z = x*yt1-y*xt1;
	  xt2 = tangentOri2(i,j,k).x;
	  yt2 = tangentOri2(i,j,k).y;
	  zt2 = tangentOri2(i,j,k).z;
	  if(absl(xt2)<0.0001) xt2 = 0;
	  if(absl(yt2)<0.0001) yt2 = 0;
	  if(absl(zt2)<0.0001) zt2 = 0;
	  //use upper space 
	  if(zt2<0)
	    {
	      xt2 = -xt2;
	      yt2 = -yt2;
	      zt2 = -zt2;
	    }
	  /*////////////////////////////////////////////////////////////////////////////
	  if((i==specX)&&(j==specY)&&(k==specZ))
	    std::cout<<xt2<<" "<<yt2<<" "<<zt2<<std::endl;
	  */////////////////////////////////////////////////////////////////////////////
	  norm = sqrt(xt2*xt2+yt2*yt2+zt2*zt2);
	  tangentOri2(i,j,k).x = xt2/norm;
	  tangentOri2(i,j,k).y = yt2/norm;
	  tangentOri2(i,j,k).z = zt2/norm;	  
	}
} 



/*

  //----------------------Smooth the vector field ----------------------------------------------
  //Gradient magnitude weighted (center line attenuation)
  int x,y,z;
  int smoothSize=5; //neighborhood size: 2n+1 * 2n+1 * 2n+1
  VOXEL *smoothTangOri1;
  smoothTangOri1 = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);
  VOXEL *smoothTangOri2;
  smoothTangOri2 = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);

  //First,  record orientation pairs in matrix
  //Second, add additional 0s of the neighborhood size to the matrix
  //Third,  compute SVD for orientation result
  vnl_matrix<float> graPairs((smoothSize*2+1)*(smoothSize*2+1)*(smoothSize*2+1)*2,3,0.0);  
  float tempLength;
  for(i=0;i<pRow;i++)
    {
      for(j=0;j<pCol;j++)
	for(k=0;k<pSli;k++)
	  {
	    for(t=0;t<(smoothSize*2+1)*(smoothSize*2+1)*(smoothSize*2+1)*2;t++)
	      {
		graPairs(t,0)=0.0;
		graPairs(t,1)=0.0;
		graPairs(t,2)=0.0;
	      }
	    t=0;
	    //compute gradient magnitude and record in pairs
	    for(x=-smoothSize;x<=smoothSize;x++)
	      for(y=-smoothSize;y<=smoothSize;y++)	  
		for(z=-smoothSize;z<=smoothSize;z++)	  
		  if(((i+x)>=0)&&((i+x)<pRow)&&((j+y)>=0)&&((j+y)<pCol)&&((k+z)>=0)&&((k+z)<pSli))
		    {
		      tempLength = GradientDTX(i+x,j+y,k+z)*GradientDTX(i+x,j+y,k+z);
		      tempLength = tempLength + GradientDTY(i+x,j+y,k+z)*GradientDTY(i+x,j+y,k+z);
		      tempLength = tempLength + GradientDTZ(i+x,j+y,k+z)*GradientDTZ(i+x,j+y,k+z);
		      graPairs(t,0) = GradientDTX(i+x,j+y,k+z)*tempLength;
		      graPairs(t,1) = GradientDTY(i+x,j+y,k+z)*tempLength;
		      graPairs(t,2) = GradientDTZ(i+x,j+y,k+z)*tempLength;
		      t = t+1;
		    }
	    //SVD
	    vnl_matrix<float> SVDV = vnl_svd<float> (graPairs).V(); 
	    SmoothGraDTX(i,j,k) = SVDV(0,0);
	    SmoothGraDTY(i,j,k) = SVDV(1,0);
	    SmoothGraDTZ(i,j,k) = SVDV(2,0);
	    if(signf(SmoothGraDTX(i,j,k))!=signf(GradientDTX(i,j,k))) 
	      SmoothGraDTX(i,j,k) = - SmoothGraDTX(i,j,k);
	    if(signf(SmoothGraDTY(i,j,k))!=signf(GradientDTY(i,j,k))) 
	      SmoothGraDTY(i,j,k) = - SmoothGraDTY(i,j,k);
	    if(signf(SmoothGraDTZ(i,j,k))!=signf(GradientDTZ(i,j,k))) 
	      SmoothGraDTZ(i,j,k) = - SmoothGraDTZ(i,j,k);
	  }
      std::cout<<i<<std::endl;
    }








	  if(xt1==0)
	    {
	      if(yt1*zt1>=0)
		{
		  yt1 = absl(yt1);
		  zt1 = absl(zt1);
		}
	      else
		{
		  yt1 = absl(yt1);
		  zt1 = -absl(zt1);
		}
	    }
	  else if(yt1==0)
	    {
	      if(xt1*zt1>=0)
		{
		  xt1 = absl(xt1);
		  zt1 = absl(zt1);
		}
	      else
		{
		  xt1 = absl(xt1);
		  zt1 = -absl(zt1);
		}
	    }
	  else if(zt1==0)
	    {
	      if(xt1*yt1>=0)
		{
		  xt1 = absl(xt1);
		  yt1 = absl(yt1);
		}
	      else
		{
		  xt1 = absl(xt1);
		  yt1 = -absl(yt1);
		}
	    }
	  else 




	  if(xt2==0)
	    {
	      if(yt2*zt2>0)
		{
		  yt2 = absl(yt2);
		  zt2 = absl(zt2);
		}
	      else
		{
		  yt2 = absl(yt2);
		  zt2 = -absl(zt2);
		}
	    }
	  else if(yt2==0)
	    {
	      if(xt2*zt2>0)
		{
		  xt2 = absl(xt2);
		  zt2 = absl(zt2);
		}
	      else
		{
		  xt2 = absl(xt2);
		  zt2 = -absl(zt2);
		}
	    }
	  else if(zt2==0)
	    {
	      if(xt2*yt2>0)
		{
		  xt2 = absl(xt2);
		  yt2 = absl(yt2);
		}
	      else
		{
		  xt2 = absl(xt2);
		  yt2 = -absl(yt2);
		}
	    }
	  else 
*/
