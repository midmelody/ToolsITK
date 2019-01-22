void FindTangentOri(float *DT, VOXEL *vecOri, VOXEL *tangentOri1, VOXEL *tangentOri2)
{
  int i,j,k;//index
  float x,y,z;//orientation vector (z axis after rotation)
  float xx,xy,xz;//x axis after rotation
  float yx,yy,yz;//y axis after rotation
  float kVx, kVy, kVz; //vertical component
  float kHx, kHy, kHz; //horizontal component
  float xt1,yt1,zt1,xt2,yt2,zt2;//tangent orientation vector
  float dotProduct; //get cosine orientation
  float norm;//length of the vector
  float sin1,sin2,cos1,cos2;
  int a,b;//grid index
  int c;
  int gridSize;
  VOXEL grid[7][7]; //converted coordinate
  VOXEL orthoSampleVec[7][7];
  vnl_matrix<float> samplePCA(8*2,3,0.0);
  float sampleValue[49];
  const int neighborSize = 3; //neighborhood Size*2+1
  int cornerX, cornerY, cornerZ;
  float weightX, weightY, weightZ;
  int t; //1d counter
  float graX,graY;
  float temp;
  std::cout<<"Find direction of principle curvature: "<<std::endl;
  for(i=0;i<pRow;i++)
    {
      std::cout<<"\r";
      std::cout<<int((i+1)*100/pRow)<<"%"<<std::flush;
      for(j=0;j<pCol;j++)
	for(k=0;k<pSli;k++)
	  { 
	    //Shape normal vector info
	    x = vecOri(i,j,k).x;
	    y = vecOri(i,j,k).y;
	    z = vecOri(i,j,k).z;
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
	    for(a=-neighborSize;a<=neighborSize;a++)
	      for(b=-neighborSize;b<=neighborSize;b++)
		{
		  //Get coordinate after rotation
		  grid[a+neighborSize][b+neighborSize].x = cos2*a;
		  grid[a+neighborSize][b+neighborSize].y = -sin1*sin2*a+cos1*b;
		  grid[a+neighborSize][b+neighborSize].z = -cos1*sin2*a-sin1*b;
		  grid[a+neighborSize][b+neighborSize].x = grid[a+neighborSize][b+neighborSize].x+i;
		  grid[a+neighborSize][b+neighborSize].y = grid[a+neighborSize][b+neighborSize].y+j;
		  grid[a+neighborSize][b+neighborSize].z = grid[a+neighborSize][b+neighborSize].z+k;
		  //PCA from 8 neighbors to get value at the location - distance weighted
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
		}
	    //From projection in the sample plane, get inplane orientation as graX and graY
	    //Constract matrix for PCA
	    vnl_matrix<float> planePCA(49*2,2,0.0);
	    float tempX, tempY, tempZ; 
	    //x and y axis after rotation (unit vector)
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
	    tempX = 0;
	    tempY = 1;
	    tempZ = 0;
	    yx = cos2*tempX;
	    yy = -sin1*sin2*tempX+cos1*tempY;
	    yz = -cos1*sin2*tempX-sin1*tempY;
	    norm = yx*yx + yy*yy + yz*yz;
	    norm = sqrt(norm);
	    yx = yx/norm;
	    yy = yy/norm;
	    yz = yz/norm;
	    //Get projection X and Y direction
	    t=0; //1d counter
	    for(a=-neighborSize;a<=neighborSize;a++)
	      for(b=-neighborSize;b<=neighborSize;b++)
		{
		  //Projection to orthogonal plane
		  dotProduct = orthoSampleVec[a+neighborSize][b+neighborSize].x * x;
		  dotProduct = dotProduct +  orthoSampleVec[a+neighborSize][b+neighborSize].y * y;
		  dotProduct = dotProduct +  orthoSampleVec[a+neighborSize][b+neighborSize].z * z;
		  kVx = dotProduct * x;
		  kVy = dotProduct * y;
		  kVz = dotProduct * z;
		  kHx = orthoSampleVec[a+neighborSize][b+neighborSize].x - kVx;
		  kHy = orthoSampleVec[a+neighborSize][b+neighborSize].y - kVy;
		  kHz = orthoSampleVec[a+neighborSize][b+neighborSize].z - kVz;
		  //X and Y component
		  //assign sample value
		  planePCA(t,0) = kHx*xx + kHy*xy + kHz*xz;   
		  planePCA(t,1) = kHx*yx + kHy*yy + kHz*yz;
		  planePCA(49+t,0) = -planePCA(t,0);
		  planePCA(49+t,1) = -planePCA(t,1);
		  t = t+1;
		}	    
	    //Analyse the direction in orthogonal plane
	    vnl_matrix<float> SVDplane = vnl_svd<float> (planePCA).V(); 
	    if(SVDplane(0,0)*SVDplane(1,0)>=0)
	      {
		graX = absl(SVDplane(0,0));
		graY = absl(SVDplane(1,0));
	      } 
	    else
	      {
		graX = -absl(SVDplane(0,0));
		graY = absl(SVDplane(1,0)); 		    
	      }
	    //Apply 3D rotation for the orientation vector (graX, graY, 0)
	    //save the normalized result in tangentOri1
	    tangentOri1(i,j,k).x = cos2*graX;
	    tangentOri1(i,j,k).y = -sin1*sin2*graX+cos1*graY;
	    tangentOri1(i,j,k).z = -cos1*sin2*graX-sin1*graY;
	    xt1 = tangentOri1(i,j,k).x;
	    yt1 = tangentOri1(i,j,k).y;
	    zt1 = tangentOri1(i,j,k).z;
	    norm = sqrt(xt1*xt1+yt1*yt1+zt1*zt1);
	    tangentOri1(i,j,k).x = xt1/norm;
	    tangentOri1(i,j,k).y = yt1/norm;
	    tangentOri1(i,j,k).z = zt1/norm;
	  }
    }
  std::cout<<std::endl; 
 
  //Apply PCA based filtering to tangentOri1
  int smoothSize=3; //neighborhood size: 2n+1 * 2n+1 * 2n+1
  VOXEL *smoothTangOri1;
  smoothTangOri1 = (VOXEL *)malloc(sizeof(VOXEL) * pRow * pCol * pSli);
  //First,  record orientation pairs in matrix
  //Second, add additional 0s of the neighborhood size to the matrix
  //Third,  compute SVD for orientation result
  vnl_matrix<float> tangOriPairs((smoothSize*2+1)*(smoothSize*2+1)*(smoothSize*2+1)*2,3,0.0);  
  int neiSize = (smoothSize*2+1)*(smoothSize*2+1)*(smoothSize*2+1);

  std::cout<<"Smoothing the orientation of principal curvature: "<<std::endl;
  for(i=0;i<pRow;i++)
    {
      std::cout<<"\r";
      std::cout<<int((i+1)*100/pRow)<<"%"<<std::flush;
      for(j=0;j<pCol;j++)
	for(k=0;k<pSli;k++)
	  {
	    for(t=0;t<(smoothSize*2+1)*(smoothSize*2+1)*(smoothSize*2+1)*2;t++)
	      {
		tangOriPairs(t,0)=0.0;
		tangOriPairs(t,1)=0.0;
		tangOriPairs(t,2)=0.0;
	      }
	    t=0;
	    //compute gradient magnitude and record in pairs
	    for(a=-smoothSize;a<=smoothSize;a++)
	      for(b=-smoothSize;b<=smoothSize;b++)	  
		for(c=-smoothSize;c<=smoothSize;c++)	  
		  if(((i+a)>=0)&&((i+a)<pRow)&&((j+b)>=0)&&((j+b)<pCol)&&((k+c)>=0)&&((k+c)<pSli))
		    {
		      tangOriPairs(t,0) = tangentOri1(i+a,j+b,k+c).x;
		      tangOriPairs(t,1) = tangentOri1(i+a,j+b,k+c).y;
		      tangOriPairs(t,2) = tangentOri1(i+a,j+b,k+c).z;
		      tangOriPairs(neiSize+t,0) = -tangentOri1(i+a,j+b,k+c).x;
		      tangOriPairs(neiSize+t,1) = -tangentOri1(i+a,j+b,k+c).y;
		      tangOriPairs(neiSize+t,2) = -tangentOri1(i+a,j+b,k+c).z;
		      t = t+1;
		    }
	    //SVD
	    vnl_matrix<float> SVDV = vnl_svd<float> (tangOriPairs).V(); 
	    smoothTangOri1(i,j,k).x = SVDV(0,0);
	    smoothTangOri1(i,j,k).y = SVDV(1,0);
	    smoothTangOri1(i,j,k).z = SVDV(2,0);	  
	  }
    }
  std::cout<<std::endl; 
 
  //Project filtered tangentOri1 back to orthogonal plane
  VOXEL ortho, paral; //decomposition
  float x1,y1,z1;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  x = vecOri(i,j,k).x;
	  y = vecOri(i,j,k).y;
	  z = vecOri(i,j,k).z;
	  x1 = smoothTangOri1(i,j,k).x;
	  y1 = smoothTangOri1(i,j,k).y;
	  z1 = smoothTangOri1(i,j,k).z;	
	  dotProduct = x*x1+y*y1+z*z1;
	  norm = sqrt(x1*x1+y1*y1+z1*z1);
	  ortho.x = x * dotProduct/norm;
	  ortho.y = y * dotProduct/norm;
	  ortho.z = z * dotProduct/norm;
	  paral.x = x1 - ortho.x;
	  paral.y = y1 - ortho.y;
	  paral.z = z1 - ortho.z;
	  norm = sqrt(paral.x*paral.x+paral.y*paral.y+paral.z*paral.z);
	  /*
	  //keep sign
	  if(signf(paral.x)!=signf(tangentOri1(i,j,k).x)) 
	    paral.x = - paral.x;
	  if(signf(paral.y)!=signf(tangentOri1(i,j,k).y)) 
	    paral.y = - paral.y;
	  if(signf(paral.z)!=signf(tangentOri1(i,j,k).z)) 
	    paral.z = - paral.z;
	  //use upper space	  
	  if(paral.z<0)
	    {
	      paral.x = - paral.x;
	      paral.y = - paral.y;
	      paral.z = - paral.z;	      
	    }
	  tangentOri1(i,j,k).x = paral.x/norm;
	  tangentOri1(i,j,k).y = paral.y/norm;
	  tangentOri1(i,j,k).z = paral.z/norm;
	  */
	  //record the normalized cross product in tangentOri2
	  tangentOri2(i,j,k).x = y*tangentOri1(i,j,k).z-z*tangentOri1(i,j,k).y;
	  tangentOri2(i,j,k).y = z*tangentOri1(i,j,k).x-x*tangentOri1(i,j,k).z;
	  tangentOri2(i,j,k).z = x*tangentOri1(i,j,k).y-y*tangentOri1(i,j,k).x;
	  xt2 = tangentOri2(i,j,k).x;
	  yt2 = tangentOri2(i,j,k).y;
	  zt2 = tangentOri2(i,j,k).z;  
	  if(absl(xt2)<0.0001) xt2 = 0;
	  if(absl(yt2)<0.0001) yt2 = 0;
	  if(absl(zt2)<0.0001) zt2 = 0;
	  /*
	  //use upper space 
	  if(zt2<0)
	    {
	      xt2 = -xt2;
	      yt2 = -yt2;
	      zt2 = -zt2;
	    }
	  */	
	  norm = sqrt(xt2*xt2+yt2*yt2+zt2*zt2);
	  tangentOri2(i,j,k).x = xt2/norm;
	  tangentOri2(i,j,k).y = yt2/norm;
	  tangentOri2(i,j,k).z = zt2/norm;
	}
} 
