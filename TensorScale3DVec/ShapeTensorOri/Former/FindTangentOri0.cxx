void FindTangentOri(float *DT, VOXEL *vecOri, VOXEL *tangentOri1, VOXEL *tangentOri2)
{
  int i,j,k;//index
  float x,y,z;//orientation vector
  float xt1,yt1,zt1,xt2,yt2,zt2;//tangent orientation vector
  float norm;//length of the vector
  float sin1,sin2,cos1,cos2;
  int a,b;//3*3 grid index
  VOXEL grid[3][3]; //converted coordinate
  float DTSample[3][3];
  int cornerX, cornerY, cornerZ;
  float weightX, weightY, weightZ;
  int t; //1d counter
  float i1,i2,j1,j2,w1,w2; //temp for trilinear 
  float graX,graY;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{ 
	  x = vecOri(i,j,k).x;
	  y = vecOri(i,j,k).y;
	  z = vecOri(i,j,k).z;
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

	  //if((i==11)&&(j==50)&&(k==74))
	  //std::cout<<x<<" "<<y<<" "<<z<<" "<<std::endl;
	  t = 0;
	  //Get coordinate after rotation
	  for(a=-1;a<=1;a++)
	    for(b=-1;b<=1;b++)
	      {
		grid[a+1][b+1].x = cos2*a;
		grid[a+1][b+1].y = -sin1*sin2*a+cos1*b;
		grid[a+1][b+1].z = -cos1*sin2*a-sin1*b;

		grid[a+1][b+1].x = grid[a+1][b+1].x+i;
		grid[a+1][b+1].y = grid[a+1][b+1].y+j;
		grid[a+1][b+1].z = grid[a+1][b+1].z+k;


		//Trilinear interpolation
		cornerX = floor(grid[a+1][b+1].x);
		cornerY = floor(grid[a+1][b+1].y);
		cornerZ = floor(grid[a+1][b+1].z);
		weightX = grid[a+1][b+1].x - cornerX;
		weightY = grid[a+1][b+1].y - cornerY;
		weightZ = grid[a+1][b+1].z - cornerZ;
		/*
		if(cornerX<0) cornerX=abs(cornerX);
		if(cornerX>pRow-2) cornerX=pRow-1-(cornerX-pRow-2);
		if(cornerY<0) cornerY=abs(cornerY);
		if(cornerY>pCol-2) cornerY=pCol-1-(cornerY-pCol-2);
		if(cornerZ<0) cornerZ=abs(cornerZ);
		if(cornerZ>pSli-2) cornerZ=pSli-1-(cornerZ-pSli-2);
		*/
		if(cornerX<0) cornerX=0;
		if(cornerX>pRow-2) cornerX=pRow-2;
		if(cornerY<0) cornerY=0;
		if(cornerY>pCol-2) cornerY=pCol-2;
		if(cornerZ<0) cornerZ=0;
		if(cornerZ>pSli-2) cornerZ=pSli-2;
		if(cornerX>=0&&(cornerX<pRow-1)&&(cornerY>=0)&&(cornerY<pCol-1)&&(cornerZ>=0)&&(cornerZ<pSli-1))
		  {
		    //3D trilinear interpolation for sample value
		    i1 = DT(cornerX,cornerY,cornerZ)*(1-weightZ) + DT(cornerX,cornerY,cornerZ+1)*weightZ;
		    i2 = DT(cornerX,cornerY+1,cornerZ)*(1-weightZ) + DT(cornerX,cornerY+1,cornerZ+1)*weightZ;
		    j1 = DT(cornerX+1,cornerY,cornerZ)*(1-weightZ) + DT(cornerX+1,cornerY,cornerZ+1)*weightZ;
		    j2 = DT(cornerX+1,cornerY+1,cornerZ)*(1-weightZ) + DT(cornerX+1,cornerY+1,cornerZ+1)*weightZ;
		    w1 = i1*(1-weightY) + i2*weightY;
		    w2 = j1*(1-weightY) + j2*weightY;
		    DTSample[a+1][b+1] = w1*(1-weightX) + w2*weightX;     
		  }
		else
		  DTSample[a+1][b+1] = -10000;
		//if((i==29)&&(j==29)&&(k==76))
		//std::cout<<DTSample[a+1][b+1]<<std::endl;

	      }
	  //Apply sobel operator to sample value
	  graX = DTSample[2][0]+2.0*DTSample[2][1]+DTSample[2][2];
	  graX = graX - (DTSample[0][0]+2.0*DTSample[0][1]+DTSample[0][2]);
	  graX = graX/8.0;
	  graY = DTSample[0][2]+2.0*DTSample[1][2]+DTSample[2][2];	  
	  graY = graY - (DTSample[0][0]+2.0*DTSample[1][0]+DTSample[2][0]);
	  graY = graY/8.0;

	  if(graX<0.00001) graX=0;
	  if(graY<0.00001) graY=0;
	  if((graX==0)&&(graY==0))
	    {
	      graX = 1;
	      graY = 0;
	    }
	  //if((i==11)&&(j==50)&&(k==74))
	  //std::cout<<graX<<" "<<graY<<std::endl;  
	  //if((i==29)&&(j==29)&&(k==76))
	  //std::cout<<graX<<" "<<graY<<std::endl;
	  //Apply 3D rotation for the orientation vector (graX, graY, 0)
	  //save the normalized result in tangentOri1
	  tangentOri1(i,j,k).x = cos2*graX;
	  tangentOri1(i,j,k).y = -sin1*sin2*graX+cos1*graY;
	  tangentOri1(i,j,k).z = -cos1*sin2*graX-sin1*graY;
	  xt1 = tangentOri1(i,j,k).x;
	  yt1 = tangentOri1(i,j,k).y;
	  zt1 = tangentOri1(i,j,k).z;
	  //if((i==11)&&(j==50)&&(k==74))
	  //std::cout<<xt1<<" "<<yt1<<" "<<zt1<<std::endl;
	  norm = sqrt(xt1*xt1+yt1*yt1+zt1*zt1);
	  tangentOri1(i,j,k).x = tangentOri1(i,j,k).x/norm;
	  tangentOri1(i,j,k).y = tangentOri1(i,j,k).y/norm;
	  tangentOri1(i,j,k).z = tangentOri1(i,j,k).z/norm;
	  //record the normalized cross product in tangentOri2
	  tangentOri2(i,j,k).x = y*zt1-z*yt1;
	  tangentOri2(i,j,k).y = z*xt1-x*zt1;
	  tangentOri2(i,j,k).z = x*yt1-y*xt1;
	  xt2 = tangentOri2(i,j,k).x;
	  yt2 = tangentOri2(i,j,k).y;
	  zt2 = tangentOri2(i,j,k).z;
	  norm = sqrt(xt2*xt2+yt2*yt2+zt2*zt2);
	  tangentOri2(i,j,k).x = tangentOri2(i,j,k).x/norm;
	  tangentOri2(i,j,k).y = tangentOri2(i,j,k).y/norm;
	  tangentOri2(i,j,k).z = tangentOri2(i,j,k).z/norm;

	  //if((i==11)&&(j==50)&&(k==74))
	  //std::cout<<xt1<<" "<<yt1<<" "<<zt1<<" "<<xt2<<" "<<yt2<<" "<<zt2<<std::endl;


	}
} 
