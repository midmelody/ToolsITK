//Compute one step filtering result using shape tensor
void shapeTensorFiltering(DIFFUSIZE *diffuSize, float *image, float maxScale, float sigmaNoise)
{   
  int i,j,k,x,y,z;
  const float diffusionConst = 1.0/27.0; //diffusion constant
  float flow, count; //flow value and neighbor number
  float localGradient; //intensity gradient along neighborhood direction
  float conductance; //conductance function value 
  float sigma; //control parameter determining the degree of ﬁltering
  float ellipsePLength,ellipseQLength; //ellipse radius of p/q along pq
  //image after filtering
  float * imageFiltered = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 

  //filtering
  for(i=0;i<pRow;i++)
    {
      std::cout<<"\r";
      std::cout<<int((i+1)*100/pRow)<<"%"<<std::flush;
      for(j=0;j<pCol;j++)
	for(k=0;k<pSli;k++)
	  {
	    flow = 0;
	    count = 1;
	    for(x=-1;x<=1;x++)
	      for(y=-1;y<=1;y++)
		for(z=-1;z<=1;z++)
		  if((i+x)>=0&&((i+x)<pRow-1)&&((j+y)>=0)&&((j+y)<pCol-1)&&((k+z)>=0)&&((k+z)<pSli-1))
		    if((x!=0)||(y!=0)||(z!=0))
		      {
			localGradient = image(i,j,k)-image(i+x,j+y,k+z);
			localGradient = localGradient/sqrt(float(x*x+y*y+z*z));
			if(localGradient<0.00001)
			  localGradient = 1;
			if(((x==0)&&(y==0)&&(z==-1))||((x==0)&&(y==0)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upD;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upD;	    
			  }
			if(((x==0)&&(y==-1)&&(z==-1))||((x==0)&&(y==1)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upN;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upN;	    
			  }
			if(((x==0)&&(y==1)&&(z==-1))||((x==0)&&(y==-1)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upS;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upS;	    
			  }
			if(((x==-1)&&(y==0)&&(z==-1))||((x==1)&&(y==0)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upW;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upW;	    
			  }
			if(((x==1)&&(y==0)&&(z==-1))||((x==-1)&&(y==0)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upE;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upE;	    
			  }
			if(((x==-1)&&(y==-1)&&(z==-1))||((x==1)&&(y==1)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upNW;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upNW;	    
			  }
			if(((x==1)&&(y==-1)&&(z==-1))||((x==-1)&&(y==1)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upNE;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upNE;	    
			  }
			if(((x==-1)&&(y==1)&&(z==-1))||((x==1)&&(y==-1)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upSW;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upSW;	    
			  }
			if(((x==1)&&(y==1)&&(z==-1))||((x==-1)&&(y==-1)&&(z==1)))
			  {
			    ellipsePLength = diffuSize(i,j,k).upSE;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).upSE;	    
			  }
			if(((x==0)&&(y==-1)&&(z==0))||((x==0)&&(y==1)&&(z==0)))
			  {
			    ellipsePLength = diffuSize(i,j,k).cuN;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).cuN;	    
			  }
			if(((x==1)&&(y==0)&&(z==0))||((x==-1)&&(y==0)&&(z==0)))
			  {
			    ellipsePLength = diffuSize(i,j,k).cuE;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).cuE;	    
			  }
			if(((x==-1)&&(y==-1)&&(z==0))||((x==1)&&(y==1)&&(z==0)))
			  {
			    ellipsePLength = diffuSize(i,j,k).cuNW;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).cuNW;	    
			  }
			if(((x==1)&&(y==-1)&&(z==0))||((x==-1)&&(y==1)&&(z==0)))
			  {
			    ellipsePLength = diffuSize(i,j,k).cuNE;
			    ellipseQLength = diffuSize(i+x,j+y,k+z).cuNE;	    
			  }
    
			//sigma controlling flitering degree
			if(ellipsePLength>ellipseQLength)
			  sigma = (1+ellipsePLength)/(1+maxScale)*sigmaNoise;
			else
			  sigma = (1+ellipseQLength)/(1+maxScale)*sigmaNoise;
			
			//conductance
			conductance = exp(-fabs(localGradient*localGradient)/(2*sigma*sigma));

			//add flow and count	
			flow = flow+conductance*localGradient;
			count = count+1;
			//std::cout<<ellipsePLength<<" "<<ellipseQLength<<" "<<sigmaNoise<<" "<<maxScale<<" "<<sigma<<" "<<localGradient<<" "<<flow<<" "<<count<<std::endl;
		      }
	    imageFiltered(i,j,k) = image(i,j,k)-flow/count;
	  }
    }
  //assign value back
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	image(i,j,k) = imageFiltered(i,j,k);
}
