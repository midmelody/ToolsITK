/*
 * =====================================================================================
 *
 *       Filename:  ArrayProcess.h
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  09/15/2011 07:46:31 PM CDT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */


void CreateBlut2( float* arrayBasis0, float* arrayBasis1, float* arrayBasis2, int xgrid, int grid_space_x )
{
    //100 100 100
    double xres=1.0/xgrid;
    //std::cout<<xres<<std::endl;
    double xgd;
    int x;
    int i;
    int index;

    for (x=0; x<xgrid; x++)
    {
        xgd=static_cast<double>(x)*xres;
        
        index=x*4;

        for (i=-1; i<3; i++)
        {
            arrayBasis0[index]=static_cast<float>(dcubic_bsp(xgd-i));
            arrayBasis1[index]=static_cast<float>(dcubic_dbsp(xgd-i));
            arrayBasis2[index]=static_cast<float>(dcubic_d2bsp(xgd-i));

            index++;
        }

    }

}



void CreateLocation4( short int* arrayLocation, float* u_x_existing, float* u_y_existing, float* u_z_existing, int xdim, int ydim, int zdim, int grid_space_x, int grid_space_y, int grid_space_z, int xgrid, int ygrid, int zgrid, unsigned char* mask )
{
    int x, y, z;
    double xe, ye, ze;
    double xBasis, yBasis, zBasis;
    int i0, j0, k0;
    int index;
    int startLocation;
    int maskID;

    double xdist, ydist, zdist;
    int xg, yg, zg;

    int planeSize=xdim*ydim;
    index=0;
    maskID=0;
    for (z=0; z<zdim; z++)
    {
        for (y=0; y<ydim; y++)
        {
            for(x=0; x<xdim; x++)
            {
                index=z*planeSize+y*xdim+x;
                
                if (mask[index]>0)
                {
                    //if (index==562 ||index==618 || index==649 || index==676)
                    //std::cout<<index<<" location: "<<x<<", "<<y<<", "<<z<<std::endl;

                    startLocation=maskID*9;

                    xe=static_cast<double>(x+u_x_existing[index]);
                    ye=static_cast<double>(y+u_y_existing[index]);
                    ze=static_cast<double>(z+u_z_existing[index]);

                    xBasis=xe/grid_space_x;
                    yBasis=ye/grid_space_y;
                    zBasis=ze/grid_space_z;

                    i0=static_cast<int>(xBasis)-1;
                    if (xBasis<0)
                    {
                        i0-=1;
                    }

                    j0=static_cast<int>(yBasis)-1;
                    if (yBasis<0)
                    {
                        j0-=1;
                    }

                    k0=static_cast<int>(zBasis)-1;
                    if (zBasis<0)
                    {
                        k0-=1;
                    }

                    xdist=xBasis-(i0+1);
                    ydist=yBasis-(j0+1);
                    zdist=zBasis-(k0+1);
                    xg=round(xdist*xgrid);
                    if(xg==xgrid)
                    {
                        xg=0;
                        i0++;
                    }
                    yg=round(ydist*ygrid);
                    if(yg==ygrid)
                    {
                        yg=0;
                        j0++;
                    }
                    zg=round(zdist*zgrid);
                    if(zg==zgrid)
                    {
                        zg=0;
                        k0++;
                    }
                    //idg=zg*gridPlaneSize+yg*xgrid+xg;

                    arrayLocation[startLocation]=x;
                    arrayLocation[startLocation+1]=y;
                    arrayLocation[startLocation+2]=z;
                    //arrayLocation[startLocation+3

                    arrayLocation[startLocation+3]=i0;
                    arrayLocation[startLocation+4]=j0;
                    arrayLocation[startLocation+5]=k0;

                    arrayLocation[startLocation+6]=xg;
                    arrayLocation[startLocation+7]=yg;
                    arrayLocation[startLocation+8]=zg;
                    //arrayLocation[startLocation+4]=idg;
                    
                
                    maskID++;
                }

            }
        }
    }


}


void differentiate(float* u, float* ux, float* uy, float* uz, int xdim, int ydim, int zdim)
{
    int planeSize=xdim*ydim;
    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
			    if (x==0)
                {
                    ux[index]=u[index+1]-u[index];//u(x+1,y,z)-u(x,y,z);
                }
                else if (x==xdim-1)
                {
                    ux[index]=u[index]-u[index-1];//u(x,y,z)-u(x-1,y,z);
                }
                else
                {
			        ux[index]=(u[index+1]-u[index-1])/2;//(u(x+1,y,z)-u(x-1,y,z))/2;
                }
			    
                if (y==0)
                {
                    uy[index]=u[index+xdim]-u[index];//u(x,y+1,z)-u(x,y,z);
                }
                else if (y==ydim-1)
                {
                    uy[index]=u[index]-u[index-xdim];//u(x,y,z)-u(x,y-1,z);
                }
                else 
                {
			        uy[index]=(u[index+xdim]-u[index-xdim])/2;//(u(x,y+1,z)-u(x,y-1,z))/2;
                }

                if (z==0)
                {
                    uz[index]=u[index+planeSize]-u[index];//u(x,y,z+1)-u(x,y,z);
                }
                else if (z==zdim-1)
                {
                    uz[index]=u[index]-u[index-planeSize];//(x,y,z)-u(x,y,z-1);
                }
                else
                {
			        uz[index]=(u[index+planeSize]-u[index-planeSize])/2;//(u(x,y,z+1)-u(x,y,z-1))/2;
                }

                index++;
			}
}

void differentiate2(float* u, float* ux, float* uy, float* uz, int xdim, int ydim, int zdim)
{
    int planeSize=xdim*ydim;
    float u_curr;
    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                u_curr=u[index];
			    if (x==0)
                {
                    ux[index]=u[index+1]-u_curr;//u(x+1,y,z)-u(x,y,z);
                }
                else if (x==xdim-1)
                {
                    ux[index]=-u_curr+u[index-1];//-u(x,y,z)+u(x-1,y,z);
                }
                else
                {
			        ux[index]=u[index+1]-2*u_curr+u[index-1];//u(x+1,y,z)-2*u(x,y,z)+u(x-1,y,z);
                }
			    
                if (y==0)
                {
                    uy[index]=u[index+xdim]-u_curr;//u(x,y+1,z)-u(x,y,z);
                }
                else if (y==ydim-1)
                {
                    uy[index]=-u_curr+u[index-xdim];//-u(x,y,z)+u(x,y-1,z);
                }
                else 
                {
			        uy[index]=u[index+xdim]-2*u_curr+u[index-xdim];//u(x,y+1,z)-2*u(x,y,z)+u(x,y-1,z);
                }

                if (z==0)
                {
                    uz[index]=u[index+planeSize]-u_curr;//u(x,y,z+1)-u(x,y,z);
                }
                else if (z==zdim-1)
                {
                    uz[index]=-u[index]+u[index-planeSize];//-u(x,y,z)+u(x,y,z-1);
                }
                else
                {
			        uz[index]=u[index+planeSize]-2*u_curr+u[index-planeSize];//u(x,y,z+1)-2*u(x,y,z)+u(x,y,z-1);
                }

                index++;
			}
}

void differentiate2_crossX(float* u_x, float* tmp, float* u_x_xy, float* u_x_xz, int xdim, int ydim, int zdim)
{
    int planeSize=xdim*ydim;
    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
			    if (x==0)
                {
                    tmp[index]=u_x[index+1]-u_x[index];//u(x+1,y,z)-u(x,y,z);
                }
                else if (x==xdim-1)
                {
                    tmp[index]=u_x[index]-u_x[index-1];//u(x,y,z)-u(x-1,y,z);
                }
                else
                {
			        tmp[index]=(u_x[index+1]-u_x[index-1])/2;//(u(x+1,y,z)-u(x-1,y,z))/2;
                }

                index++;
            }

    index=0;
    for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (y==0)
                {
                    u_x_xy[index]=tmp[index+xdim]-tmp[index];//u(x,y+1,z)-u(x,y,z);
                }
                else if (y==ydim-1)
                {
                    u_x_xy[index]=tmp[index]-tmp[index-xdim];//u(x,y,z)-u(x,y-1,z);
                }
                else 
                {
			        u_x_xy[index]=(tmp[index+xdim]-tmp[index-xdim])/2;//(u(x,y+1,z)-u(x,y-1,z))/2;
                }

                index++;
            }


    index=0;
    for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (z==0)
                {
                    u_x_xz[index]=tmp[index+planeSize]-tmp[index];//u(x,y,z+1)-u(x,y,z);
                }
                else if (z==zdim-1)
                {
                    u_x_xz[index]=tmp[index]-tmp[index-planeSize];//(x,y,z)-u(x,y,z-1);
                }
                else
                {
			        u_x_xz[index]=(tmp[index+planeSize]-tmp[index-planeSize])/2;//(u(x,y,z+1)-u(x,y,z-1))/2;
                }

                index++;
			}


}

                
void differentiate2_crossY(float* u_y, float* tmp, float* u_y_xy, float* u_y_yz, int xdim, int ydim, int zdim)
{
    int planeSize=xdim*ydim;
    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (y==0)
                {
                    tmp[index]=u_y[index+xdim]-u_y[index];//u(x,y+1,z)-u(x,y,z);
                }
                else if (y==ydim-1)
                {
                    tmp[index]=u_y[index]-u_y[index-xdim];//u(x,y,z)-u(x,y-1,z);
                }
                else 
                {
			        tmp[index]=(u_y[index+xdim]-u_y[index-xdim])/2;//(u(x,y+1,z)-u(x,y-1,z))/2;
                }

                index++;
            }


    index=0;
    for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (x==0)
                {
                    u_y_xy[index]=tmp[index+1]-tmp[index];//u(x+1,y,z)-u(x,y,z);
                }
                else if (x==xdim-1)
                {
                    u_y_xy[index]=tmp[index]-tmp[index-1];//u(x,y,z)-u(x-1,y,z);
                }
                else
                {
			        u_y_xy[index]=(tmp[index+1]-tmp[index-1])/2;//(u(x+1,y,z)-u(x-1,y,z))/2;
                }

                index++;
            }


    index=0;
    for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (z==0)
                {
                    u_y_yz[index]=tmp[index+planeSize]-tmp[index];//u(x,y,z+1)-u(x,y,z);
                }
                else if (z==zdim-1)
                {
                    u_y_yz[index]=tmp[index]-tmp[index-planeSize];//(x,y,z)-u(x,y,z-1);
                }
                else
                {
			        u_y_yz[index]=(tmp[index+planeSize]-tmp[index-planeSize])/2;//(u(x,y,z+1)-u(x,y,z-1))/2;
                }

                index++;
			}


}


void differentiate2_crossZ(float* u_z, float* tmp, float* u_z_xz, float* u_z_yz, int xdim, int ydim, int zdim)
{
    int planeSize=xdim*ydim;
    int index=0;
	for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (z==0)
                {
                    tmp[index]=u_z[index+planeSize]-u_z[index];//u(x,y,z+1)-u(x,y,z);
                }
                else if (z==zdim-1)
                {
                    tmp[index]=u_z[index]-u_z[index-planeSize];//(x,y,z)-u(x,y,z-1);
                }
                else
                {
			        tmp[index]=(u_z[index+planeSize]-u_z[index-planeSize])/2;//(u(x,y,z+1)-u(x,y,z-1))/2;
                }

                index++;
            }


    index=0;
    for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (y==0)
                {
                    u_z_xz[index]=tmp[index+xdim]-tmp[index];//u(x,y+1,z)-u(x,y,z);
                }
                else if (y==ydim-1)
                {
                    u_z_xz[index]=tmp[index]-tmp[index-xdim];//u(x,y,z)-u(x,y-1,z);
                }
                else 
                {
			        u_z_xz[index]=(tmp[index+xdim]-tmp[index-xdim])/2;//(u(x,y+1,z)-u(x,y-1,z))/2;
                }

                index++;
            }


    index=0;
    for (int z=0;z<zdim;z++)
		for (int y=0;y<ydim;y++)
			for (int x=0;x<xdim;x++)
			{
                if (x==0)
                {
                    u_z_yz[index]=tmp[index+1]-tmp[index];//u(x+1,y,z)-u(x,y,z);
                }
                else if (x==xdim-1)
                {
                    u_z_yz[index]=tmp[index]-tmp[index-1];//u(x,y,z)-u(x-1,y,z);
                }
                else
                {
			        u_z_yz[index]=(tmp[index+1]-tmp[index-1])/2;//(u(x+1,y,z)-u(x-1,y,z))/2;
                }
                
                index++;
			}


}


void interpolate_uxyz( float* uxx, float* uxy, float* uxz, float* uyx, float* uyy, float* uyz, float* uzx, float* uzy, float* uzz, float* cx, float* cy, float* cz, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int grid_space_x, int grid_space_y, int grid_space_z, float* precalcux, float* precalcuy, float* precalcuz)
{
    int planeSize1=16*grid_space_x*grid_space_y;
    int lineSize1=4*grid_space_x;
    int planeSize2=xdim*ydim;
    int lineSize2=xdim;
    for ( int k=0;k<xdim*ydim*zdim;k++)
    {
        uxx[k]=0;
        uxy[k]=0;
        uxz[k]=0;
        uyx[k]=0;
        uyy[k]=0;
        uyz[k]=0;
        uzx[k]=0;
        uzy[k]=0;
        uzz[k]=0;
    }
    int index=0;
    for (int k=-1;k<num_element_z-1;k++)
    {
        for (int j=-1;j<num_element_y-1;j++)
        {
            for (int i=-1;i<num_element_x-1;i++)
            {
                for (int z=0;z<4*grid_space_z;z++)
                {
                    int curr_z=(k-2)*grid_space_z+z;
                    if (curr_z<0 || curr_z>=zdim)
                        continue;
                    for (int y=0;y<4*grid_space_y;y++)
                    {
                        int curr_y=(j-2)*grid_space_y+y;
                        if (curr_y<0 || curr_y>=ydim)
                            continue;
                        for (int x=0;x<4*grid_space_x;x++)
                        {
                            int curr_x=(i-2)*grid_space_x+x;
                            if (curr_x<0 || curr_x>=xdim)
                                continue;
                            int index1=z*planeSize1+y*lineSize1+x;
                            int index2=curr_z*planeSize2+curr_y*lineSize2+curr_x;
                            uxx[index2]+=cx[index]*precalcux[index1];
                            uxy[index2]+=cx[index]*precalcuy[index1];
                            uxz[index2]+=cx[index]*precalcuz[index1];

                            uyx[index2]+=cy[index]*precalcux[index1];
                            uyy[index2]+=cy[index]*precalcuy[index1];
                            uyz[index2]+=cy[index]*precalcuz[index1];

                            uzx[index2]+=cz[index]*precalcux[index1];
                            uzy[index2]+=cz[index]*precalcuy[index1];
                            uzz[index2]+=cz[index]*precalcuz[index1];
                        }
                    }
                }
                index++;
            }
        }
    }
}

//interpolate disp from coeff
//note: coeff must also be on the image grid!
void interpolate_u( float* ux, float* uy, float* uz, float* cx, float* cy, float* cz, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int grid_space_x, int grid_space_y, int grid_space_z, float* precalcu )
{
    int planeSize1=16*grid_space_x*grid_space_y;
    int lineSize1=4*grid_space_x;
    int planeSize2=xdim*ydim;
    int lineSize2=xdim;
    for ( int k=0;k<xdim*ydim*zdim;k++)
    {
        ux[k]=0;
        uy[k]=0;
        uz[k]=0;
    }
    int index=0;
    for (int k=-1;k<num_element_z-1;k++)
    {
        for (int j=-1;j<num_element_y-1;j++)
        {
            for (int i=-1;i<num_element_x-1;i++)
            {
                //std::cout<<i<<","<<j<<","<<k<<std::endl;
                for (int z=0;z<4*grid_space_z;z++)
                {
                    int curr_z=(k-2)*grid_space_z+z;
                    if (curr_z<0 || curr_z>=zdim)
                        continue;
                    for (int y=0;y<4*grid_space_y;y++)
                    {
                        int curr_y=(j-2)*grid_space_y+y;
                        if (curr_y<0 || curr_y>=ydim)
                            continue;
                        for (int x=0;x<4*grid_space_x;x++)
                        {
                            int curr_x=(i-2)*grid_space_x+x;
                            if (curr_x<0 || curr_x>=xdim)
                                continue;
                            int index1=z*planeSize1+y*lineSize1+x;
                            int index2=curr_z*planeSize2+curr_y*lineSize2+curr_x;
                            ux[index2]+=cx[index]*precalcu[index1];
                            uy[index2]+=cy[index]*precalcu[index1];
                            uz[index2]+=cz[index]*precalcu[index1];
                        }
                    }
                }
                index++;
            }
        }
    }
}


//interpolate disp from coeff
//note: coeff must also be on the image grid!
bool BSplineCoeff2FullDispImage( float* c_x, float* c_y, float* c_z, float* u_x, float* u_y, float* u_z, int grid_space_x, int grid_space_y, int grid_space_z, int scale_x_ratio, int scale_y_ratio, int scale_z_ratio, int num_element_x, int num_element_y, int num_element_z, int xdim, int ydim, int zdim )
{
    int grid_space_x_scaled = grid_space_x*scale_x_ratio;
    int grid_space_y_scaled = grid_space_y*scale_y_ratio;
    int grid_space_z_scaled = grid_space_z*scale_z_ratio;

    int arraySize = 4*4*4*grid_space_x_scaled*grid_space_y_scaled*grid_space_z_scaled;
    float* precalcud = new float[arraySize];

    int planeSizeArray = 16*grid_space_x_scaled*grid_space_y_scaled;
    int lineSizeArray = 4*grid_space_x_scaled;

    for (int z=0;z<4*grid_space_z_scaled;z++)
    {
        for (int y=0;y<4*grid_space_y_scaled;y++)
        {
            for (int x=0;x<4*grid_space_x_scaled;x++)
            {
                int id = z*planeSizeArray+y*lineSizeArray+x;
                precalcud[id]=cubic_bsp(static_cast<float>(x)/grid_space_x_scaled-2)*cubic_bsp(static_cast<float>(y)/grid_space_y_scaled-2)*cubic_bsp(static_cast<float>(z)/grid_space_z_scaled-2);
            }
        }
    }

    int planeSize1=16*grid_space_x_scaled*grid_space_y_scaled;
    int lineSize1=4*grid_space_x_scaled;
    int planeSize2=xdim*ydim;
    int lineSize2=xdim;
    for ( int k=0;k<xdim*ydim*zdim;k++)
    {
        u_x[k]=0;
        u_y[k]=0;
        u_z[k]=0;
    }
    int index=0;
    for (int k=-1;k<num_element_z-1;k++)
    {
        for (int j=-1;j<num_element_y-1;j++)
        {
            for (int i=-1;i<num_element_x-1;i++)
            {
                //std::cout<<i<<","<<j<<","<<k<<std::endl;
                for (int z=0;z<4*grid_space_z_scaled;z++)
                {
                    int curr_z=(k-2)*grid_space_z_scaled+z;
                    if (curr_z<0 || curr_z>=zdim)
                        continue;
                    for (int y=0;y<4*grid_space_y_scaled;y++)
                    {
                        int curr_y=(j-2)*grid_space_y_scaled+y;
                        if (curr_y<0 || curr_y>=ydim)
                            continue;
                        for (int x=0;x<4*grid_space_x_scaled;x++)
                        {
                            int curr_x=(i-2)*grid_space_x_scaled+x;
                            if (curr_x<0 || curr_x>=xdim)
                                continue;
                            int index1=z*planeSize1+y*lineSize1+x;
                            int index2=curr_z*planeSize2+curr_y*lineSize2+curr_x;
                            u_x[index2]+=c_x[index]*precalcud[index1];
                            u_y[index2]+=c_y[index]*precalcud[index1];
                            u_z[index2]+=c_z[index]*precalcud[index1];
                        }
                    }
                }
                index++;
            }
        }
    }

    for ( int k=0;k<xdim*ydim*zdim;k++)
    {
        u_x[k]*=scale_x_ratio;
        u_y[k]*=scale_y_ratio;
        u_z[k]*=scale_z_ratio;
    }

    delete [] precalcud;
    return true;
}


void calcJac( float* Jac, float* u_x, float* u_y, float* u_z, int xdim, int ydim, int zdim )
{
    int imgSize=xdim*ydim*zdim;
    float* u_x_x=new float[imgSize];
    float* u_x_y=new float[imgSize];
    float* u_x_z=new float[imgSize];
    float* u_y_x=new float[imgSize];
    float* u_y_y=new float[imgSize];
    float* u_y_z=new float[imgSize];
    float* u_z_x=new float[imgSize];
    float* u_z_y=new float[imgSize];
    float* u_z_z=new float[imgSize];

    differentiate(u_x, u_x_x, u_x_y, u_x_z, xdim, ydim, zdim);
	differentiate(u_y, u_y_x, u_y_y, u_y_z, xdim, ydim, zdim);
	differentiate(u_z, u_z_x, u_z_y, u_z_z, xdim, ydim, zdim);

    for (int i=0; i<imgSize; i++)
        Jac[i]=(1+u_x_x[i])*(1+u_y_y[i])*(1+u_z_z[i])+u_x_y[i]*u_y_z[i]*u_z_x[i]+u_x_z[i]*u_y_x[i]*u_z_y[i]-u_x_z[i]*(1+u_y_y[i])*u_z_x[i]-u_y_z[i]*u_z_y[i]*(1+u_x_x[i])-u_x_y[i]*u_y_x[i]*(1+u_z_z[i]);

    delete [] u_x_x;
    delete [] u_x_y;
    delete [] u_x_z;
    delete [] u_y_x;
    delete [] u_y_y;
    delete [] u_y_z;
    delete [] u_z_x;
    delete [] u_z_y;
    delete [] u_z_z;

}

void composeJac( float* Jac_total, float* Jac1, float* Jac2_def, int imgSize )
{
    for (int i=0; i<imgSize; i++)
        Jac_total[i]=Jac1[i]*Jac2_def[i];
}


