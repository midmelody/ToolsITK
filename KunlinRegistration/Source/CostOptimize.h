/*
 * =====================================================================================
 *
 *       Filename:  CostOptimize.h
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  12/16/2010 04:00:14 PM CST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////
//define cost and optimization functions

void lap_calc(short int* arrayLocation, float* uxx, float* uxy, float* uxz, float* uyx, float* uyy, float* uyz, float* uzx, float* uzy, float* uzz, int xdim, int ydim, int zdim, float& cost, unsigned char* mask, int maskCT2)
{
    int planeSize=xdim*ydim;
    float lap_x,lap_y,lap_z;
    int index;
    int startLocation;
    int maskID;
    int x,y,z;

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        if (mask[index]>0)
        {
            lap_x=uxx[index]+uxy[index]+uxz[index];
            lap_y=uyx[index]+uyy[index]+uyz[index];
            lap_z=uzx[index]+uzy[index]+uzz[index];
            cost+=lap_x*lap_x+lap_y*lap_y+lap_z*lap_z;
        }

    }

}


//optimization for 2nd order derivative of u
void lap_optimize( ap::real_1d_array& grad, float* arrayBasis0, float* arrayBasis2, short int* arrayLocation, float* uxx, float* uxy, float* uxz, float* uyx, float* uyy, float* uyz, float* uzx, float* uzy, float* uzz, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, float lap, float& cost, unsigned char* mask, int maskCT2, int forrev, int grid_space_x, int grid_space_y, int grid_space_z )
{
    int planeSize=xdim*ydim;
    int coeffPlaneSize=num_element_x*num_element_y;
    int basisSize=num_element_x*num_element_y*num_element_z;
    int bias=forrev*3*basisSize;
    float lapx2=lap*2;
    float lap_x,lap_y,lap_z;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    float temp_basis;
    int idbx,idby,idbz;
    float basisProduct1, basisProduct2, basisProduct3;
    int xg, yg, zg;
    int x, y, z;

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        i0=arrayLocation[startLocation+3];
        j0=arrayLocation[startLocation+4];
        k0=arrayLocation[startLocation+5];

        xg=arrayLocation[startLocation+6];
        yg=arrayLocation[startLocation+7];
        zg=arrayLocation[startLocation+8];

        if (mask[index]>0)
        {
            //change addition order or multiplication order will change the
            //results
            lap_x=uxx[index]+uxy[index]+uxz[index];
            lap_y=uyx[index]+uyy[index]+uyz[index];
            lap_z=uzx[index]+uzy[index]+uzz[index];
            /*float lap_y2=uyy[index]+uyx[index]+uyz[index];
            float lap_z2=uzz[index]+uzx[index]+uzy[index];
            cost+=lap_x*lap_x+lap_y*lap_y+lap_z*lap_z;            

            if (lap_y-lap_y2>0)
                std::cout<<"Y ("<<x<<", "<<y<<", "<<y<<"): "<<uyx[index]<<" + "<<uyy[index]<<" + "<<uyz[index]<<" = "<<lap_y<<", "<<lap_y2<<" diff = "<<lap_y-lap_y2<<std::endl;

            if (lap_z-lap_z2>0)
                std::cout<<"Z ("<<x<<", "<<y<<", "<<z<<"): "<<uzx[index]<<" + "<<uzy[index]<<" + "<<uzz[index]<<" = "<<lap_z<<", "<<lap_z2<<" diff = "<<lap_z-lap_z2<<std::endl;*/

            idbz=zg*4;
            for (k=k0; k<k0+4; k++)
            {
                idby=yg*4;
                for (j=j0; j<j0+4; j++)
                {
                    idbx=xg*4;
                    for (i=i0; i<i0+4; i++)
                    {
                        //std::cout<<"id="<<id<<std::endl;
                        if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                        {
                            basisProduct1=arrayBasis2[idbx]*arrayBasis0[idby]*arrayBasis0[idbz]/(grid_space_x*grid_space_x);
                            basisProduct2=arrayBasis0[idbx]*arrayBasis2[idby]*arrayBasis0[idbz]/(grid_space_y*grid_space_y);
                            basisProduct3=arrayBasis0[idbx]*arrayBasis0[idby]*arrayBasis2[idbz]/(grid_space_z*grid_space_z);

                            temp_basis=basisProduct1+basisProduct2+basisProduct3;

                            //coeff global index (storage form)
                            coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);

                            grad(coeffIndex+1+bias)+=lapx2*lap_x*temp_basis;
                            grad(coeffIndex+basisSize+1+bias)+=lapx2*lap_y*temp_basis;
                            grad(coeffIndex+2*basisSize+1+bias)+=lapx2*lap_z*temp_basis;
                        }
                        idbx++;
                    }
                    idby++;
                }
                idbz++;
            }


        }

    }


}


//optimization for linear elastic smoothness (2nd order derivative of u)
void smooth_clac( short int* arrayLocation, float* uxx, float* uxy, float* uxz, float* uyx, float* uyy, float* uyz, float* uzx, float* uzy, float* uzz, float* ux_xy, float* ux_xz, float* uy_xy, float* uy_yz, float* uz_xz, float* uz_yz, float* ux, float* uy, float* uz, int xdim, int ydim, int zdim, float alpha, float beta, float gamma, float& cost, unsigned char* mask, int maskCT2 )
{
    int planeSize=xdim*ydim;
    float smooth_x,smooth_y,smooth_z;
    int index;
    int startLocation;
    int maskID;
    int x, y, z;

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;


        //has difference, but do not affect result, since smooth_weight = 0
        if (mask[index]>0)
        {
            smooth_x = (alpha+beta)*uxx[index] + alpha*uxy[index] + alpha*uxz[index] 
                        + beta*uy_xy[index] + beta*uz_xz[index] + gamma*ux[index]; 
            smooth_y = (alpha+beta)*uyy[index] + alpha*uyx[index] + alpha*uyz[index]
                        + beta*ux_xy[index] + beta*uz_yz[index] + gamma*uy[index];
            smooth_z = (alpha+beta)*uzz[index] + alpha*uzx[index] + alpha*uzy[index]
                        + beta*ux_xz[index] + beta*uy_yz[index] + gamma*uz[index];

            /*float smooth_y2 = alpha*uyx[index] + (alpha+beta)*uyy[index] + alpha*uyz[index]
                        + beta*ux_xy[index] + beta*uz_yz[index] + gamma*uy[index];
            float smooth_z2 = alpha*uzx[index] + alpha*uzy[index] + (alpha+beta)*uzz[index] 
                        + beta*ux_xz[index] + beta*uy_yz[index] + gamma*uz[index];

            if (smooth_y-smooth_y2>0)
                std::cout<<"Y ("<<x<<", "<<y<<", "<<y<<"): "<<uyx[index]<<" + "<<uyy[index]<<" + "<<uyz[index]<<" = "<<smooth_y<<", "<<smooth_y2<<" diff = "<<smooth_y-smooth_y2<<std::endl;

            if (smooth_z-smooth_z2>0)
                std::cout<<"Z ("<<x<<", "<<y<<", "<<z<<"): "<<uzx[index]<<" + "<<uzy[index]<<" + "<<uzz[index]<<" = "<<smooth_z<<", "<<smooth_z2<<" diff = "<<smooth_z-smooth_z2<<std::endl;
            */


            cost += smooth_x*smooth_x + smooth_y*smooth_y + smooth_z*smooth_z;
                
            //lap_x=uxx[index]+uxy[index]+uxz[index];
            //lap_y=uyx[index]+uyy[index]+uyz[index];
            //lap_z=uzx[index]+uzy[index]+uzz[index];
            //cost+=lap_x*lap_x+lap_y*lap_y+lap_z*lap_z;            

        }

    }


}

//optimization for linear elastic smoothness (2nd order derivative of u)
void smooth_optimize( ap::real_1d_array& grad, float* arrayBasis0, float* arrayBasis1, float* arrayBasis2, short int* arrayLocation, float* uxx, float* uxy, float* uxz, float* uyx, float* uyy, float* uyz, float* uzx, float* uzy, float* uzz, float* ux_xy, float* ux_xz, float* uy_xy, float* uy_yz, float* uz_xz, float* uz_yz, float* ux, float* uy, float* uz, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, float smooth, float alpha, float beta, float gamma, float& cost, unsigned char* mask, int maskCT2, int forrev, int grid_space_x, int grid_space_y, int grid_space_z )
{
    /*smooth=0.01;
    alpha=1.0;
    beta=0.0;
    gamma=0.0;*/

    int planeSize=xdim*ydim;
    int coeffPlaneSize=num_element_x*num_element_y;
    int basisSize=num_element_x*num_element_y*num_element_z;
    int bias=forrev*3*basisSize;
    float smoothx2=smooth*2;
    float smooth_x,smooth_y,smooth_z;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    float temp_basis;
    int idbx,idby,idbz;
    float basisProduct1, basisProduct2, basisProduct3;
    float basisProduct_xy, basisProduct_xz, basisProduct_yz;
    float basisProduct0;
    int xg, yg, zg;
    int x, y, z;

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        i0=arrayLocation[startLocation+3];
        j0=arrayLocation[startLocation+4];
        k0=arrayLocation[startLocation+5];

        xg=arrayLocation[startLocation+6];
        yg=arrayLocation[startLocation+7];
        zg=arrayLocation[startLocation+8];

        if (mask[index]>0)
        {
            smooth_x = (alpha+beta)*uxx[index] + alpha*uxy[index] + alpha*uxz[index] 
                        + beta*uy_xy[index] + beta*uz_xz[index] + gamma*ux[index]; 
            smooth_y = (alpha+beta)*uyy[index] + alpha*uyx[index] + alpha*uyz[index]
                        + beta*ux_xy[index] + beta*uz_yz[index] + gamma*uy[index];
            smooth_z = (alpha+beta)*uzz[index] + alpha*uzx[index] + alpha*uzy[index]
                        + beta*ux_xz[index] + beta*uy_yz[index] + gamma*uz[index];
              
            /*smooth_x = (alpha+beta)*uxx[index] + alpha*uxy[index] + alpha*uxz[index]
                        + beta*uy_xy[index] + beta*uz_xz[index] + gamma*ux[index]; 
            smooth_y = alpha*uyx[index] + (alpha+beta)*uyy[index] + alpha*uyz[index]
                        + beta*ux_xy[index] + beta*uz_yz[index] + gamma*uy[index];
            smooth_z = alpha*uzx[index] + alpha*uzy[index] + (alpha+beta)*uzz[index]
                        + beta*ux_xz[index] + beta*uy_yz[index] + gamma*uz[index];*/

            /*smooth_x=uxx[index]+uxy[index]+uxz[index];
            smooth_y=uyy[index]+uyx[index]+uyz[index];
            smooth_z=uzz[index]+uzx[index]+uzy[index];

            float smooth_y2=uyx[index]+uyy[index]+uyz[index];
            float smooth_z2=uzx[index]+uzy[index]+uzz[index];

            if (smooth_y-smooth_y2>0)
                std::cout<<"Y ("<<x<<", "<<y<<", "<<y<<"): "<<uyx[index]<<" + "<<uyy[index]<<" + "<<uyz[index]<<" = "<<smooth_y<<", "<<smooth_y2<<" diff = "<<smooth_y-smooth_y2<<std::endl;

            if (smooth_z-smooth_z2>0)
                std::cout<<"Z ("<<x<<", "<<y<<", "<<z<<"): "<<uzx[index]<<" + "<<uzy[index]<<" + "<<uzz[index]<<" = "<<smooth_z<<", "<<smooth_z2<<" diff = "<<smooth_z-smooth_z2<<std::endl;
            */


            

            /*smooth_x=uxx[index]+uxy[index]+uxz[index];
            smooth_y=uyx[index]+uyy[index]+uyz[index];
            smooth_z=uzx[index]+uzy[index]+uzz[index];*/


            cost += smooth_x*smooth_x + smooth_y*smooth_y + smooth_z*smooth_z;
                
            //lap_x=uxx[index]+uxy[index]+uxz[index];
            //lap_y=uyx[index]+uyy[index]+uyz[index];
            //lap_z=uzx[index]+uzy[index]+uzz[index];
            //cost+=lap_x*lap_x+lap_y*lap_y+lap_z*lap_z;            

            //index through 64 basis to 1 point
            idbz=zg*4;
            for (k=k0; k<k0+4; k++)
            {
                idby=yg*4;
                for (j=j0; j<j0+4; j++)
                {
                    idbx=xg*4;
                    for (i=i0; i<i0+4; i++)
                    {
                        //std::cout<<"id="<<id<<std::endl;
                        if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                        {

                            basisProduct_xy=arrayBasis1[idbx]*arrayBasis1[idby]*arrayBasis0[idbz]/(grid_space_x*grid_space_y);
                            basisProduct_xz=arrayBasis1[idbx]*arrayBasis0[idby]*arrayBasis1[idbz]/(grid_space_x*grid_space_z);
                            basisProduct_yz=arrayBasis0[idbx]*arrayBasis1[idby]*arrayBasis1[idbz]/(grid_space_y*grid_space_z);

                            basisProduct0=arrayBasis0[idbx]*arrayBasis0[idby]*arrayBasis0[idbz];

                            basisProduct1=arrayBasis2[idbx]*arrayBasis0[idby]*arrayBasis0[idbz]/(grid_space_x*grid_space_x);
                            basisProduct2=arrayBasis0[idbx]*arrayBasis2[idby]*arrayBasis0[idbz]/(grid_space_y*grid_space_y);
                            basisProduct3=arrayBasis0[idbx]*arrayBasis0[idby]*arrayBasis2[idbz]/(grid_space_z*grid_space_z);

                            temp_basis=basisProduct1+basisProduct2+basisProduct3;

                            //coeff global index (storage form)
                            coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);
    
                            grad(coeffIndex+1+bias) += smoothx2*(
                                    smooth_x*( alpha*temp_basis +beta*basisProduct1+gamma*basisProduct0 ) 
                                    + smooth_y*beta*basisProduct_xy 
                                    + smooth_z*beta*basisProduct_xz 
                                    );
                            grad(coeffIndex+basisSize+1+bias) += smoothx2*(
                                    smooth_x*beta*basisProduct_xy
                                    + smooth_y*( alpha*temp_basis +beta*basisProduct2+gamma*basisProduct0 )
                                    + smooth_z*beta*basisProduct_yz
                                    );
                            grad(coeffIndex+2*basisSize+1+bias) += smoothx2*(
                                    smooth_x*beta*basisProduct_xz
                                    + smooth_y*beta*basisProduct_yz 
                                    + smooth_z*( alpha*temp_basis +beta*basisProduct3+gamma*basisProduct0 )
                                    );
                                    
                            /*grad(coeffIndex+1+bias) += smoothx2* smooth_x*temp_basis;
                                    //smooth_x*( alpha*temp_basis )//+beta*basisProduct1+gamma*basisProduct0 ) 
                                    //+ smooth_y*beta*basisProduct_xy 
                                    //+ smooth_z*beta*basisProduct_xz 
                                    //);
                            grad(coeffIndex+basisSize+1+bias) += smoothx2* smooth_y*temp_basis;
                                    //smooth_x*beta*basisProduct_xy
                                    //+ 
                                    //smooth_y*( alpha*temp_basis )//+beta*basisProduct2+gamma*basisProduct0 )
                                    //+ smooth_z*beta*basisProduct_yz
                                    //);
                            grad(coeffIndex+2*basisSize+1+bias) += smoothx2* smooth_z*temp_basis;
                                    //smooth_x*beta*basisProduct_xz
                                    //+ smooth_y*beta*basisProduct_yz 
                                    //+ 
                                    //smooth_z*( alpha*temp_basis )//+beta*basisProduct3+gamma*basisProduct0 )
                                    //);
                                    //
                            */

                            /*grad(coeffIndex+1+bias)+=smoothx2*smooth_x*temp_basis;
                            grad(coeffIndex+basisSize+1+bias)+=smoothx2*smooth_y*temp_basis;
                            grad(coeffIndex+2*basisSize+1+bias)+=smoothx2*smooth_z*temp_basis;*/


                            /*temp_basis=basisProduct1+basisProduct2+basisProduct3;

                            //coeff global index (storage form)
                            coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);

                            grad(coeffIndex+1+bias)+=lapx2*lap_x*temp_basis;
                            grad(coeffIndex+basisSize+1+bias)+=lapx2*lap_y*temp_basis;
                            grad(coeffIndex+2*basisSize+1+bias)+=lapx2*lap_z*temp_basis;
                            */
                        }
                        idbx++;
                    }
                    idby++;
                }
                idbz++;
            }


        }

    }


}



//optimization for linear elastic smoothness (2nd order derivative of u)
void smooth_optimize_debug( ap::real_1d_array& grad, float* arrayBasis0, float* arrayBasis1, float* arrayBasis2, short int* arrayLocation, float* uxx, float* uxy, float* uxz, float* uyx, float* uyy, float*uyz, float* uzx, float* uzy, float* uzz, float* ux_xy, float* ux_xz, float* uy_xy, float* uy_yz, float* uz_xz, float* uz_yz, float* ux, float* uy, float* uz, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, float smooth, float alpha, float beta, float gamma, float& cost, unsigned char* mask, int maskCT2, int forrev, int grid_space_x, int grid_space_y, int grid_space_z )
{
    /*smooth=0.01;
    alpha=1.0;
    beta=0.0;
    gamma=0.0;*/

    int planeSize=xdim*ydim;
    float smooth_x,smooth_y,smooth_z;
    int index; 
    int startLocation;
    int maskID;
    int x, y, z;

    /*int coeffPlaneSize=num_element_x*num_element_y;
    int basisSize=num_element_x*num_element_y*num_element_z;
    int bias=forrev*3*basisSize;
    float smoothx2=smooth*2;
    int coeffIndex;
    int i0, j0, k0;
    int i, j, k;
    float temp_basis;
    int idbx,idby,idbz;
    float basisProduct1, basisProduct2, basisProduct3;
    float basisProduct_xy, basisProduct_xz, basisProduct_yz;
    float basisProduct0;
    int xg, yg, zg;
    */
    

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        if (mask[index]>0)
        {
            /*smooth_x = (alpha+beta)*uxx[index] + alpha*uxy[index] + alpha*uxz[index] 
                        + beta*uy_xy[index] + beta*uz_xz[index] + gamma*ux[index]; 
            smooth_y = (alpha+beta)*uyy[index] + alpha*uyx[index] + alpha*uyz[index]
                        + beta*ux_xy[index] + beta*uz_yz[index] + gamma*uy[index];
            smooth_z = (alpha+beta)*uzz[index] + alpha*uzx[index] + alpha*uzy[index]
                        + beta*ux_xz[index] + beta*uy_yz[index] + gamma*uz[index];*/
              
            smooth_x=uxx[index]+uxy[index]+uxz[index];
            smooth_y=uyy[index]+uyx[index]+uyz[index];
            smooth_z=uzz[index]+uzx[index]+uzy[index];

            float smooth_y2=uyx[index]+uyy[index]+uyz[index];
            float smooth_z2=uzx[index]+uzy[index]+uzz[index];

            if (smooth_y-smooth_y2>0)
                std::cout<<"Y ("<<x<<", "<<y<<", "<<y<<"): "<<uyx[index]<<" + "<<uyy[index]<<" + "<<uyz[index]<<" = "<<smooth_y<<", "<<smooth_y2<<" diff = "<<smooth_y-smooth_y2<<std::endl;

            if (smooth_z-smooth_z2>0)
                std::cout<<"Z ("<<x<<", "<<y<<", "<<z<<"): "<<uzx[index]<<" + "<<uzy[index]<<" + "<<uzz[index]<<" = "<<smooth_z<<", "<<smooth_z2<<" diff = "<<smooth_z-smooth_z2<<std::endl;

            /*smooth_x = (alpha+beta)*uxx[index] + alpha*uxy[index] + alpha*uxz[index]
                        + beta*uy_xy[index] + beta*uz_xz[index] + gamma*ux[index]; 
            smooth_y = alpha*uyx[index] + (alpha+beta)*uyy[index] + alpha*uyz[index]
                        + beta*ux_xy[index] + beta*uz_yz[index] + gamma*uy[index];
            smooth_z = alpha*uzx[index] + alpha*uzy[index] + (alpha+beta)*uzz[index]
                        + beta*ux_xz[index] + beta*uy_yz[index] + gamma*uz[index];*/

            /*smooth_x=uxx[index]+uxy[index]+uxz[index];
            smooth_y=uyx[index]+uyy[index]+uyz[index];
            smooth_z=uzx[index]+uzy[index]+uzz[index];*/


            cost += smooth_x*smooth_x + smooth_y*smooth_y + smooth_z*smooth_z;
                
            //lap_x=uxx[index]+uxy[index]+uxz[index];
            //lap_y=uyx[index]+uyy[index]+uyz[index];
            //lap_z=uzx[index]+uzy[index]+uzz[index];
            //cost+=lap_x*lap_x+lap_y*lap_y+lap_z*lap_z;            

        }

    }


}


void lmk_calc( int xdim, int ydim, int zdim, int num_lmk, int* xc1, int* yc1, int* zc1, int* xc2, int* yc2, int* zc2, float* u_x, float* u_y, float* u_z, float spacing_x, float spacing_y, float spacing_z, float& cost )
{
    //float lmk_dist=0;
    int planeSize=xdim*ydim;
    int l, x, y, z;
    int index;
    float xdeformed, ydeformed, zdeformed;
    float disp_diff_x, disp_diff_y, disp_diff_z;
    
    for (l=0; l<num_lmk; l++)
    {
        x=xc2[l];
        y=yc2[l];
        z=zc2[l];
        index=z*planeSize+y*xdim+x;
        xdeformed=x+u_x[index];
        ydeformed=y+u_y[index];
        zdeformed=z+u_z[index];
        disp_diff_x=(xdeformed-xc1[l])*spacing_x;
        disp_diff_y=(ydeformed-yc1[l])*spacing_y;
        disp_diff_z=(zdeformed-zc1[l])*spacing_z;
        cost+=sqrt(disp_diff_x*disp_diff_x+disp_diff_y*disp_diff_y+disp_diff_z*disp_diff_z);
    }

    cost/=num_lmk;

}


//gradient descent optimization for lmk disp
void lmk_optimize( ap::real_1d_array& grad, float* arrayBasis, short int* arrayLocation, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int num_lmk, int* xc1, int* yc1, int* zc1, int* xc2, int* yc2, int* zc2,  float* u_x, float* u_y, float* u_z, float param, float& cost, int maskCT2, int forrev )
{
    int planeSize=xdim*ydim;
    int coeffPlaneSize=num_element_x*num_element_y;
    int basisSize=num_element_x*num_element_y*num_element_z;
    int bias=forrev*3*basisSize;
    float paramx2=param*2;
    //float lmk_dist=0;
    int l, x, y, z;
    int index, pntIndex, maskID;
    float xdeformed, ydeformed, zdeformed;
    float disp_diff_x, disp_diff_y, disp_diff_z;
    int i0, j0, k0;
    int i, j, k;
    int coeffIndex;
    int startLocation;
    int idbx,idby,idbz;
    float basisProduct;
    int xg, yg, zg;
    int xp, yp, zp;

    for (l=0; l<num_lmk; l++)
    {
        x=xc2[l];
        y=yc2[l];
        z=zc2[l];
        index=z*planeSize+y*xdim+x;
        xdeformed=x+u_x[index];
        ydeformed=y+u_y[index];
        zdeformed=z+u_z[index];
        disp_diff_x=xdeformed-xc1[l];
        disp_diff_y=ydeformed-yc1[l];
        disp_diff_z=zdeformed-zc1[l];
        cost+=(disp_diff_x*disp_diff_x+disp_diff_y*disp_diff_y+disp_diff_z*disp_diff_z);

        for(maskID=0; maskID<maskCT2; maskID++)
        {
            startLocation=maskID*9;

            xp=arrayLocation[startLocation];
            yp=arrayLocation[startLocation+1];
            zp=arrayLocation[startLocation+2];
            pntIndex=zp*planeSize+yp*xdim+xp;

            if (index==pntIndex)
            {
                i0=arrayLocation[startLocation+3];
                j0=arrayLocation[startLocation+4];
                k0=arrayLocation[startLocation+5];

                xg=arrayLocation[startLocation+6];
                yg=arrayLocation[startLocation+7];
                zg=arrayLocation[startLocation+8];

                idbz=zg*4;
                for (k=k0; k<k0+4; k++)
                {
                    idby=yg*4;
                    for (j=j0; j<j0+4; j++)
                    {
                        idbx=xg*4;
                        for (i=i0; i<i0+4; i++)
                        {
                            if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                            {
                                //tmpID=idg*64+id;
                                basisProduct=arrayBasis[idbx]*arrayBasis[idby]*arrayBasis[idbz];

                                //coeff global index (storage form)
                                coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);

                                grad(coeffIndex+1+bias)+=paramx2*disp_diff_x*basisProduct;
                                grad(coeffIndex+basisSize+1+bias)+=paramx2*disp_diff_y*basisProduct;
                                grad(coeffIndex+2*basisSize+1+bias)+=paramx2*disp_diff_z*basisProduct;
                            }

                            idbx++;
                        }
                        idby++;
                    }
                    idbz++;
                }

                break;
            }
        }

    }

    cost/=num_lmk;

}




//calculate cost and gradients---Jx, Jy, Jz
void ved_calc( short int* arrayLocation, float* image2, float* image1, float* u_x_total, float* u_y_total, float* u_z_total, int xdim, int ydim, int zdim, float& cost, unsigned char* mask, int maskCT2 )
{
    int planeSize=xdim*ydim;
    //int imgSize=xdim*ydim*zdim;
    int index;
    int startLocation;
    int maskID;

    float diffVal;
    float image1_defVal;
    float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    int x, y, z;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        if (mask[index]>0)
        {
            image1_defVal=0.0;
            //index2coord(index, x, y, z, xdim, ydim, zdim);
            hx=x+u_x_total[index];
			hy=y+u_y_total[index];
			hz=z+u_z_total[index];

            x1=static_cast<int>(hx);
            if (hx<0)
                x1-=1;
			x2=hx-x1;

			y1=static_cast<int>(hy);
            if (hy<0)
                y1-=1;
			y2=hy-y1;
			
            z1=static_cast<int>(hz);
            if (hz<0)
                z1-=1;
			z2=hz-z1;

            //x3,y3,z3---last point
            x3=x1+1;
            y3=y1+1;
            z3=z1+1;

            if (x1>=0 && x3<xdim && y1>=0 && y3<ydim && z1>=0 && z3<zdim)
            {
                id0=z1*planeSize+y1*xdim+x1;
                id1=id0+1;
                id2=id0+xdim;
                id3=id2+1;
                id4=id0+planeSize;
                id5=id4+1;
                id6=id4+xdim;
                id7=id6+1;

                dist0=(1-x2)*(1-y2)*(1-z2);
                dist1=x2*(1-y2)*(1-z2);
                dist2=(1-x2)*y2*(1-z2);
                dist3=x2*y2*(1-z2);
                dist4=(1-x2)*(1-y2)*z2;
                dist5=x2*(1-y2)*z2;
                dist6=(1-x2)*y2*z2;
                dist7=x2*y2*z2;

				image1_defVal=dist0*image1[id0]+dist1*image1[id1]+dist2*image1[id2]+dist3*image1[id3]+dist4*image1[id4]+dist5*image1[id5]+dist6*image1[id6]+dist7*image1[id7];
            }

            diffVal=image1_defVal-image2[index];
            cost+=diffVal*diffVal;
        }

    }


}



//calculate cost and gradients---Jx, Jy, Jz
void ved_optimize_ratio( ap::real_1d_array& grad, float* arrayBasis, short int* arrayLocation, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, float* image2, float* image1, float* grad_x_T, float* grad_y_T, float* grad_z_T, float* u_x_total, float* u_y_total, float* u_z_total, float ved_weight, float& cost, float& vedcoeff, unsigned char* mask, int maskCT2, int forrev )
{
    int planeSize=xdim*ydim;
    //int imgSize=xdim*ydim*zdim;
    int coeffPlaneSize=num_element_x*num_element_y;
    int basisSize=num_element_x*num_element_y*num_element_z;
    int bias=forrev*3*basisSize;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int idbx, idby, idbz;
    float basisProduct;
    int xg, yg, zg;

    float diffVal;

    float image1_defVal,grad_x_T_defVal,grad_y_T_defVal,grad_z_T_defVal;
    float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    int x, y, z;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        if (mask[index]>0)
        {
            image1_defVal=0.0;
            //index2coord(index, x, y, z, xdim, ydim, zdim);
            hx=x+u_x_total[index];
			hy=y+u_y_total[index];
			hz=z+u_z_total[index];

            x1=static_cast<int>(hx);
            if (hx<0)
                x1-=1;
			x2=hx-x1;

			y1=static_cast<int>(hy);
            if (hy<0)
                y1-=1;
			y2=hy-y1;
			
            z1=static_cast<int>(hz);
            if (hz<0)
                z1-=1;
			z2=hz-z1;

            //x3,y3,z3---last point
            x3=x1+1;
            y3=y1+1;
            z3=z1+1;

            if (x1>=0 && x3<xdim && y1>=0 && y3<ydim && z1>=0 && z3<zdim)
            {
                id0=z1*planeSize+y1*xdim+x1;
                id1=id0+1;
                id2=id0+xdim;
                id3=id2+1;
                id4=id0+planeSize;
                id5=id4+1;
                id6=id4+xdim;
                id7=id6+1;

                dist0=(1-x2)*(1-y2)*(1-z2);
                dist1=x2*(1-y2)*(1-z2);
                dist2=(1-x2)*y2*(1-z2);
                dist3=x2*y2*(1-z2);
                dist4=(1-x2)*(1-y2)*z2;
                dist5=x2*(1-y2)*z2;
                dist6=(1-x2)*y2*z2;
                dist7=x2*y2*z2;

				image1_defVal=dist0*image1[id0]+dist1*image1[id1]+dist2*image1[id2]+dist3*image1[id3]+dist4*image1[id4]+dist5*image1[id5]+dist6*image1[id6]+dist7*image1[id7];
            }

            diffVal=image1_defVal-image2[index];
            cost+=diffVal*diffVal;
        }
    }

    vedcoeff=ved_weight/cost;
    float vedx2=vedcoeff*2;
    

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        i0=arrayLocation[startLocation+3];
        j0=arrayLocation[startLocation+4];
        k0=arrayLocation[startLocation+5];

        xg=arrayLocation[startLocation+6];
        yg=arrayLocation[startLocation+7];
        zg=arrayLocation[startLocation+8];

        if (mask[index]>0)
        {
            image1_defVal=0.0;
            grad_x_T_defVal=0.0;
            grad_y_T_defVal=0.0;
            grad_z_T_defVal=0.0;

            //index2coord(index, x, y, z, xdim, ydim, zdim);
            hx=x+u_x_total[index];
			hy=y+u_y_total[index];
			hz=z+u_z_total[index];

            x1=static_cast<int>(hx);
            if (hx<0)
                x1-=1;
			x2=hx-x1;

			y1=static_cast<int>(hy);
            if (hy<0)
                y1-=1;
			y2=hy-y1;
			
            z1=static_cast<int>(hz);
            if (hz<0)
                z1-=1;
			z2=hz-z1;

            //x3,y3,z3---last point
            x3=x1+1;
            y3=y1+1;
            z3=z1+1;

            if (x1>=0 && x3<xdim && y1>=0 && y3<ydim && z1>=0 && z3<zdim)
            {
                id0=z1*planeSize+y1*xdim+x1;
                id1=id0+1;
                id2=id0+xdim;
                id3=id2+1;
                id4=id0+planeSize;
                id5=id4+1;
                id6=id4+xdim;
                id7=id6+1;

                dist0=(1-x2)*(1-y2)*(1-z2);
                dist1=x2*(1-y2)*(1-z2);
                dist2=(1-x2)*y2*(1-z2);
                dist3=x2*y2*(1-z2);
                dist4=(1-x2)*(1-y2)*z2;
                dist5=x2*(1-y2)*z2;
                dist6=(1-x2)*y2*z2;
                dist7=x2*y2*z2;

				image1_defVal=dist0*image1[id0]+dist1*image1[id1]+dist2*image1[id2]+dist3*image1[id3]+dist4*image1[id4]+dist5*image1[id5]+dist6*image1[id6]+dist7*image1[id7];

				grad_x_T_defVal=dist0*grad_x_T[id0]+dist1*grad_x_T[id1]+dist2*grad_x_T[id2]+dist3*grad_x_T[id3]+dist4*grad_x_T[id4]+dist5*grad_x_T[id5]+dist6*grad_x_T[id6]+dist7*grad_x_T[id7];
				grad_y_T_defVal=dist0*grad_y_T[id0]+dist1*grad_y_T[id1]+dist2*grad_y_T[id2]+dist3*grad_y_T[id3]+dist4*grad_y_T[id4]+dist5*grad_y_T[id5]+dist6*grad_y_T[id6]+dist7*grad_y_T[id7];
				grad_z_T_defVal=dist0*grad_z_T[id0]+dist1*grad_z_T[id1]+dist2*grad_z_T[id2]+dist3*grad_z_T[id3]+dist4*grad_z_T[id4]+dist5*grad_z_T[id5]+dist6*grad_z_T[id6]+dist7*grad_z_T[id7];
            }

            diffVal=image1_defVal-image2[index];

            idbz=zg*4;
            for (k=k0; k<k0+4; k++)
            {
                idby=yg*4;
                for (j=j0; j<j0+4; j++)
                {
                    idbx=xg*4;
                    for (i=i0; i<i0+4; i++)
                    {
                        if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                        {
                            basisProduct=arrayBasis[idbx]*arrayBasis[idby]*arrayBasis[idbz];
                            //coeff global index (storage form)
                            coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);

                            grad(coeffIndex+1+bias)+=vedx2*diffVal*grad_x_T_defVal*basisProduct;
                            grad(coeffIndex+basisSize+1+bias)+=vedx2*diffVal*grad_y_T_defVal*basisProduct;
                            grad(coeffIndex+2*basisSize+1+bias)+=vedx2*diffVal*grad_z_T_defVal*basisProduct;
                        }
                        idbx++;
                    }
                    idby++;
                }
                idbz++;
            }

        }

    }


}


//calculate cost and gradients---Jx, Jy, Jz
void ved_optimize( ap::real_1d_array& grad, float* arrayBasis, short int* arrayLocation, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, float* image2, float* image1, float* grad_x_T, float* grad_y_T, float* grad_z_T, float* u_x_total, float* u_y_total, float* u_z_total, float ved, float& cost, unsigned char* mask, int maskCT2, int forrev )
{
    int planeSize=xdim*ydim;
    //int imgSize=xdim*ydim*zdim;
    int coeffPlaneSize=num_element_x*num_element_y;
    int basisSize=num_element_x*num_element_y*num_element_z;
    int bias=forrev*3*basisSize;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int idbx, idby, idbz;
    float basisProduct;
    int xg, yg, zg;

    float diffVal;

    float image1_defVal,grad_x_T_defVal,grad_y_T_defVal,grad_z_T_defVal;
    float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    int x, y, z;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;


    float vedx2=ved*2;
    

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        i0=arrayLocation[startLocation+3];
        j0=arrayLocation[startLocation+4];
        k0=arrayLocation[startLocation+5];

        xg=arrayLocation[startLocation+6];
        yg=arrayLocation[startLocation+7];
        zg=arrayLocation[startLocation+8];

        if (mask[index]>0)
        {
            image1_defVal=0.0;
            grad_x_T_defVal=0.0;
            grad_y_T_defVal=0.0;
            grad_z_T_defVal=0.0;

            //index2coord(index, x, y, z, xdim, ydim, zdim);
            hx=x+u_x_total[index];
			hy=y+u_y_total[index];
			hz=z+u_z_total[index];

            x1=static_cast<int>(hx);
            if (hx<0)
                x1-=1;
			x2=hx-x1;

			y1=static_cast<int>(hy);
            if (hy<0)
                y1-=1;
			y2=hy-y1;
			
            z1=static_cast<int>(hz);
            if (hz<0)
                z1-=1;
			z2=hz-z1;

            //x3,y3,z3---last point
            x3=x1+1;
            y3=y1+1;
            z3=z1+1;

            if (x1>=0 && x3<xdim && y1>=0 && y3<ydim && z1>=0 && z3<zdim)
            {
                id0=z1*planeSize+y1*xdim+x1;
                id1=id0+1;
                id2=id0+xdim;
                id3=id2+1;
                id4=id0+planeSize;
                id5=id4+1;
                id6=id4+xdim;
                id7=id6+1;

                dist0=(1-x2)*(1-y2)*(1-z2);
                dist1=x2*(1-y2)*(1-z2);
                dist2=(1-x2)*y2*(1-z2);
                dist3=x2*y2*(1-z2);
                dist4=(1-x2)*(1-y2)*z2;
                dist5=x2*(1-y2)*z2;
                dist6=(1-x2)*y2*z2;
                dist7=x2*y2*z2;

				image1_defVal=dist0*image1[id0]+dist1*image1[id1]+dist2*image1[id2]+dist3*image1[id3]+dist4*image1[id4]+dist5*image1[id5]+dist6*image1[id6]+dist7*image1[id7];

				grad_x_T_defVal=dist0*grad_x_T[id0]+dist1*grad_x_T[id1]+dist2*grad_x_T[id2]+dist3*grad_x_T[id3]+dist4*grad_x_T[id4]+dist5*grad_x_T[id5]+dist6*grad_x_T[id6]+dist7*grad_x_T[id7];
				grad_y_T_defVal=dist0*grad_y_T[id0]+dist1*grad_y_T[id1]+dist2*grad_y_T[id2]+dist3*grad_y_T[id3]+dist4*grad_y_T[id4]+dist5*grad_y_T[id5]+dist6*grad_y_T[id6]+dist7*grad_y_T[id7];
				grad_z_T_defVal=dist0*grad_z_T[id0]+dist1*grad_z_T[id1]+dist2*grad_z_T[id2]+dist3*grad_z_T[id3]+dist4*grad_z_T[id4]+dist5*grad_z_T[id5]+dist6*grad_z_T[id6]+dist7*grad_z_T[id7];
            }

            diffVal=image1_defVal-image2[index];
            cost+=diffVal*diffVal;

            idbz=zg*4;
            for (k=k0; k<k0+4; k++)
            {
                idby=yg*4;
                for (j=j0; j<j0+4; j++)
                {
                    idbx=xg*4;
                    for (i=i0; i<i0+4; i++)
                    {
                        if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                        {
                            basisProduct=arrayBasis[idbx]*arrayBasis[idby]*arrayBasis[idbz];
                            //coeff global index (storage form)
                            coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);

                            grad(coeffIndex+1+bias)+=vedx2*diffVal*grad_x_T_defVal*basisProduct;
                            grad(coeffIndex+basisSize+1+bias)+=vedx2*diffVal*grad_y_T_defVal*basisProduct;
                            grad(coeffIndex+2*basisSize+1+bias)+=vedx2*diffVal*grad_z_T_defVal*basisProduct;
                        }
                        idbx++;
                    }
                    idby++;
                }
                idbz++;
            }

        }

    }


}

//calculate cost and gradients---Jx, Jy, Jz
void ssd_optimize( ap::real_1d_array& grad, float* arrayBasis, short int* arrayLocation, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, float* image2, float* image1, float* grad_x_T, float* grad_y_T, float* grad_z_T, float* u_x_total, float* u_y_total, float* u_z_total, float ssd, float& cost, unsigned char* mask, int maskCT2, int forrev )
{
    int planeSize=xdim*ydim;
    //int imgSize=xdim*ydim*zdim;
    int coeffPlaneSize=num_element_x*num_element_y;
    int basisSize=num_element_x*num_element_y*num_element_z;
    int bias=forrev*3*basisSize;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int idbx, idby, idbz;
    float basisProduct;
    int xg, yg, zg;

    float diffVal;

    float image1_defVal,grad_x_T_defVal,grad_y_T_defVal,grad_z_T_defVal;
    float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    int x, y, z;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;


    float ssdx2=ssd*2;
    

    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        i0=arrayLocation[startLocation+3];
        j0=arrayLocation[startLocation+4];
        k0=arrayLocation[startLocation+5];

        xg=arrayLocation[startLocation+6];
        yg=arrayLocation[startLocation+7];
        zg=arrayLocation[startLocation+8];

        if (mask[index]>0)
        {
            image1_defVal=0.0;
            grad_x_T_defVal=0.0;
            grad_y_T_defVal=0.0;
            grad_z_T_defVal=0.0;

            //index2coord(index, x, y, z, xdim, ydim, zdim);
            hx=x+u_x_total[index];
			hy=y+u_y_total[index];
			hz=z+u_z_total[index];

            x1=static_cast<int>(hx);
            if (hx<0)
                x1-=1;
			x2=hx-x1;

			y1=static_cast<int>(hy);
            if (hy<0)
                y1-=1;
			y2=hy-y1;
			
            z1=static_cast<int>(hz);
            if (hz<0)
                z1-=1;
			z2=hz-z1;

            //x3,y3,z3---last point
            x3=x1+1;
            y3=y1+1;
            z3=z1+1;

            if (x1>=0 && x3<xdim && y1>=0 && y3<ydim && z1>=0 && z3<zdim)
            {
                id0=z1*planeSize+y1*xdim+x1;
                id1=id0+1;
                id2=id0+xdim;
                id3=id2+1;
                id4=id0+planeSize;
                id5=id4+1;
                id6=id4+xdim;
                id7=id6+1;

                dist0=(1-x2)*(1-y2)*(1-z2);
                dist1=x2*(1-y2)*(1-z2);
                dist2=(1-x2)*y2*(1-z2);
                dist3=x2*y2*(1-z2);
                dist4=(1-x2)*(1-y2)*z2;
                dist5=x2*(1-y2)*z2;
                dist6=(1-x2)*y2*z2;
                dist7=x2*y2*z2;

				image1_defVal=dist0*image1[id0]+dist1*image1[id1]+dist2*image1[id2]+dist3*image1[id3]+dist4*image1[id4]+dist5*image1[id5]+dist6*image1[id6]+dist7*image1[id7];

				grad_x_T_defVal=dist0*grad_x_T[id0]+dist1*grad_x_T[id1]+dist2*grad_x_T[id2]+dist3*grad_x_T[id3]+dist4*grad_x_T[id4]+dist5*grad_x_T[id5]+dist6*grad_x_T[id6]+dist7*grad_x_T[id7];
				grad_y_T_defVal=dist0*grad_y_T[id0]+dist1*grad_y_T[id1]+dist2*grad_y_T[id2]+dist3*grad_y_T[id3]+dist4*grad_y_T[id4]+dist5*grad_y_T[id5]+dist6*grad_y_T[id6]+dist7*grad_y_T[id7];
				grad_z_T_defVal=dist0*grad_z_T[id0]+dist1*grad_z_T[id1]+dist2*grad_z_T[id2]+dist3*grad_z_T[id3]+dist4*grad_z_T[id4]+dist5*grad_z_T[id5]+dist6*grad_z_T[id6]+dist7*grad_z_T[id7];
            }

            diffVal=image1_defVal-image2[index];
            cost+=diffVal*diffVal;

            idbz=zg*4;
            for (k=k0; k<k0+4; k++)
            {
                idby=yg*4;
                for (j=j0; j<j0+4; j++)
                {
                    idbx=xg*4;
                    for (i=i0; i<i0+4; i++)
                    {
                        if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                        {
                            basisProduct=arrayBasis[idbx]*arrayBasis[idby]*arrayBasis[idbz];
                            //coeff global index (storage form)
                            coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);

                            grad(coeffIndex+1+bias)+=ssdx2*diffVal*grad_x_T_defVal*basisProduct;
                            grad(coeffIndex+basisSize+1+bias)+=ssdx2*diffVal*grad_y_T_defVal*basisProduct;
                            grad(coeffIndex+2*basisSize+1+bias)+=ssdx2*diffVal*grad_z_T_defVal*basisProduct;
                        }
                        idbx++;
                    }
                    idby++;
                }
                idbz++;
            }

        }

    }


}

//calculate cost and gradients---Jx, Jy, Jz
void sstvd_optimize( ap::real_1d_array& grad, float* arrayBasis0, float* arrayBasis1, short int* arrayLocation, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, float* image2, float* image1, float* grad_x_T, float* grad_y_T, float* grad_z_T, float* u_x_total, float* u_y_total, float* u_z_total, float* u_x_x, float* u_x_y, float* u_x_z, float* u_y_x, float* u_y_y, float* u_y_z, float* u_z_x, float* u_z_y, float* u_z_z, float sstvd, float& cost, unsigned char* mask, int maskCT2, int forrev, float* Je, int grid_space_x, int grid_space_y, int grid_space_z, float& JacMin, float& JacMax )
{
    //int imgSize=xdim*ydim*zdim;
    int planeSize=xdim*ydim;
    int coeffPlaneSize=num_element_x*num_element_y;
    int basisSize=num_element_x*num_element_y*num_element_z;
    int bias=forrev*3*basisSize;
    float Jx,Jy,Jz;
    float sstvdx2=sstvd*2;
    int i0, j0, k0;
    int i, j, k;
    int index, coeffIndex;
    int startLocation;
    int maskID;
    int idbx,idby,idbz;
    float basisProduct0, basisProduct1, basisProduct2, basisProduct3;
    int xg, yg, zg;

    float JcVal;
    float JacVal;
    float diffVal;

    float image1_defVal,grad_x_T_defVal,grad_y_T_defVal,grad_z_T_defVal;
    float dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
    int id0, id1, id2, id3, id4, id5, id6, id7;
    int x, y, z;
    float x2, y2, z2;
    float hx, hy, hz;
    int x1, y1, z1;
    int x3, y3, z3;

    /*int id;
    float zero=0.0;
    float image111,image112,image121,image122,image211,image212,image221,image222;
    float gx111,gx112,gx121,gx122,gx211,gx212,gx221,gx222;
    float gy111,gy112,gy121,gy122,gy211,gy212,gy221,gy222;
    float gz111,gz112,gz121,gz122,gz211,gz212,gz221,gz222;*/
    for (maskID=0; maskID<maskCT2; maskID++)
    {
        startLocation=maskID*9;

        x=arrayLocation[startLocation];
        y=arrayLocation[startLocation+1];
        z=arrayLocation[startLocation+2];
        index=z*planeSize+y*xdim+x;

        i0=arrayLocation[startLocation+3];
        j0=arrayLocation[startLocation+4];
        k0=arrayLocation[startLocation+5];

        xg=arrayLocation[startLocation+6];
        yg=arrayLocation[startLocation+7];
        zg=arrayLocation[startLocation+8];

        if (mask[index]>0)
        {
            image1_defVal=0.0;
            grad_x_T_defVal=0.0;
            grad_y_T_defVal=0.0;
            grad_z_T_defVal=0.0;

            //index2coord(index, x, y, z, xdim, ydim, zdim);
            hx=x+u_x_total[index];
			hy=y+u_y_total[index];
			hz=z+u_z_total[index];

            //x1,y1,z1---starting point; 
            //x2,y2,z2---dist from starting point
            x1=static_cast<int>(hx);
            if (hx<0)
                x1-=1;
			x2=hx-x1;

			y1=static_cast<int>(hy);
            if (hy<0)
                y1-=1;
			y2=hy-y1;
			
            z1=static_cast<int>(hz);
            if (hz<0)
                z1-=1;
			z2=hz-z1;

            //x3,y3,z3---last point
            x3=x1+1;
            y3=y1+1;
            z3=z1+1;

            /*
            image111=zero;
			image211=zero;
			image121=zero;
			image221=zero;
			image112=zero;
			image212=zero;
			image122=zero;
			image222=zero;

            gx111=zero;
			gx211=zero;
			gx121=zero;
			gx221=zero;
			gx112=zero;
			gx212=zero;
			gx122=zero;
			gx222=zero;

            gy111=zero;
			gy211=zero;
			gy121=zero;
			gy221=zero;
			gy112=zero;
			gy212=zero;
			gy122=zero;
			gy222=zero;

            gz111=zero;
			gz211=zero;
			gz121=zero;
			gz221=zero;
			gz112=zero;
			gz212=zero;
			gz122=zero;
			gz222=zero;

            //linear interpolating
            id=z1*planeSize+y1*xdim+x1;
            if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
            {
			    image111=image1[id];//(x1,y1,z1);
			    gx111=grad_x_T[id];
			    gy111=grad_y_T[id];
			    gz111=grad_z_T[id];
            }
            if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z1>=0 && z1<zdim)
            {
			    image211=image1[id+1];//(x3,y1,z1);
			    gx211=grad_x_T[id+1];
			    gy211=grad_y_T[id+1];
			    gz211=grad_z_T[id+1];
            }
            if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
            {
			    image121=image1[id+xdim];//(x1,y3,z1);
			    gx121=grad_x_T[id+xdim];
			    gy121=grad_y_T[id+xdim];
			    gz121=grad_z_T[id+xdim];
            }
            if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z1>=0 && z1<zdim)
            {
			    image221=image1[id+xdim+1];//(x3,y3,z1);
			    gx221=grad_x_T[id+xdim+1];
			    gy221=grad_y_T[id+xdim+1];
			    gz221=grad_z_T[id+xdim+1];
            }
            if (x1>=0 && x1<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
            {
			    image112=image1[id+planeSize];//(x1,y1,z3);
			    gx112=grad_x_T[id+planeSize];
			    gy112=grad_y_T[id+planeSize];
			    gz112=grad_z_T[id+planeSize];
            }
            if (x3>=0 && x3<xdim && y1>=0 && y1<ydim && z3>=0 && z3<zdim)
            {
			    image212=image1[id+planeSize+1];//(x3,y1,z3);
			    gx212=grad_x_T[id+planeSize+1];
			    gy212=grad_y_T[id+planeSize+1];
			    gz212=grad_z_T[id+planeSize+1];
            }
            if (x1>=0 && x1<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
            {
			    image122=image1[id+planeSize+xdim];//(x1,y3,z3);
			    gx122=grad_x_T[id+planeSize+xdim];
			    gy122=grad_y_T[id+planeSize+xdim];
			    gz122=grad_z_T[id+planeSize+xdim];
            }
            if (x3>=0 && x3<xdim && y3>=0 && y3<ydim && z3>=0 && z3<zdim)
            {
			    image222=image1[id+planeSize+xdim+1];//(x3,y3,z3);
			    gx222=grad_x_T[id+planeSize+xdim+1];
			    gy222=grad_y_T[id+planeSize+xdim+1];
			    gz222=grad_z_T[id+planeSize+xdim+1];
            }

			image1_defVal=(1-x2)*(1-y2)*(1-z2)*image111+x2*(1-y2)*(1-z2)*image211+(1-x2)*y2*(1-z2)*image121+x2*y2*(1-z2)*image221+(1-x2)*(1-y2)*z2*image112+x2*(1-y2)*z2*image212+(1-x2)*y2*z2*image122+x2*y2*z2*image222;
			grad_x_T_defVal=(1-x2)*(1-y2)*(1-z2)*gx111+x2*(1-y2)*(1-z2)*gx211+(1-x2)*y2*(1-z2)*gx121+x2*y2*(1-z2)*gx221+(1-x2)*(1-y2)*z2*gx112+x2*(1-y2)*z2*gx212+(1-x2)*y2*z2*gx122+x2*y2*z2*gx222;
			grad_y_T_defVal=(1-x2)*(1-y2)*(1-z2)*gy111+x2*(1-y2)*(1-z2)*gy211+(1-x2)*y2*(1-z2)*gy121+x2*y2*(1-z2)*gy221+(1-x2)*(1-y2)*z2*gy112+x2*(1-y2)*z2*gy212+(1-x2)*y2*z2*gy122+x2*y2*z2*gy222;
			grad_z_T_defVal=(1-x2)*(1-y2)*(1-z2)*gz111+x2*(1-y2)*(1-z2)*gz211+(1-x2)*y2*(1-z2)*gz121+x2*y2*(1-z2)*gz221+(1-x2)*(1-y2)*z2*gz112+x2*(1-y2)*z2*gz212+(1-x2)*y2*z2*gz122+x2*y2*z2*gz222;
            */


            
            if (x1>=0 && x3<xdim && y1>=0 && y3<ydim && z1>=0 && z3<zdim)
            {
                id0=z1*planeSize+y1*xdim+x1;
                id1=id0+1;
                id2=id0+xdim;
                id3=id2+1;
                id4=id0+planeSize;
                id5=id4+1;
                id6=id4+xdim;
                id7=id6+1;

                dist0=(1-x2)*(1-y2)*(1-z2);
                dist1=x2*(1-y2)*(1-z2);
                dist2=(1-x2)*y2*(1-z2);
                dist3=x2*y2*(1-z2);
                dist4=(1-x2)*(1-y2)*z2;
                dist5=x2*(1-y2)*z2;
                dist6=(1-x2)*y2*z2;
                dist7=x2*y2*z2;

				image1_defVal=dist0*image1[id0]+dist1*image1[id1]+dist2*image1[id2]+dist3*image1[id3]+dist4*image1[id4]+dist5*image1[id5]+dist6*image1[id6]+dist7*image1[id7];

				grad_x_T_defVal=dist0*grad_x_T[id0]+dist1*grad_x_T[id1]+dist2*grad_x_T[id2]+dist3*grad_x_T[id3]+dist4*grad_x_T[id4]+dist5*grad_x_T[id5]+dist6*grad_x_T[id6]+dist7*grad_x_T[id7];
				grad_y_T_defVal=dist0*grad_y_T[id0]+dist1*grad_y_T[id1]+dist2*grad_y_T[id2]+dist3*grad_y_T[id3]+dist4*grad_y_T[id4]+dist5*grad_y_T[id5]+dist6*grad_y_T[id6]+dist7*grad_y_T[id7];
				grad_z_T_defVal=dist0*grad_z_T[id0]+dist1*grad_z_T[id1]+dist2*grad_z_T[id2]+dist3*grad_z_T[id3]+dist4*grad_z_T[id4]+dist5*grad_z_T[id5]+dist6*grad_z_T[id6]+dist7*grad_z_T[id7];
            }   
            


            JcVal=(1+u_x_x[index])*(1+u_y_y[index])*(1+u_z_z[index])+u_x_y[index]*u_y_z[index]*u_z_x[index]+u_x_z[index]*u_y_x[index]*u_z_y[index]-u_x_z[index]*(1+u_y_y[index])*u_z_x[index]-u_y_z[index]*u_z_y[index]*(1+u_x_x[index])-u_x_y[index]*u_y_x[index]*(1+u_z_z[index]);
            JacVal=Je[index]*JcVal;
            if (JacVal>JacMax)
                JacMax=JacVal;
            if (JacVal<JacMin)
                JacMin=JacVal;

            diffVal=(JacVal*(image1_defVal+1000.0)-(image2[index]+1000.0))/1055.0;
            cost+=diffVal*diffVal;
            //if (temp_cost != temp_cost)
            /*if(index==618)
            //if (temp_cost>10)
            {
                //std::cout<<"if in mask: "<<static_cast<int>(mask[index])<<std::endl;
                //std::cout<<index<<": "<<x<<", "<<y<<", "<<z<<std::endl;
                //std::cout<<x+y*xdim+z*xdim*ydim<<std::endl;
                std::cout<<index<<": "<<image1_defVal<<" "<<temp_cost<<" "<<Jac[index]<<std::endl;
                std::cout<<"derivatives: "<<grad_x_T_defVal<<" "<<grad_y_T_defVal<<" "<<grad_z_T_defVal<<std::endl;
                //std::cout<<u_x_x[index]<<", "<<u_x_y[index]<<", "<<u_x_z[index]<<std::endl;
                //std::cout<<u_y_x[index]<<", "<<u_y_y[index]<<", "<<u_y_z[index]<<std::endl;
                //std::cout<<u_z_x[index]<<", "<<u_z_y[index]<<", "<<u_z_z[index]<<std::endl;
            }*/
        

            idbz=zg*4;
            for (k=k0; k<k0+4; k++)
            {
                idby=yg*4;
                for (j=j0; j<j0+4; j++)
                {
                    idbx=xg*4;
                    for (i=i0; i<i0+4; i++)
                    {
                        //std::cout<<"id="<<id<<std::endl;
                        if (i>=-1 && i<num_element_x-1 && j>=-1 && j<num_element_y-1 && k>=-1 && k<num_element_z-1)
                        {
                            //tmpID=idg*64+id;
                            basisProduct0=arrayBasis0[idbx]*arrayBasis0[idby]*arrayBasis0[idbz];
                            basisProduct1=arrayBasis1[idbx]*arrayBasis0[idby]*arrayBasis0[idbz]/grid_space_x;
                            basisProduct2=arrayBasis0[idbx]*arrayBasis1[idby]*arrayBasis0[idbz]/grid_space_y;
                            basisProduct3=arrayBasis0[idbx]*arrayBasis0[idby]*arrayBasis1[idbz]/grid_space_z;

                            Jx=basisProduct1*(1+u_y_y[index])*(1+u_z_z[index])+basisProduct2*u_y_z[index]*u_z_x[index]+basisProduct3*u_y_x[index]*u_z_y[index]-basisProduct3*(1+u_y_y[index])*u_z_x[index]-basisProduct1*u_y_z[index]*u_z_y[index]-basisProduct2*u_y_x[index]*(1+u_z_z[index]);

                            Jy=(1+u_x_x[index])*basisProduct2*(1+u_z_z[index])+u_x_y[index]*basisProduct3*u_z_x[index]+u_x_z[index]*basisProduct1*u_z_y[index]-u_x_z[index]*basisProduct2*u_z_x[index]-(1+u_x_x[index])*basisProduct3*u_z_y[index]-u_x_y[index]*basisProduct1*(1+u_z_z[index]);

                            Jz=(1+u_x_x[index])*(1+u_y_y[index])*basisProduct3+u_x_y[index]*u_y_z[index]*basisProduct1+u_x_z[index]*u_y_x[index]*basisProduct2-u_x_z[index]*(1+u_y_y[index])*basisProduct1-(1+u_x_x[index])*u_y_z[index]*basisProduct2-u_x_y[index]*u_y_x[index]*basisProduct3;

                            /*if (index==618)
                            {
                                std::cout<<index<<" derivative: "<<std::endl;
                                std::cout<<basisProduct0<<", "<<basisProduct1<<", "<<basisProduct2<<", "<<basisProduct3<<std::endl;
                                std::cout<<Jx<<", "<<Jy<<", "<<Jz<<std::endl;
                            }*/


                            //coeff global index (storage form)
                            coeffIndex=(k+1)*coeffPlaneSize+(j+1)*num_element_x+(i+1);

                            grad(coeffIndex+1+bias)+=sstvdx2*diffVal*((image1_defVal+1000.0)*Jx*Je[index]+JacVal*grad_x_T_defVal*basisProduct0)/1055.0;
                            grad(coeffIndex+basisSize+1+bias)+=sstvdx2*diffVal*((image1_defVal+1000.0)*Jy*Je[index]+JacVal*grad_y_T_defVal*basisProduct0)/1055.0;
                            grad(coeffIndex+2*basisSize+1+bias)+=sstvdx2*diffVal*((image1_defVal+1000.0)*Jz*Je[index]+JacVal*grad_z_T_defVal*basisProduct0)/1055.0;

                        }
                        idbx++;
                    }
                    idby++;
                }
                idbz++;
            }


        }

    }


}


