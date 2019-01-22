/*
 * =====================================================================================
 *
 *       Filename:  LmkProcess.h
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  11/03/2009 03:31:16 PM CST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

int ReadLmkFile( char* filename, int* lmk_x, int* lmk_y, int* lmk_z )
{
    std::cout<<"Reading Lmk File: "<<filename<<std::endl;
    FILE *fp;
    fp=fopen(filename, "r");
    int num_lmk;
    fscanf(fp, "%d\n",&num_lmk);

    for ( int l=0; l<num_lmk; l++ )
    {
        fscanf(fp, "%d\t%d\t%d\n",&lmk_x[l],&lmk_y[l],&lmk_z[l]);
    }

    fclose(fp);
    return num_lmk;
}

int ReadLmkFile( const char* filename, int* lmk_x, int* lmk_y, int* lmk_z )
{
    std::cout<<"Reading Lmk File: "<<filename<<std::endl;
    FILE *fp;
    fp=fopen(filename, "r");
    int num_lmk;
    fscanf(fp, "%d\n",&num_lmk);

    for ( int l=0; l<num_lmk; l++ )
    {
        fscanf(fp, "%d\t%d\t%d\n",&lmk_x[l],&lmk_y[l],&lmk_z[l]);
    }

    fclose(fp);
    return num_lmk;
}

float CalcLmkErrorBefore( int num_lmk, int* x1, int* y1, int* z1, int* x2, int* y2, int* z2, float spacing_x, float spacing_y, float spacing_z, float* dist12 )
{
    float dx, dy, dz;
    float mean12 = 0;
    for ( int l=0; l<num_lmk; l++ )
    {
        dx = (x2[l]-x1[l])*spacing_x;
        dy = (y2[l]-y1[l])*spacing_y;
        dz = (z2[l]-z1[l])*spacing_z;

        dist12[l] = sqrt(static_cast<float>(dx*dx+dy*dy+dz*dz));
        mean12 += dist12[l];
    }

    mean12 /= num_lmk;
    return mean12;
}

float CalcLmkErrorAfter( int num_lmk, int* x1, int* y1, int* z1, int* x2, int* y2, int* z2, int xdim, int ydim, int zdim, float* disp12x, float* disp12y, float* disp12z, float spacing_x, float spacing_y, float spacing_z,float* dist12 )
{
    int planeSize=xdim*ydim;
    int index;
    float dx, dy, dz;
    float mean12 = 0;
    for ( int l=0; l<num_lmk; l++ )
    {
        index = z2[l]*planeSize+y2[l]*xdim+x2[l];
        dx = (x2[l]+disp12x[index] -x1[l])*spacing_x;
        dy = (y2[l]+disp12y[index] -y1[l])*spacing_y;
        dz = (z2[l]+disp12z[index] -z1[l])*spacing_z;

        dist12[l] = sqrt(static_cast<float>(dx*dx+dy*dy+dz*dz));
        mean12 += dist12[l];
    }

    //std::cout<<"total lmkerror: "<<mean12<<std::endl;
    mean12 /= num_lmk;
    return mean12;
}

bool CalcWriteLmkError( int num_lmk, int* x1, int* y1, int* z1, int* x2, int* y2, int* z2, int xdim, int ydim, int zdim, float* disp12x, float* disp12y, float* disp12z, float spacing_x, float spacing_y, float spacing_z, FILE* fp1, FILE* fp2 )
{
    int planeSize=xdim*ydim;
    int index;
    float dx, dy, dz;
    float mean = 0;
    float dx12, dy12, dz12;
    float mean12 = 0;

    float* dist=new float[num_lmk];
    float* dist12=new float[num_lmk];

    float* disp1=new float[num_lmk];
    float* disp2=new float[num_lmk];
    float* disp3=new float[num_lmk];

    for ( int l=0; l<num_lmk; l++ )
    {
        dx = (x2[l]-x1[l])*spacing_x;
        dy = (y2[l]-y1[l])*spacing_y;
        dz = (z2[l]-z1[l])*spacing_z;
        dist[l] = sqrt(static_cast<float>(dx*dx+dy*dy+dz*dz));
        mean += dist[l];

        index = z2[l]*planeSize+y2[l]*xdim+x2[l];      
        dx12 = (x2[l]+disp12x[index] -x1[l])*spacing_x;
        dy12 = (y2[l]+disp12y[index] -y1[l])*spacing_y;
        dz12 = (z2[l]+disp12z[index] -z1[l])*spacing_z;
        dist12[l] = sqrt(static_cast<float>(dx12*dx12+dy12*dy12+dz12*dz12));
        mean12 += dist12[l];

        disp1[l] = disp12x[index];
        disp2[l] = disp12y[index];
        disp3[l] = disp12z[index];
    }
    mean /= num_lmk;
    mean12 /= num_lmk;

    float tmp;
    float dev_sum=0;
    for (int l=0; l<num_lmk; l++)
    {
        tmp = dist12[l]-mean12;
        dev_sum += tmp*tmp;
    }
    float standard_dev = sqrt( dev_sum/num_lmk );

    float d1, d2, d3, dEuler;
    for ( int l=0; l<num_lmk; l++ )
    {
        d1 = x2[l]+disp1[l] -x1[l];
        d2 = y2[l]+disp2[l] -y1[l];
        d3 = z2[l]+disp3[l] -z1[l];
        dEuler = sqrt(static_cast<float>(d1*d1+d2*d2+d3*d3));
        fprintf(fp1,"%3d (%3d, %3d, %3d) + (%7.2f %7.2f %7.2f) - (%3d, %3d, %3d) = (%7.2f %7.2f %7.2f)\t%7.2f\n",l+1,x2[l],y2[l],z2[l],disp1[l],disp2[l],disp3[l],x1[l],y1[l],z1[l],d1,d2,d3,dEuler);
    }


    for ( int l=0; l<num_lmk; l++ )
    {
       fprintf(fp2,"%3d %7.2f %7.2f\n",l+1,dist[l],dist12[l]);
    }
    std::cout<<"Avg lmkerror before and after registration on "<<num_lmk<<" lmks:"<<std::endl;
    printf("%7.2f\t%7.2f\t%7.2f\n",mean,mean12,standard_dev);
    fprintf(fp2,"Avg lmkerror before and after registration on %d lmks:\n", num_lmk);
    fprintf(fp2,"%7.2f\t%7.2f\t%7.2f\n",mean,mean12,standard_dev);


    delete dist;
    delete dist12;

    delete disp1;
    delete disp2;
    delete disp3;

    return true;
}



bool WriteLmkDef( int num_lmk, int* x2, int* y2, int* z2, int xdim, int ydim, int zdim, float* disp12x, float* disp12y, float* disp12z, float spacing_x, float spacing_y, float spacing_z, FILE* fp2 )
{
    fprintf(fp2,"%d\n", num_lmk);

    int planeSize=xdim*ydim;
    int index;
    float x2d, y2d, z2d;
    for ( int l=0; l<num_lmk; l++ )
    {
        index = z2[l]*planeSize+y2[l]*xdim+x2[l];      
        x2d = x2[l]+disp12x[index];
        y2d = y2[l]+disp12y[index];
        z2d = z2[l]+disp12z[index];
        fprintf(fp2,"%6.2f   %6.2f   %6.2f\n",x2d,y2d,z2d);
    }

    return true;
}

