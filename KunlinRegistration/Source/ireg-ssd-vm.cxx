#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <time.h>
#include "lbfgsbcpp.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkBinaryThresholdImageFilter.h"

//#include "cubic_BSpline.h"
#include "CostOptimize.h"
#include "ImageType.h"
#include "IMGIO.h"
#include "BasicProcess.h"
#include "LmkProcess.h"
#include "ArrayProcess.h"

//////////////////////////////
#include "ConfigReader.h"
//////////////////////////////

using namespace std;

#define bgValue -1000
#define dbgValue 0
#define debugTag 0


#define precalcu(x,y,z) precalcu[(z)*16*grid_space_x*grid_space_y+(y)*4*grid_space_x+x]
#define precalcux(x,y,z) precalcux[(z)*16*grid_space_x*grid_space_y+(y)*4*grid_space_x+x]
#define precalcuy(x,y,z) precalcuy[(z)*16*grid_space_x*grid_space_y+(y)*4*grid_space_x+x]
#define precalcuz(x,y,z) precalcuz[(z)*16*grid_space_x*grid_space_y+(y)*4*grid_space_x+x]
#define precalcu2x(x,y,z) precalcu2x[(z)*16*grid_space_x*grid_space_y+(y)*4*grid_space_x+x]
#define precalcu2y(x,y,z) precalcu2y[(z)*16*grid_space_x*grid_space_y+(y)*4*grid_space_x+x]
#define precalcu2z(x,y,z) precalcu2z[(z)*16*grid_space_x*grid_space_y+(y)*4*grid_space_x+x]


//////////////////////////////////////////////////////////////////////////////////////////////////////
//This function should return cost and gradients 
void funcgrad(const ap::real_1d_array& inputx, float& f, ap::real_1d_array& grad, float* arrayBasis0, float* arrayBasis1, float* arrayBasis2, short int* arrayLocation, int xdim, int ydim, int zdim, int num_element_x, int num_element_y, int num_element_z, int grid_space_x, int grid_space_y, int grid_space_z, float* imagein1, float* imagein2, float* grad_x_T, float* grad_y_T, float* grad_z_T, float* c_x, float* c_y, float* c_z, float* u_x_existing, float* u_y_existing, float* u_z_existing, float* u_x_total, float* u_y_total, float* u_z_total, unsigned char* mask, int maskCT2, int num_lmk, int* xc_1, int* yc_1, int* zc_1, int* xc_2, int* yc_2, int* zc_2, float spacing_x, float spacing_y, float spacing_z, float sstvd, float param, float smooth, float alpha, float beta, float gamma, float& cost_sstvd, float& cost_lmk, float& cost_smooth, float* vm1, float* vm2, float* vm_grad_x_T, float* vm_grad_y_T, float* vm_grad_z_T, float ved, float& cost_ved, float& vedcoeff, float* precalcu, float* precalcux, float* precalcuy, float* precalcuz, float* precalcu2x, float* precalcu2y, float* precalcu2z, float& JacMin, float& JacMax)
{
    //std::cout<<"here 0"<<std::endl;

    int basisSize=num_element_x*num_element_y*num_element_z;
    int imgSize=xdim*ydim*zdim;
    
    /*float* vdeformed=new float[imgSize];
    float* grad_x_Tt=new float[imgSize];
	float* grad_y_Tt=new float[imgSize];
	float* grad_z_Tt=new float[imgSize];*/

    float* u_x_x=new float[imgSize];
    float* u_x_y=new float[imgSize];
    float* u_x_z=new float[imgSize];
    float* u_y_x=new float[imgSize];
    float* u_y_y=new float[imgSize];
    float* u_y_z=new float[imgSize];
    float* u_z_x=new float[imgSize];
    float* u_z_y=new float[imgSize];
    float* u_z_z=new float[imgSize];

    float* u_x_xy=new float[imgSize];
    float* u_x_xz=new float[imgSize];
    float* u_y_xy=new float[imgSize];
    float* u_y_yz=new float[imgSize];
    float* u_z_xz=new float[imgSize];
    float* u_z_yz=new float[imgSize];
    

    float* tmp=new float[imgSize];
    //for (int ii=0; ii<imgSize; ii++)
      //  tmp[ii]=0;
    /*float* u_x_x_def=new float[imgSize];
    float* u_x_y_def=new float[imgSize];
    float* u_x_z_def=new float[imgSize];
    float* u_y_x_def=new float[imgSize];
    float* u_y_y_def=new float[imgSize];
    float* u_y_z_def=new float[imgSize];
    float* u_z_x_def=new float[imgSize];
    float* u_z_y_def=new float[imgSize];
    float* u_z_z_def=new float[imgSize];*/

    /*float* uc_x=new float[imgSize];
    float* uc_y=new float[imgSize];
    float* uc_z=new float[imgSize];*/

    /*float* vm_vdeformed=new float[imgSize];
    float* vm_grad_x_Tt=new float[imgSize];
	float* vm_grad_y_Tt=new float[imgSize];
	float* vm_grad_z_Tt=new float[imgSize];*/


    int index=0;
    //std::cout<<"func started"<<std::endl;

    for (int k=0;k<num_element_z;k++)
        for (int j=0;j<num_element_y;j++)
            for (int i=0;i<num_element_x;i++)
            {
                c_x[index]=inputx(index+1);
                c_y[index]=inputx(index+basisSize+1);
                c_z[index]=inputx(index+2*basisSize+1);
                index++;
            }

    //initialize gradients and cost
	for (int i=0;i<3*basisSize;i++)
        grad(i+1)=0;

    //set each cost=0;
    cost_sstvd=0;
    cost_ved=0;
    cost_lmk=0;
    //cost_lap=0;
    cost_smooth=0;
    JacMin=1;//100.0;
    JacMax=1;//0.0;
    vedcoeff=ved;



    compose_u9(u_x_total, u_y_total, u_z_total, u_x_existing, u_y_existing, u_z_existing, arrayBasis0, arrayLocation, c_x, c_y, c_z, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, maskCT2);

    /*for(int ii=0; ii<imgSize; ii++)
    {
        u_x_total[ii]=0.0;
        u_y_total[ii]=0.0;
        u_z_total[ii]=0.0;
    }*/

    ///////////////////////////////////////////////////////////////////////////////////

    //compute 1st order derivatives (step 1) --- current disp
    //interpolate_uxyz(u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, c_x, c_y, c_z, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, grid_space_x, grid_space_y, grid_space_z, precalcux, precalcuy, precalcuz);
    
    //deform derivatives (step 2)
    //deformed0_derivatives(u_x_x_def, u_x_y_def, u_x_z_def, u_y_x_def, u_y_y_def, u_y_z_def, u_z_x_def, u_z_y_def, u_z_z_def, u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);

    /*deformed1_inner(u_x_x, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);
    deformed1_inner(u_x_y, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);
    deformed1_inner(u_x_z, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);

    deformed1_inner(u_y_x, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);
    deformed1_inner(u_y_y, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);
    deformed1_inner(u_y_z, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);

    deformed1_inner(u_z_x, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);
    deformed1_inner(u_z_y, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);
    deformed1_inner(u_z_z, tmp, xdim, ydim, zdim, u_x_existing, u_y_existing, u_z_existing, dbgValue);*/


    ///////////////////////////////////////////////////////////////////////////////////
    //these operations have been moved into optimize func to save memory
    //
    //image1 deformed using total disp
    //deformed0(vdeformed, imagein1, xdim, ydim, zdim, u_x_total, u_y_total, u_z_total, bgValue);

    //image1_gradient deformed using total disp
    //deformed0(grad_x_Tt, grad_x_T, xdim, ydim, zdim, u_x_total, u_y_total, u_z_total, dbgValue);
    //deformed0(grad_y_Tt, grad_y_T, xdim, ydim, zdim, u_x_total, u_y_total, u_z_total, dbgValue);
    //deformed0(grad_z_Tt, grad_z_T, xdim, ydim, zdim, u_x_total, u_y_total, u_z_total, dbgValue);

    ///////////////////////////////////////////////////////////////////////////////////
    //ssd optimize
    ssd_optimize(grad, arrayBasis0, arrayLocation, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, imagein2, imagein1, grad_x_T, grad_y_T, grad_z_T, u_x_total, u_y_total, u_z_total, sstvd, cost_sstvd, mask, maskCT2, 0);
    //sstvd_optimize(grad, arrayBasis0, arrayBasis1, arrayLocation, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, imagein2, imagein1, grad_x_T, grad_y_T, grad_z_T, u_x_total, u_y_total, u_z_total, u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, sstvd, cost_sstvd, mask, maskCT2, 0, Je, grid_space_x, grid_space_y, grid_space_z, JacMin, JacMax);

    //std::cout<<"here 1"<<std::endl;
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    //ved-optimize

    //vm1 deformed using total disp
    //deformed0(vm_vdeformed, vm1, xdim, ydim, zdim, u_x_total, u_y_total, u_z_total, 0);

    if (fabs(ved)>1e-8)
    {
        //vm1_gradient deformed using total disp
        //deformed0(vm_grad_x_Tt, vm_grad_x_T, xdim, ydim, zdim, u_x_total, u_y_total, u_z_total, 0);
        //deformed0(vm_grad_y_Tt, vm_grad_y_T, xdim, ydim, zdim, u_x_total, u_y_total, u_z_total, 0);
        //deformed0(vm_grad_z_Tt, vm_grad_z_T, xdim, ydim, zdim, u_x_total, u_y_total, u_z_total, 0);

        //float ved_weight=ved*cost_sstvd;
        //ved_optimize_ratio(grad, arrayBasis0, arrayLocation, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, vm2, vm1, vm_grad_x_T, vm_grad_y_T, vm_grad_z_T, u_x_total, u_y_total, u_z_total, ved_weight, cost_ved, vedcoeff, mask, maskCT2, 0);
        ved_optimize(grad, arrayBasis0, arrayLocation, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, vm2, vm1, vm_grad_x_T, vm_grad_y_T, vm_grad_z_T, u_x_total, u_y_total, u_z_total, ved, cost_ved, mask, maskCT2, 0);
    }
    else
    {
        ved_calc(arrayLocation, vm2, vm1, u_x_total, u_y_total, u_z_total, xdim, ydim, zdim, cost_ved, mask, maskCT2);
    }

    //std::cout<<"here 2"<<std::endl;

    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    //lap optimize

    //revise: use 2nd order interpolation
    differentiate2(u_x_total, u_x_x, u_x_y, u_x_z, xdim, ydim, zdim);
    differentiate2(u_y_total, u_y_x, u_y_y, u_y_z, xdim, ydim, zdim);
    differentiate2(u_z_total, u_z_x, u_z_y, u_z_z, xdim, ydim, zdim);

    differentiate2_crossX(u_x_total, tmp, u_x_xy, u_x_xz, xdim, ydim, zdim);
    differentiate2_crossY(u_y_total, tmp, u_y_xy, u_y_yz, xdim, ydim, zdim);
    differentiate2_crossZ(u_z_total, tmp, u_z_xz, u_z_yz, xdim, ydim, zdim);
    

    /*if (fabs(lap)>1e-8)
    {
        lap_optimize(grad, arrayBasis0, arrayBasis2, arrayLocation, u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, lap, cost_lap, mask, maskCT2, 0, grid_space_x, grid_space_y, grid_space_z);
    }
    else
    {
        lap_calc(arrayLocation, u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, xdim, ydim, zdim, cost_lap, mask, maskCT2);
    }*/


    //if (fabs(param)>1e-8)
    if (fabs(smooth)>1e-8)
    {
        //smooth=0;
        //std::cout<<"here smooth_optimize"<<std::endl;
        smooth_optimize(grad, arrayBasis0, arrayBasis1, arrayBasis2, arrayLocation, u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, u_x_xy, u_x_xz, u_y_xy, u_y_yz, u_z_xz, u_z_yz, u_x_total, u_y_total, u_z_total, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, smooth, alpha, beta, gamma, cost_smooth, mask, maskCT2, 0, grid_space_x, grid_space_y, grid_space_z);
        //lap_optimize(grad, arrayBasis0, arrayBasis2, arrayLocation, u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, smooth, cost_smooth, mask, maskCT2, 0, grid_space_x, grid_space_y, grid_space_z);
    }
    else
    {
        //std::cout<<"here smooth_calc"<<std::endl;
        smooth_clac(arrayLocation, u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, u_x_xy, u_x_xz, u_y_xy, u_y_yz, u_z_xz, u_z_yz, u_x_total, u_y_total, u_z_total, xdim, ydim, zdim, alpha, beta, gamma, cost_smooth, mask, maskCT2);
        //lap_calc(arrayLocation, u_x_x, u_x_y, u_x_z, u_y_x, u_y_y, u_y_z, u_z_x, u_z_y, u_z_z, xdim, ydim, zdim, cost_smooth, mask, maskCT2);
    }

    //std::cout<<"here 3: "<<cost_smooth<<std::endl;
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    //lmk optimize
    if (fabs(param)>1e-8)
    {
        lmk_optimize(grad, arrayBasis0, arrayLocation, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, num_lmk, xc_1, yc_1, zc_1, xc_2, yc_2, zc_2, u_x_total, u_y_total, u_z_total, param, cost_lmk, maskCT2, 0);
    }
    else
    {
        lmk_calc(xdim, ydim, zdim, num_lmk, xc_1, yc_1, zc_1, xc_2, yc_2, zc_2, u_x_total, u_y_total, u_z_total, spacing_x, spacing_y, spacing_z, cost_lmk);
    }


    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    //calc total cost

    f = sstvd*cost_sstvd + vedcoeff*cost_ved + smooth*cost_smooth + param*cost_lmk;
    //std::cout<<f<<std::endl;


    /*delete [] vdeformed;
    delete [] grad_x_Tt;
    delete [] grad_y_Tt;
    delete [] grad_z_Tt;*/

    delete [] u_x_x;
    delete [] u_x_y;
    delete [] u_x_z;
    delete [] u_y_x;
    delete [] u_y_y;
    delete [] u_y_z;
    delete [] u_z_x;
    delete [] u_z_y;
    delete [] u_z_z;
    //std::cout<<"here 4"<<std::endl;

    delete [] u_x_xy;
    delete [] u_x_xz;
    delete [] u_y_xy;
    delete [] u_y_yz;
    delete [] u_z_xz;
    delete [] u_z_yz;
    //std::cout<<"here 5"<<std::endl;

    delete [] tmp;
    /*delete [] u_x_x_def;
    delete [] u_x_y_def;
    delete [] u_x_z_def;
    delete [] u_y_x_def;
    delete [] u_y_y_def;
    delete [] u_y_z_def;
    delete [] u_z_x_def;
    delete [] u_z_y_def;
    delete [] u_z_z_def;*/

    /*delete [] uc_x;
    delete [] uc_y;
    delete [] uc_z;*/

    /*delete [] vm_vdeformed;
    delete [] vm_grad_x_Tt;
    delete [] vm_grad_y_Tt;
    delete [] vm_grad_z_Tt;*/

}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{
    //reading paramters
    Configuration config;
    //std::cout<<"here 0"<<std::endl;
    config.ReadConfiguration( argv[1] );
    //std::cout<<"here 1"<<std::endl;
    //exit(0);

    std::cout<<"################################################################"<<std::endl;
    std::cout<<"                     Reading Parameters Done                    "<<std::endl;
    std::cout<<"################################################################"<<std::endl<<std::endl;


///////////////////////////////////////////////////////////////////////////////////
    
    //std::cout<<"size of int is: "<<sizeof(int)<<std::endl; //4
    //std::cout<<"size of short int is: "<<sizeof(short int)<<std::endl; //2
    //std::cout<<"size of float is: "<<sizeof(float)<<std::endl; //4

    FILE *fpr, *fp, *fpc;
    fpr=fopen(argv[1],"r");

    int ifRead = config.if_read_disp;
    string bfr = config.initial_disp_dir;
    string sufr = config.initial_disp_suffix + ".mha";

    int lowerThreshold = config.mask_lower_thr;
    int upperThreshold = config.mask_upper_thr;
    std::cout<<"lowerThreshold: "<<lowerThreshold<<std::endl;
    std::cout<<"upperThreshold: "<<upperThreshold<<std::endl;

    string img1Name = config.moving_image_file;
    string img2Name = config.fixed_image_file;
    string mask1Name = config.moving_mask_file;
    string mask2Name = config.fixed_mask_file;
    string vm1Name = config.moving_vm_file;
    string vm2Name = config.fixed_vm_file;

    int readVM = config.readVM;

    const char* dirName = config.result_dir.c_str();
    const char* fileNameResult = config.log_filename.c_str();
    const char* fileNameCoeff = config.coeff_filename.c_str();

    const char* fileNameLmk1 = config.moving_lmk_file.c_str();
    const char* fileNameLmk2 = config.fixed_lmk_file.c_str();
    const char* fileNameLMKResult = config.lmkerror_filename.c_str();

    
    int total_level = config.total_level_number;
    int savefull = config.if_save_fullsize_disp;
    int writeChoice = config.if_write_internal_image;

    float wt_tv = config.wt_tv;
    float wt_vm = config.wt_vm;
	//float wt_lap = config.wt_lap;
    float wt_smooth = config.wt_smooth;
    float alpha = config.alpha;
    float beta = config.beta;
    float gamma = config.gamma;
	float wt_lmk = config.wt_lmk;

    int blutScale = config.blut_scale;
    int xgrid=blutScale;
    int ygrid=blutScale;
    int zgrid=blutScale;

    float sigmaMin_global = config.vm_sigma_min_global;
    float sigmaMax_global = config.vm_sigma_max_global;
    int nScales_global = config.vm_iteration_global;
    ////////////////////////////////////////////////////////////////////////////////////////////
    //Read HU and mask images
    std::cout<<"Reading images ..."<<std::endl;


    std::cout<<img1Name<<std::endl;
    std::cout<<img2Name<<std::endl;
    std::cout<<mask1Name<<std::endl;
    std::cout<<mask2Name<<std::endl;

    ImageType::Pointer iimagein1 = ImageIO.ReadImage< ImageType >( img1Name );
    ImageType::Pointer iimagein2 = ImageIO.ReadImage< ImageType >( img2Name );
    
    MaskImageType::Pointer imaskin1 = ImageIO.ReadImage< MaskImageType >( mask1Name );
    MaskImageType::Pointer imaskin2 = ImageIO.ReadImage< MaskImageType >( mask2Name );

    ////////////////////////////////////////////////////////////////////////////////////////////
    //threshold masks, and generate 0-1 masks
    std::cout<<"Thresholding masks ..."<<std::endl;

    const MaskPixelType insideValue = 1;
    const MaskPixelType outsideValue = 0;

    const MaskPixelType dilateValue = 1;

    MaskImageType::Pointer imask1 = ThresholdImage< MaskImageType >( imaskin1, lowerThreshold, upperThreshold, insideValue, outsideValue );
    MaskImageType::Pointer imask2 = ThresholdImage< MaskImageType >( imaskin2, lowerThreshold, upperThreshold, insideValue, outsideValue );

    imaskin1 = NULL;
    imaskin2 = NULL;

    ////////////////////////////////////////////////////////////////////
    //mask input images, and write out masked images (registration input HU images)
    std::cout<<"Masking HU images ..."<<std::endl;

    const MaskImageType::PixelType unmarkValue = outsideValue;
    ImageType::Pointer iimage1_o = MaskImage< ImageType, MaskImageType >( iimagein1, imask1, bgValue, unmarkValue );
    ImageType::Pointer iimage2 = MaskImage< ImageType, MaskImageType >( iimagein2, imask2, bgValue, unmarkValue );

    iimagein1 = NULL;
    iimagein2 = NULL;

    ////////////////////////////////////////////////////////////////////
    //histmatching template to target HU range

    ImageType::Pointer iimage1 = HistMatching_bg< ImageType, MaskImageType >( iimage1_o, iimage2, imask1, imask2, bgValue, unmarkValue );
    std::cout<<"HistMatching Image1 to Image2 Done"<<std::endl;

    //ImageType::Pointer iimage1 = iimage1_o;

    ////////////////////////////////////////////////////////////////////
    //write masked HU images

    string bf = config.result_dir;
    std::cout<<bf<<std::endl;
    
    //if ( writeChoice==1 )
    {
        std::cout<<"Writing Template and Target images ..."<<std::endl; 

        string wfn1_o = bf + "Image1_o.mha";
        std::cout<<wfn1_o<<std::endl;
        ImageIO.CastWriteImage< ImageType, HUImageType >( iimage1_o, wfn1_o );

        string wfn1 = bf + "Image1.mha";
        string wfn2 = bf + "Image2.mha";
        std::cout<<wfn1<<std::endl;
        std::cout<<wfn2<<std::endl;

        ImageIO.CastWriteImage< ImageType, HUImageType >( iimage1, wfn1 );
        ImageIO.CastWriteImage< ImageType, HUImageType >( iimage2, wfn2 );
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //get image basic information

    ImageType::DirectionType direction = iimage1->GetDirection();
    //std::cout<<iimage1->GetDirection()<<std::endl;

    ImageType::SizeType sizeFull = iimage1->GetRequestedRegion().GetSize();
    //std::cout<<iimage1->GetRequestedRegion().GetSize()<<std::endl;
    ImageType::SizeType size;
	ImageType::SpacingType spacingFull = iimage1->GetSpacing();
	ImageType::SpacingType spacing;

    ImageType::IndexType iStart;
    iStart.Fill(0);
    
    ImageType::PointType origin = iimage1->GetOrigin();
	/*origin[0]=0.0;
	origin[1]=0.0;
	origin[2]=0.0;*/
	std::cout<<"Origin: "<<origin<<std::endl;

	ImageType::RegionType regionFull;
	regionFull.SetIndex(iStart);
	regionFull.SetSize(sizeFull);

    int xdimo=sizeFull[0];
    int ydimo=sizeFull[1];
    int zdimo=sizeFull[2];
    int imgSizeFull=xdimo*ydimo*zdimo;
    //int planeSizeFull=xdimo*ydimo;
	
    float spacingFull_x = static_cast<float>(spacingFull[0]);
    float spacingFull_y = static_cast<float>(spacingFull[1]);
    float spacingFull_z = static_cast<float>(spacingFull[2]);
    
    /////////////////////////////////////////////////////
    //initialize arrays
    std::cout<<"Initializing arrays ..."<<std::endl; 

    float* image1HUo_o=new float[imgSizeFull];
    float* image1HUo=new float[imgSizeFull];
	float* image2HUo=new float[imgSizeFull];
    
    unsigned char* mask1o=new unsigned char[imgSizeFull];
    unsigned char* mask1o_def=new unsigned char[imgSizeFull];
    unsigned char* mask2o=new unsigned char[imgSizeFull];
    unsigned char* maskUniono=new unsigned char[imgSizeFull];
    //used to output deformed image
	float* vdeformed1=new float[imgSizeFull];

    
    float sigmaMin=2.0;
    float sigmaMax=3.0;
    int nScales=10;
    //float sigmaMin=sigmaMin_global;
    //float sigmaMax=sigmaMax_global;
    //int nScales=nScales_global;
    //ImageType::Pointer iVMImage1 = CalcVesselnessMeasure< ImageType, ImageType, MaskImageType >( iimage1, imask1, sigmaMin, sigmaMax, nScales ); 
    //ImageType::Pointer iVMImage2 = CalcVesselnessMeasure< ImageType, ImageType, MaskImageType >( iimage2, imask2, sigmaMin, sigmaMax, nScales );
    ImageType::Pointer iVMImage1;
    ImageType::Pointer iVMImage2;
    if (readVM>0)
    {
        std::cout<<"############################"<<std::endl;
        std::cout<<vm1Name<<std::endl;
        std::cout<<"############################"<<std::endl;
        iVMImage1 = ImageIO.ReadImage< ImageType >( vm1Name );
        iVMImage2 = ImageIO.ReadImage< ImageType >( vm2Name );
    }
    else
    {
        iVMImage1 = CalcVesselnessMeasure< ImageType, ImageType, MaskImageType >( iimage1, imask1, sigmaMin, sigmaMax, nScales ); 
        iVMImage2 = CalcVesselnessMeasure< ImageType, ImageType, MaskImageType >( iimage2, imask2, sigmaMin, sigmaMax, nScales );
    }
    float* vm1o=new float[imgSizeFull];
	float* vm2o=new float[imgSizeFull];
        
    //if ( writeChoice==1 )
    {
        std::cout<<"Writing Template and Target images Vesselness Measurement ..."<<std::endl; 

        string wfn1 = bf + "Image1VM.mha";
        string wfn2 = bf + "Image2VM.mha";
        std::cout<<wfn1<<std::endl;
        std::cout<<wfn2<<std::endl;

        ImageIO.WriteImage< ImageType >( iVMImage1, wfn1 );
        ImageIO.WriteImage< ImageType >( iVMImage2, wfn2 );
    }


	/*iimage1=ImageIO.CastImage<MaskImageType, ImageType>( imask1 );
	iimage2=ImageIO.CastImage<MaskImageType, ImageType>( imask2 );

	string wfn1 = bf + "Image1.mha";
        string wfn2 = bf + "Image2.mha";
	ImageIO.CastWriteImage< ImageType, HUImageType >( iimage1, wfn1 );
        ImageIO.CastWriteImage< ImageType, HUImageType >( iimage2, wfn2 );
    */



    /*float* image1HU=new float[imgSizeFull];
	float* image2HU=new float[imgSizeFull];
    unsigned char* mask1=new unsigned char[imgSizeFull];
    unsigned char* mask2=new unsigned char[imgSizeFull];
    unsigned char* maskUnion=new unsigned char[imgSizeFull];

    float* grad_x_T=new float[imgSizeFull];
	float* grad_y_T=new float[imgSizeFull];
	float* grad_z_T=new float[imgSizeFull];

    /////////////////////////////////////////////
    //for vesselness images
    float* vm1=new float[imgSizeFull];
	float* vm2=new float[imgSizeFull];
    float* vm1def=new float[imgSizeFull];
    float* vm_grad_x_T=new float[imgSizeFull];
	float* vm_grad_y_T=new float[imgSizeFull];
	float* vm_grad_z_T=new float[imgSizeFull];*/
    /////////////////////////////////////////////

    float* u_x_for_total=new float[imgSizeFull];
    float* u_y_for_total=new float[imgSizeFull];
    float* u_z_for_total=new float[imgSizeFull];

    float* u_x_for=new float[imgSizeFull];
    float* u_y_for=new float[imgSizeFull];
    float* u_z_for=new float[imgSizeFull];

	float* u_x_for_existing_fullsize=new float[imgSizeFull];
    float* u_y_for_existing_fullsize=new float[imgSizeFull];
    float* u_z_for_existing_fullsize=new float[imgSizeFull];

    float* u_x_for_existing=new float[imgSizeFull];
    float* u_y_for_existing=new float[imgSizeFull];
    float* u_z_for_existing=new float[imgSizeFull];

	time_t start,end;
    int* tdiff=new int[total_level+1];
    int total_time=0;
	
    ///////////////////////////////////////////////////////////////////////////////////////////
    //initialize Jac related images

    float* Jac_curr=new float[imgSizeFull];
    float* Jac_curr_def=new float[imgSizeFull];
    float* Jac_exis=new float[imgSizeFull];
    float* Jac_exis_updated=new float[imgSizeFull];

    //float* Je=new float[imgSizeFull];

    ///////////////////////////////////////////////////////////////////////////////////////////

    int arraySize=xgrid*4;
    float* arrayBasis0 = new float[arraySize];
    float* arrayBasis1 = new float[arraySize];
    float* arrayBasis2 = new float[arraySize];

    ///////////////////////////////////////////////////////////////////////////////////////////
    //convert input images to arrays

	std::cout<<"Converting HU ITKImages to Arrays done"<<std::endl;

    ITKImage2Array< ImageType >( iimage1_o, image1HUo_o );
    iimage1_o = NULL;
    ITKImage2Array< ImageType >( iimage1, image1HUo );
    ITKImage2Array< ImageType >( iimage2, image2HUo );
    iimage1 = NULL;
    iimage2 = NULL;

    ITKImage2ArrayUCHAR< MaskImageType >( imask1, mask1o );
    ITKImage2ArrayUCHAR< MaskImageType >( imask2, mask2o );
    imask1 = NULL;
    imask2 = NULL;

    ITKImage2Array< ImageType >( iVMImage1, vm1o );
    ITKImage2Array< ImageType >( iVMImage2, vm2o );
    iVMImage1 = NULL;
    iVMImage2 = NULL;

    ///////////////////////////////////////////////////////////////////////////////////////////
    //read landmarks

    int num_lmk;
    fp=fopen(fileNameLmk1,"r");
    fscanf(fp,"%d\n",&num_lmk);
    fclose(fp);
    
    int* xc_1o=new int[num_lmk];
    int* yc_1o=new int[num_lmk];
    int* zc_1o=new int[num_lmk];
    ReadLmkFile(fileNameLmk1, xc_1o, yc_1o, zc_1o );
    int* xc_2o=new int[num_lmk];
    int* yc_2o=new int[num_lmk];
    int* zc_2o=new int[num_lmk];
    if ( num_lmk != ReadLmkFile(fileNameLmk2, xc_2o, yc_2o, zc_2o ) )
    {
        std::cout<<"Lmk numbers are different!"<<std::endl;
        return false;
    }

    float* dist=new float[num_lmk];
    float avglmkdist = CalcLmkErrorBefore( num_lmk, xc_1o, yc_1o, zc_1o, xc_2o, yc_2o, zc_2o, spacingFull_x, spacingFull_y, spacingFull_z, dist );
    std::cout<<"Avg lmkerror before registration..."<<std::endl;
    std::cout<<avglmkdist<<std::endl;
    
    float* dist12=new float[num_lmk];

    int* xc_1=new int[num_lmk];
    int* yc_1=new int[num_lmk];
    int* zc_1=new int[num_lmk];
    int* xc_2=new int[num_lmk];
    int* yc_2=new int[num_lmk];
    int* zc_2=new int[num_lmk];

    ///////////////////////////////////////////////////////////////////////////////////////////
    //read existing disps

	if ( ifRead == 1 )
    {
		/*string bfr(dispRead);
		//std::cout<<"Read in: "<<dispRead<<std::endl;
		std::cout<<"Read in: "<<bfr<<std::endl;

		DispFieldType::Pointer initial_dfield = ImageIO.ReadImage< DispFieldType >(bfr);
        DispField2Arrays< DispFieldType >( initial_dfield, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize );*/

        std::cout<<"Read in: "<<sufr<<std::endl;
        string fn3 = bfr + "Disp12_x" + sufr;
        string fn4 = bfr + "Disp12_y" + sufr;
        string fn5 = bfr + "Disp12_z" + sufr;

        DispImageType::Pointer imagein3 = ImageIO.ReadImage< DispImageType >( fn3 );
        DispImageType::Pointer imagein4 = ImageIO.ReadImage< DispImageType >( fn4 );
        DispImageType::Pointer imagein5 = ImageIO.ReadImage< DispImageType >( fn5 );
		
        DispITKImage2Array< DispImageType >( imagein3, imagein4, imagein5, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize );

        imagein3 = NULL;
        imagein4 = NULL;
        imagein5 = NULL;
    }
    else
    {
        int id=0;
        for ( int z=0; z<zdimo; z++ )
            for ( int y=0; y<ydimo; y++ )
                for ( int x=0; x<xdimo; x++ )
                {
                    u_x_for_existing_fullsize[id]=0;
                    u_y_for_existing_fullsize[id]=0;
	    			u_z_for_existing_fullsize[id]=0;
                    id++;
                }
        if ( id!=imgSizeFull )
        {
            std::cout<<"Initializing disp problem! Exiting..."<<std::endl;
            exit(3);
        }
    }

    //////////////////////////////////////////////////////////////////////////

    /*typedef itk::MinimumMaximumImageCalculator<ImageType> RangeCalculatorType;
    RangeCalculatorType::Pointer calculator = RangeCalculatorType::New(); 
	calculator->SetImage( iimage1 );
	calculator->Compute();
	//int MaxValue = calculator->GetMaximum();
	int MinValue = calculator->GetMinimum();

    ImageType::PixelType DefaultPixelValue = MinValue;
	std::cout<<"DefaultPixel Value is: "<<DefaultPixelValue<<std::endl;*/
   
    /*const std::string nearestInterpolation = "n";
    const std::string linearInterpolation = "l";
    const std::string bsplineInterpolation = "b";*/

    //////////////////////////////////////////////////////////////////////////
    //for each level
    
    char buffer[255];
    strcpy(buffer,dirName);
    strcat(buffer,fileNameResult);
    fp=fopen(buffer,"w");


    char buffer2[255];
    strcpy(buffer2,dirName);
    strcat(buffer2,fileNameCoeff);
    std::cout<<buffer2<<std::endl;
    fpc=fopen(buffer2, "wb");

    fwrite(&total_level,sizeof(int),1,fpc);
    fwrite(&xdimo,sizeof(int),1,fpc);
    fwrite(&ydimo,sizeof(int),1,fpc);
    fwrite(&zdimo,sizeof(int),1,fpc);

    
    char buffer0[255];
    strcpy(buffer0,dirName);
    int res=1000;
    int res_iter=0;
    //char* tmpStr;
    string statFileName;
    FILE* fpStat;


    int if_write_disp;
	char suffix[255];

    int scale_x, scale_y, scale_z;

    int grid_space_x, grid_space_y, grid_space_z;
	
    float eps1, eps2, eps3;

	int maxIter;

    int dilateRadius;

    //float sigmaMin,sigmaMax;
    //int nScales;

    float vmConstWeight = wt_vm;

    
    float cost_sstvd,cost_ved,cost_smooth,cost_lmk;
    float JacMin, JacMax;


    //std::cout<<"here each"<<std::endl;
    fprintf(fp, "prog = ireg-ssd-vm.exe   level = %d   dir = %s\n", total_level, dirName);
    std::cout<<"wt_tv="<<wt_tv<<", "<<"wt_vm="<<wt_vm<<", "<<"wt_smooth="<<wt_smooth<<", "<<"wt_lmk="<<wt_lmk<<std::endl;
    fprintf(fp, "wt_tv=%f, wt_vm=%f, wt_smooth=%f, wt_lmk=%f\n", wt_tv, wt_vm, wt_smooth, wt_lmk);
    for (int level=0;level<total_level;level++)
	{
		time(&start);

        if_write_disp = config.level_setting[level].if_write_disp;
        scale_x = config.level_setting[level].image_size_ratio[0];
        scale_y = config.level_setting[level].image_size_ratio[1];
        scale_z = config.level_setting[level].image_size_ratio[2];
		
        grid_space_x = config.level_setting[level].spline_spacing[0];
        grid_space_y = config.level_setting[level].spline_spacing[1];
        grid_space_z = config.level_setting[level].spline_spacing[2];

        eps1 = config.level_setting[level].epsilon[0];
        eps2 = config.level_setting[level].epsilon[1];
        eps3 = config.level_setting[level].epsilon[2];
	    
        maxIter = config.level_setting[level].iteration;
        
        //fscanf(fpr,"%d",&dilateRadius);
        //dilateRadius=static_cast<int>(ceil(grid_space_x/2.479472335));
        dilateRadius = config.level_setting[level].dilate_radius;
    
        sigmaMin = config.level_setting[level].vm_sigma_min;
        sigmaMax = config.level_setting[level].vm_sigma_max;
        nScales = config.level_setting[level].vm_iteration;
		
        //if (level == 0) //level==0 no ssvmd
        /*if (level == 0 || scale_x == 8) //scale_x==8 and level==0 no ssvmd
            wt_vm = 0;
        else 
            wt_vm = vmConstWeight;*/
        wt_vm = vmConstWeight;
        std::cout<<"scale_x = "<<scale_x<<std::endl;
        std::cout<<"wt_vm = "<<wt_vm<<std::endl;

        sprintf(suffix,"_%d_%d_%d_%d_%d_%d.mha",scale_x,scale_y,scale_z,grid_space_x,grid_space_y,grid_space_z);
        string sufW(suffix); 
        std::cout<<sufW<<std::endl;

        std::cout<<"level "<<level+1<<" start: "<<sufW<<std::endl;
        std::cout<<"dilateRadius = "<<dilateRadius<<std::endl;
        std::cout<<"sigmaMin = "<<sigmaMin<<",  sigmaMax = "<<sigmaMax<<",  nScales = "<<nScales<<std::endl;
        fprintf(fp, "\nlevel %d start: %s\n", level+1, suffix);
        fprintf(fp, "dilateRadius = %d\n", dilateRadius);
        fprintf(fp, "sigmaMin = %f,  sigmaMax = %f,  nScales = %d\n", sigmaMin, sigmaMax, nScales);

        string maskImgName=bf+"mask"+sufW;
        string defImgName=bf+"Deformed12"+sufW;
        string defImgName_o=bf+"Deformed12_o"+sufW;
        string disp1ImgName=bf+"Disp12_x"+sufW;
        string disp2ImgName=bf+"Disp12_y"+sufW;
        string disp3ImgName=bf+"Disp12_z"+sufW;
        string JacImgName=bf+"Jac"+sufW;
        string HU1_filename=bf+"HU1"+sufW;
        string HU2_filename=bf+"HU2"+sufW;
        string VM1_filename=bf+"VM1"+sufW;
        string VM2_filename=bf+"VM2"+sufW;
        
        if (scale_z!=res)
        {
            if (res!=1000)
                fclose(fpStat);

            res_iter=0;

            res=scale_z;
            string tmpStr;
            stringstream tmpStr2;
            tmpStr2<<res;
            tmpStr=tmpStr2.str();
            //itoa( res, tmpStr, 10 );
            statFileName = string(buffer0)+"res"+tmpStr+".stat";
            fpStat=fopen(statFileName.c_str(),"w");

            fprintf(fpStat, "-----------------------------------------------------------------------------------------------------------------------\n" );
            fprintf(fpStat, "iter SSTVD-cost   weighted SSVMD-cost    weighted LAP-cost    total-cost   lmkerror  JacMin   JacMaxInv  SSVMD-weight\n");
            fprintf(fpStat, "-----------------------------------------------------------------------------------------------------------------------\n" );
        }
        else
            res_iter++;


        ///////////////////////////////////////////////
        //calculate current image resolution
        int xdim=(xdimo-1)/scale_x+1;
		int ydim=(ydimo-1)/scale_y+1;
		int zdim=(zdimo-1)/scale_z+1;
        int imgSize=xdim*ydim*zdim;
        //int planeSize=xdim*ydim;
        size[0]=xdim;
    	size[1]=ydim;
    	size[2]=zdim;
     
        float spacing_x = spacingFull_x*scale_x;
        float spacing_y = spacingFull_y*scale_y;
        float spacing_z = spacingFull_z*scale_z;
        spacing[0]=spacing_x;
        spacing[1]=spacing_y;
        spacing[2]=spacing_z;

        ///////////////////////////////////////////////
        //initialize coeff files
        int num_element_x=(xdim-1)/grid_space_x+4;
		int num_element_y=(ydim-1)/grid_space_y+4;
		int num_element_z=(zdim-1)/grid_space_z+4;
        int basisSize=num_element_x*num_element_y*num_element_z;
        //int coeffPlaneSize=num_element_x*num_element_y;
        float* c_x_for=new float[basisSize];
        float* c_y_for=new float[basisSize];
        float* c_z_for=new float[basisSize];
	    
        for (int k=0; k<basisSize; k++)
        {
            c_x_for[k]=0;
            c_y_for[k]=0;
            c_z_for[k]=0;
        }

        ///////////////////////////////////////////////
        float* image1HU=new float[imgSize];
    	float* image2HU=new float[imgSize];
        unsigned char* mask1=new unsigned char[imgSize];
        unsigned char* mask2=new unsigned char[imgSize];
        unsigned char* maskUnion=new unsigned char[imgSize];
    
        //used to output deformed image
        float* grad_x_T=new float[imgSize];
    	float* grad_y_T=new float[imgSize];
    	float* grad_z_T=new float[imgSize];
    
        /////////////////////////////////////////////
        //for vesselness images
        float* vm1=new float[imgSize];
    	float* vm2=new float[imgSize];
        float* vm_grad_x_T=new float[imgSize];
    	float* vm_grad_y_T=new float[imgSize];
    	float* vm_grad_z_T=new float[imgSize];

        ///////////////////////////////////////////////
        int precalcArraySize=64*grid_space_x*grid_space_y*grid_space_z;
        float* precalcu=new float[precalcArraySize];
    	float* precalcux=new float[precalcArraySize];
    	float* precalcuy=new float[precalcArraySize];
    	float* precalcuz=new float[precalcArraySize];
    	float* precalcu2x=new float[precalcArraySize];
    	float* precalcu2y=new float[precalcArraySize];
    	float* precalcu2z=new float[precalcArraySize];

        for (int z=0;z<4*grid_space_z;z++)
        {
			for (int y=0;y<4*grid_space_y;y++)
            {
				for (int x=0;x<4*grid_space_x;x++)
				{
					precalcu(x,y,z)=static_cast<float>(dcubic_bsp(static_cast<double>(x)/grid_space_x-2)*dcubic_bsp(static_cast<double>(y)/grid_space_y-2)*dcubic_bsp(static_cast<double>(z)/grid_space_z-2));
					precalcux(x,y,z)=static_cast<float>(dcubic_dbsp(static_cast<double>(x)/grid_space_x-2)*dcubic_bsp(static_cast<double>(y)/grid_space_y-2)*dcubic_bsp(static_cast<double>(z)/grid_space_z-2)/grid_space_x);
					precalcuy(x,y,z)=static_cast<float>(dcubic_bsp(static_cast<double>(x)/grid_space_x-2)*dcubic_dbsp(static_cast<double>(y)/grid_space_y-2)*dcubic_bsp(static_cast<double>(z)/grid_space_z-2)/grid_space_y);
					precalcuz(x,y,z)=static_cast<float>(dcubic_bsp(static_cast<double>(x)/grid_space_x-2)*dcubic_bsp(static_cast<double>(y)/grid_space_y-2)*dcubic_dbsp(static_cast<double>(z)/grid_space_z-2)/grid_space_z);
					precalcu2x(x,y,z)=static_cast<float>(dcubic_d2bsp(static_cast<double>(x)/grid_space_x-2)*dcubic_bsp(static_cast<double>(y)/grid_space_y-2)*dcubic_bsp(static_cast<double>(z)/grid_space_z-2)/grid_space_x/grid_space_x);
					precalcu2y(x,y,z)=static_cast<float>(dcubic_bsp(static_cast<double>(x)/grid_space_x-2)*dcubic_d2bsp(static_cast<double>(y)/grid_space_y-2)*dcubic_bsp(static_cast<double>(z)/grid_space_z-2)/grid_space_y/grid_space_y);
					precalcu2z(x,y,z)=static_cast<float>(dcubic_bsp(static_cast<double>(x)/grid_space_x-2)*dcubic_bsp(static_cast<double>(y)/grid_space_y-2)*dcubic_d2bsp(static_cast<double>(z)/grid_space_z-2)/grid_space_z/grid_space_z);
				}
            }
        }
        ///////////////////////////////////////////////
		//interpolate image to the current resolution
        downsizeImg_byAvg_OuterZero(image1HU, image1HUo, xdimo, ydimo, zdimo, scale_x, scale_y, scale_z, 1, bgValue);
		downsizeImg_byAvg_OuterZero(image2HU, image2HUo, xdimo, ydimo, zdimo, scale_x, scale_y, scale_z, 1, bgValue);
        //ImageType::Pointer iimage1c = DownsampleImage< ImageType >( iimage1, scale_x, scale_y, scale_z, linearInterpolation, bgValue);
        //ImageType::Pointer iimage2c = DownsampleImage< ImageType >( iimage2, scale_x, scale_y, scale_z, linearInterpolation, bgValue);

        downsizeImg_byAvg_OuterZero(vm1, vm1o, xdimo, ydimo, zdimo, scale_x, scale_y, scale_z, 1, 0);
		downsizeImg_byAvg_OuterZero(vm2, vm2o, xdimo, ydimo, zdimo, scale_x, scale_y, scale_z, 1, 0);

        rescaleImageIntensity( vm1, 0, 1, xdim, ydim, zdim ); 
        rescaleImageIntensity( vm2, 0, 1, xdim, ydim, zdim ); 

        ///////////////////////////////////////////////
        //interpolate u_existing_fullsize to u_existing at current resolution
		downsizeImg_byAvg_bdry(u_x_for_existing, u_x_for_existing_fullsize, xdimo, ydimo, zdimo, scale_x, scale_y, scale_z, scale_x);
		downsizeImg_byAvg_bdry(u_y_for_existing, u_y_for_existing_fullsize, xdimo, ydimo, zdimo, scale_x, scale_y, scale_z, scale_y);
		downsizeImg_byAvg_bdry(u_z_for_existing, u_z_for_existing_fullsize, xdimo, ydimo, zdimo, scale_x, scale_y, scale_z, scale_z);

        //calcJac(Je, u_x_for_existing, u_y_for_existing, u_z_for_existing, xdim, ydim, zdim);

        ///////////////////////////////////////////////
        
        deformed2_inner_mask( mask1o_def, mask1o, xdimo, ydimo, zdimo, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize, 0);

        /*MaskImageType::Pointer mask1oDef = ImageIO.initImage< MaskImageType >( iStart, sizeFull, direction, origin, spacingFull );
        ArrayUCHAR2ITKImage< MaskImageType >( mask1o_def, mask1oDef );
        ImageIO.WriteImage< MaskImageType >( mask1oDef, maskImgName );
        mask1oDef = NULL;*/

        union_masks( mask1o_def, mask2o, maskUniono, imgSizeFull );
        
        //interpolate the maskunion to the current resolution
        int maskCT=0;
        downsizeMask_wCt(maskUnion, maskUniono, xdimo, ydimo, zdimo, scale_x, scale_y, scale_z, maskCT);
        std::cout<<"maskCT = "<<maskCT<<std::endl;
        //wt_tv=1e-6;
        std::cout<<"wt_tv = "<<wt_tv<<std::endl;
        
        ///////////////////////////////////////////////
        //dilate mask (ROI) for current level
        std::cout<<"Dilating Mask ..."<<std::endl;
        
        MaskImageType::Pointer maskUc = ImageIO.initImage< MaskImageType >( iStart, size, direction, origin, spacing );
        ArrayUCHAR2ITKImage< MaskImageType >( maskUnion, maskUc );

        MaskImageType::Pointer maskUcd = DilateImage< MaskImageType >( maskUc, dilateRadius, dilateValue );
        maskUc = NULL;

        if (writeChoice==1)
            ImageIO.WriteImage< MaskImageType >( maskUcd, maskImgName );
        
        ITKImage2ArrayUCHAR< MaskImageType >( maskUcd, maskUnion );
        maskUcd = NULL;

        int maskCT1 = findMaskCT( maskUnion, imgSize );
        std::cout<<"maskCT1 = "<<maskCT1<<std::endl; //after dilating

        ///////////////////////////////////////////////
        //expand mask 1 voxel and initiate arrayBasis
        unsigned char* maskUnionExpand = new unsigned char[imgSize];
        int maskCT2=0;
        expand_mask(maskUnionExpand, maskUnion, xdim, ydim, zdim, maskCT2);
        std::cout<<"maskCT2 = "<<maskCT2<<std::endl; //after expanding by 1 outer surface

        fprintf(fp, "maskCT = %d,  maskCT1 = %d,  maskCT2 = %d\n", maskCT, maskCT1, maskCT2);

        ///////////////////////////////////////////////
        std::cout<<"initilizing arrayBasis"<<std::endl;
        for (int kk=0; kk<arraySize; kk++)
        {
            arrayBasis0[kk]=0;
            arrayBasis1[kk]=0;
            arrayBasis2[kk]=0;
        }
        CreateBlut2(arrayBasis0, arrayBasis1, arrayBasis2, xgrid, grid_space_x);


        std::cout<<"initilizing arrayLocation"<<std::endl;
        short int* arrayLocation = new short int[maskCT2*9];
        CreateLocation4(arrayLocation, u_x_for_existing, u_y_for_existing, u_z_for_existing, xdim, ydim, zdim, grid_space_x, grid_space_y, grid_space_z, xgrid, ygrid, zgrid, maskUnionExpand);


        ///////////////////////////////////////////////
        ///////////////////////////////////////////////
        //calculate derivatives of imagein1
        differentiate(image1HU, grad_x_T, grad_y_T, grad_z_T, xdim, ydim, zdim);
        ///////////////////////////////////////////////
        ///////////////////////////////////////////////



        ///////////////////////
        //write VM images
        //ITKImage2Array< ImageType >( VMImage1, vm1 );
        //ITKImage2Array< ImageType >( VMImage2, vm2 );

        if ( writeChoice == 1 )
        {
        ImageType::Pointer HUImage1 = ImageIO.initImage< ImageType >( iStart, size, direction, origin, spacing );
        Array2ITKImage< ImageType >( image1HU, HUImage1 );
        ImageIO.CastWriteImage< ImageType, HUImageType >( HUImage1, HU1_filename );

        ImageType::Pointer HUImage2 = ImageIO.initImage< ImageType >( iStart, size, direction, origin, spacing );
        Array2ITKImage< ImageType >( image2HU, HUImage2 );
        ImageIO.CastWriteImage< ImageType, HUImageType >( HUImage2, HU2_filename );

        ImageType::Pointer VMImage1 = ImageIO.initImage< ImageType >( iStart, size, direction, origin, spacing );
        Array2ITKImage< ImageType >( vm1, VMImage1 );
        ImageIO.WriteImage< ImageType >( VMImage1, VM1_filename );

        ImageType::Pointer VMImage2 = ImageIO.initImage< ImageType >( iStart, size, direction, origin, spacing );
        Array2ITKImage< ImageType >( vm2, VMImage2 );
        ImageIO.WriteImage< ImageType >( VMImage2, VM2_filename );

            /*std::cout<<"Writing VM images..."<<std::endl;
            string vm1ImgName=bf+"vm1"+sufW;
            string vm2ImgName=bf+"vm2"+sufW;
            ImageIO.WriteImage< ImageType >( VMImage1, vm1ImgName );
            ImageIO.WriteImage< ImageType >( VMImage2, vm2ImgName );*/

        HUImage1 = NULL;
        HUImage2 = NULL;
        VMImage1 = NULL;
        VMImage2 = NULL;
        }
    
        ///////////////////////////////////////////////
        ///////////////////////////////////////////////
        differentiate(vm1, vm_grad_x_T, vm_grad_y_T, vm_grad_z_T, xdim, ydim, zdim);
        ///////////////////////////////////////////////
        ///////////////////////////////////////////////


        ///////////////////////////////////////////////////////////////////////////
        //interpolate lmk at current resolution
        std::cout<<"interpolate lmk disp"<<std::endl;
	    for (int l=0; l<num_lmk; l++)
	    {
            xc_1[l]=static_cast<int>((static_cast<float>(xc_1o[l]))/scale_x+0.5);
	    	yc_1[l]=static_cast<int>(static_cast<float>(yc_1o[l])/scale_y+0.5);
	    	zc_1[l]=static_cast<int>(static_cast<float>(zc_1o[l])/scale_z+0.5);
	    	xc_2[l]=static_cast<int>((static_cast<float>(xc_2o[l]))/scale_x+0.5);
	    	yc_2[l]=static_cast<int>((static_cast<float>(yc_2o[l]))/scale_y+0.5);
	    	zc_2[l]=static_cast<int>((static_cast<float>(zc_2o[l]))/scale_z+0.5);
	    }
			
        ///////////////////////////////////////////////
		//optimize for sim
        if (true)
        {
            int basisSize1=num_element_x*num_element_y*num_element_z;
            if (basisSize1 != basisSize)
            {
                std::cout<<"basisSize1 not equal to basisSize!"<<std::endl;
                exit(3);
            }
            ap::real_1d_array phi;
            phi.setbounds(1,3*basisSize1);
            ap::integer_1d_array nbd;
            nbd.setbounds(1,3*basisSize1);
            ap::real_1d_array ll;
            ll.setbounds(1,3*basisSize1);
            ap::real_1d_array uu;
            uu.setbounds(1,3*basisSize1);
            int index=0;
            for (int k=0;k<num_element_z;k++)
                for (int j=0;j<num_element_y;j++)
                    for (int i=0;i<num_element_x;i++)
                    {
                        phi(index+1)=c_x_for[index];
                        phi(index+basisSize1+1)=c_y_for[index];
                        phi(index+2*basisSize1+1)=c_z_for[index];
                        index++;
                    }
            for (int k=0;k<3*basisSize1;k++)
            {
                nbd(k+1)=2;
            }
            for (int k=0;k<basisSize1;k++)
            {
                uu(k+1)=static_cast<float>(grid_space_x)/2.479472335;
                //std::cout<<k<<" "<<uu(k+1)<<std::endl;
                uu(k+basisSize1+1)=static_cast<float>(grid_space_y)/2.479472335;
                uu(k+2*basisSize1+1)=static_cast<float>(grid_space_z)/2.479472335;
                ll(k+1)=-uu(k+1);
                ll(k+basisSize1+1)=-uu(k+basisSize1+1);
                ll(k+2*basisSize1+1)=-uu(k+2*basisSize1+1);
            }
            //std::cout<<"ok here"<<std::endl;
            int info;
            std::cout<<"begin optimization"<<std::endl;
            printf("------------------------------------------------------------------------------------------------------------------------\n" );
            printf("iter SSTVD-cost   weighted SSVMD-cost    weighted LAP-cost    total-cost   lmkerror   JacMin   JacMaxInv  SSVMD-weight\n");
            printf("------------------------------------------------------------------------------------------------------------------------\n" );
            fprintf(fp, "-----------------------------------------------------------------------------------------------------------------------\n" );
            fprintf(fp, "iter SSTVD-cost   weighted SSVMD-cost    weighted LAP-cost    total-cost   lmkerror  JacMin   JacMaxInv  SSVMD-weight\n");
            fprintf(fp, "-----------------------------------------------------------------------------------------------------------------------\n" );


            lbfgsbminimize(3*basisSize1, 7, phi, eps1, eps2, eps3, maxIter, nbd, ll, uu, info, arrayBasis0, arrayBasis1, arrayBasis2, arrayLocation, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, grid_space_x, grid_space_y, grid_space_z, image1HU, image2HU, grad_x_T, grad_y_T, grad_z_T, c_x_for, c_y_for, c_z_for, u_x_for_existing, u_y_for_existing, u_z_for_existing, u_x_for_total, u_y_for_total, u_z_for_total, maskUnion, maskCT2, num_lmk, xc_1, yc_1, zc_1, xc_2, yc_2, zc_2, spacing_x, spacing_y, spacing_z, wt_tv, wt_lmk, wt_smooth, alpha, beta, gamma, fp, vm1, vm2, vm_grad_x_T, vm_grad_y_T, vm_grad_z_T, wt_vm, precalcu, precalcux, precalcuy, precalcuz, precalcu2x, precalcu2y, precalcu2z, fpStat, res_iter, cost_sstvd, cost_ved, cost_smooth, cost_lmk, JacMin, JacMax);
            std::cout<<"optimization ended "<<info<<std::endl;
            index=0;
            for (int k=0;k<num_element_z;k++)
                for (int j=0;j<num_element_y;j++)
                    for (int i=0;i<num_element_x;i++)
                    {
                        c_x_for[index]=phi(index+1);
                        c_y_for[index]=phi(index+basisSize1+1);
                        c_z_for[index]=phi(index+2*basisSize1+1);
                        index++;
                    }
        }

        fwrite(&level,sizeof(int),1,fpc);
        fwrite(&grid_space_x,sizeof(int),1,fpc);
        fwrite(&grid_space_y,sizeof(int),1,fpc);
        fwrite(&grid_space_z,sizeof(int),1,fpc);
        fwrite(&scale_x,sizeof(int),1,fpc);
        fwrite(&scale_y,sizeof(int),1,fpc);
        fwrite(&scale_z,sizeof(int),1,fpc);
        fwrite(&num_element_x,sizeof(int),1,fpc);
        fwrite(&num_element_y,sizeof(int),1,fpc);
        fwrite(&num_element_z,sizeof(int),1,fpc);
        fwrite(&xdim,sizeof(int),1,fpc);
        fwrite(&ydim,sizeof(int),1,fpc);
        fwrite(&zdim,sizeof(int),1,fpc);
        fwrite(c_x_for,sizeof(float),basisSize,fpc);
        fwrite(c_y_for,sizeof(float),basisSize,fpc);
        fwrite(c_z_for,sizeof(float),basisSize,fpc);

        ////////////////////////////////////////////////////////////////////////
        //update Jacobian
        BSplineCoeff2FullDispImage( c_x_for, c_y_for, c_z_for, u_x_for, u_y_for, u_z_for, grid_space_x, grid_space_y, grid_space_z, scale_x, scale_y, scale_z, num_element_x, num_element_y, num_element_z, xdimo, ydimo, zdimo );
        
        if (level==0)
        {
            calcJac(Jac_exis, u_x_for, u_y_for, u_z_for, xdimo, ydimo, zdimo);
        }
        else 
        {
            calcJac(Jac_curr, u_x_for, u_y_for, u_z_for, xdimo, ydimo, zdimo);
            deformed2_inner(Jac_curr_def, Jac_curr, xdimo, ydimo, zdimo, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize, 1 );
            composeJac(Jac_exis_updated, Jac_exis, Jac_curr_def, imgSizeFull);
            for (int i=0; i<imgSizeFull; i++)
                Jac_exis[i]=Jac_exis_updated[i];
        }
    
        //if ( writeChoice == 1 )
        {
            ImageType::Pointer JacImg = ImageIO.initImage< ImageType >( iStart, sizeFull, direction, origin, spacingFull );
            Array2ITKImage< ImageType >( Jac_exis, JacImg );
            ImageIO.WriteImage< ImageType >( JacImg, JacImgName );
        }

        ////////////////////////////////////////////////////////////////////////
        //update disps
        std::cout<<"updating displacements"<<std::endl;
        //interpolate c to current resolution         
        compose_u9(u_x_for_total, u_y_for_total, u_z_for_total, u_x_for_existing, u_y_for_existing, u_z_for_existing, arrayBasis0, arrayLocation, c_x_for, c_y_for, c_z_for, xdim, ydim, zdim, num_element_x, num_element_y, num_element_z, maskCT2);
        for (int k=0;k<imgSize;k++)
        {
            u_x_for_existing[k]=u_x_for_total[k];
            u_y_for_existing[k]=u_y_for_total[k];
            u_z_for_existing[k]=u_z_for_total[k];
        }
        
        //interpolate c to fullsize u_existing
        if (scale_x==1 && scale_y==1 && scale_z==1)
        {
            for (int k=0;k<xdim*ydim*zdim;k++)
            {
                u_x_for_existing_fullsize[k]=u_x_for_total[k];
                u_y_for_existing_fullsize[k]=u_y_for_total[k];
                u_z_for_existing_fullsize[k]=u_z_for_total[k];
            }
        }
        else
        {
            //BSplineCoeff2FullDispImage( c_x_for, c_y_for, c_z_for, u_x_for, u_y_for, u_z_for, grid_space_x, grid_space_y, grid_space_z, scale_x, scale_y, scale_z, num_element_x, num_element_y, num_element_z, xdimo, ydimo, zdimo );

            //interpolate c to fullsize u_existing
            compose_disp(u_x_for, u_y_for, u_z_for, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize, u_x_for_total, u_y_for_total, u_z_for_total, xdimo, ydimo, zdimo);

            for (int k=0;k<imgSizeFull;k++)
            {
                u_x_for_existing_fullsize[k]=u_x_for_total[k];
                u_y_for_existing_fullsize[k]=u_y_for_total[k];
                u_z_for_existing_fullsize[k]=u_z_for_total[k];
            }
        }

        ////////////////////////////////////////////////////////////////////////
        //calc lmkerror
        float mean12 = CalcLmkErrorAfter( num_lmk, xc_1o, yc_1o, zc_1o, xc_2o, yc_2o, zc_2o, xdimo, ydimo, zdimo, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize, spacingFull_x, spacingFull_y, spacingFull_z, dist12 );
        
        std::cout<<"Avg lmkerror after registration:"<<std::endl;
        printf("%2.5f\n",mean12);
        fprintf(fp,"Avg lmkerror after registration:\n%2.5f\n",mean12);
		
        ////////////////////////////////////////////////////////////////
        //generate deformed image
        std::cout<<"calculating and saving deformed image"<<std::endl;
        if (savefull == 1)
        {
            deformed2_inner(vdeformed1, image1HUo_o, xdimo, ydimo, zdimo, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize, bgValue);
            ImageType::Pointer idefImg = ImageIO.initImage< ImageType >( iStart, sizeFull, direction, origin, spacingFull );
            Array2ITKImage< ImageType >( vdeformed1, idefImg );

            ImageIO.CastWriteImage< ImageType, HUImageType >( idefImg, defImgName );
            idefImg = NULL;
        }
        else
        {
            deformed2_inner(vdeformed1, image1HU, xdim, ydim, zdim, u_x_for_existing, u_y_for_existing, u_z_for_existing, bgValue);
            ImageType::Pointer idefImg = ImageIO.initImage< ImageType >( iStart, size, direction, origin, spacing );
            Array2ITKImage< ImageType >( vdeformed1, idefImg );

            ImageIO.CastWriteImage< ImageType, HUImageType >( idefImg, defImgName );
            idefImg = NULL;
        }

        ///////////////////////////////////////////////////////////////
        //write deformed vm image
        //if (writeChoice==1 && wt_vm>0)
        if (wt_vm>0)
        {
            float* vm1def=new float[imgSize];
            deformed2_inner(vm1def, vm1, xdim, ydim, zdim, u_x_for_existing, u_y_for_existing, u_z_for_existing, 0);
            ImageType::Pointer VMImage1def = ImageIO.initImage< ImageType >( iStart, size, direction, origin, spacing );
            Array2ITKImage< ImageType >( vm1def, VMImage1def );
            string vm1defImgName=bf+"vmDef12"+sufW;
            ImageIO.WriteImage< ImageType >( VMImage1def, vm1defImgName );

            VMImage1def = NULL;
            delete [] vm1def;
        }

        ///////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////
        //write out disp images
        if (if_write_disp==1)
        {
            std::cout<<"Writing deformed _ORIGINAL_ template image"<<std::endl;
            deformed2_inner(vdeformed1, image1HUo_o, xdimo, ydimo, zdimo, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize, bgValue);
            ImageType::Pointer idefImg = ImageIO.initImage< ImageType >( iStart, sizeFull, direction, origin, spacingFull );
            Array2ITKImage< ImageType >( vdeformed1, idefImg );

            ImageIO.CastWriteImage< ImageType, HUImageType >( idefImg, defImgName_o );
            idefImg = NULL;


            std::cout<<"Writing displacements"<<std::endl;

            DispImageType::Pointer idisp1Img = ImageIO.initImage< DispImageType >( iStart, sizeFull, direction, origin, spacingFull );
            DispImageType::Pointer idisp2Img = ImageIO.initImage< DispImageType >( iStart, sizeFull, direction, origin, spacingFull );
            DispImageType::Pointer idisp3Img = ImageIO.initImage< DispImageType >( iStart, sizeFull, direction, origin, spacingFull );
            
            //Array2ITKImage< DispImageType >( u_x_for_existing_fullsize, idisp1Img );
            //Array2ITKImage< DispImageType >( u_y_for_existing_fullsize, idisp2Img );
            //Array2ITKImage< DispImageType >( u_z_for_existing_fullsize, idisp3Img );

            DispArray2ITKImage< DispImageType >( u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize, idisp1Img, idisp2Img, idisp3Img );

            ImageIO.WriteImage< DispImageType >( idisp1Img, disp1ImgName );
            ImageIO.WriteImage< DispImageType >( idisp2Img, disp2ImgName );
            ImageIO.WriteImage< DispImageType >( idisp3Img, disp3ImgName );


            /*sprintf(suffix,"_%d_%d_%d_%d_%d_%d.mhd",scale_x,scale_y,scale_z,grid_space_x,grid_space_y,grid_space_z);
            string sufW_mhd(suffix);
            string mhdName=bf+"Disp"+sufW_mhd;
            DispFieldType::Pointer mhdImg = ImageIO.initmhd< DispFieldType >( iStart, sizeFull, direction, origin, spacingFull );
            ComposeDispField< DispFieldType, ImageType >( mhdImg, idisp1Img, idisp2Img, idisp3Img );
            ImageIO.WriteImage< DispFieldType >( mhdImg, mhdName );
            mhdImg = NULL;*/

            idisp1Img = NULL;
            idisp2Img = NULL;
            idisp3Img = NULL;
        }

		///////////////////////////////////////////////////////////////////////////////////////////////////////
        //time report
        
        delete [] c_x_for;
        delete [] c_y_for;
        delete [] c_z_for;

        delete [] image1HU;
        delete [] image2HU;
        delete [] mask1;
        delete [] mask2;
        delete [] maskUnion;
        delete [] grad_x_T;
        delete [] grad_y_T;
        delete [] grad_z_T;
        delete [] vm1;
        delete [] vm2;
        delete [] vm_grad_x_T;
        delete [] vm_grad_y_T;
        delete [] vm_grad_z_T;

        delete [] maskUnionExpand;
        
        delete [] arrayLocation;

        delete [] precalcu;
        delete [] precalcux;
        delete [] precalcuy;
        delete [] precalcuz;
        delete [] precalcu2x;
        delete [] precalcu2y;
        delete [] precalcu2z;

        time(&end);
		tdiff[level]=difftime(end,start);
        std::cout<<"end of level "<<level+1<<" @ "<<tdiff[level]<<" seconds"<<std::endl;
        fprintf(fp, "end of level %d @ %d seconds\n", level+1, tdiff[level]);
        fprintf(fp, "%10.2f  %10.2f  %10.2f  %10.4f  %10.4f  %10.4f  %10.4f\n", cost_sstvd, cost_ved, cost_smooth, cost_lmk, JacMin, 1.0/JacMax, tdiff[level]/60.0);

        total_time+=tdiff[level];

	}

    
    for (int level=0; level<total_level; level++)
    {
        std::cout<<"Level "<<level+1<<" costs "<<tdiff[level]/60.0<<" min;"<<std::endl;
        fprintf(fp, "Level %d costs %6.2f min;\n", level+1, tdiff[level]/60.0);
    }
    std::cout<<"The total running time is: "<<total_time/60.0<<"min ("<<total_time/3600.0<<" hour)."<<std::endl;
    fprintf(fp,"The total running time is: %6.2f min (%4.2f hour).\n", total_time/60.0, total_time/3600.0);
    fprintf(fp, "\n%10.2f  %10.2f  %10.2f  %10.4f  %10.4f  %10.4f  %10.4f", cost_sstvd, cost_ved, cost_smooth, cost_lmk, JacMin, 1.0/JacMax, total_time/60.0);
	fclose(fpr);
    fclose(fp);
	fclose(fpc);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //final lmkerror report

    FILE *fp2;
    strcpy(buffer,dirName);
    strcat(buffer,fileNameLMKResult);
    fp2=fopen(buffer,"w");

    float mean12 = CalcLmkErrorAfter( num_lmk, xc_1o, yc_1o, zc_1o, xc_2o, yc_2o, zc_2o, xdimo, ydimo, zdimo, u_x_for_existing_fullsize, u_y_for_existing_fullsize, u_z_for_existing_fullsize, spacingFull_x, spacingFull_y, spacingFull_z, dist12 );

    fprintf(fp2,"lmkerror after registration:\n");
    for ( int l=0; l<num_lmk; l++ )
    {
       fprintf(fp2,"%2.5f\t%2.5f\n",dist[l],dist12[l]);
    }

    std::cout<<"Avg lmkerror after registration:"<<std::endl;
    printf("%2.5f\t%2.5f\n",avglmkdist,mean12);
    fprintf(fp2,"Avg lmkerror after registration:\n");
    fprintf(fp2,"%2.5f\t%2.5f\n",avglmkdist,mean12);
    fclose(fp2);


    ///////////////////////////////
    //free space
    delete [] image1HUo_o;
    delete [] image1HUo;
    delete [] image2HUo;

    delete [] mask1o;
    delete [] mask2o;
    delete [] maskUniono;
    delete [] vdeformed1;

    delete [] vm1o;
    delete [] vm2o;

    delete [] u_x_for_total;
    delete [] u_y_for_total;
    delete [] u_z_for_total;

    delete [] u_x_for;
    delete [] u_y_for;
    delete [] u_z_for;

    delete [] u_x_for_existing_fullsize;
    delete [] u_y_for_existing_fullsize;
    delete [] u_z_for_existing_fullsize;

    delete [] u_x_for_existing;
    delete [] u_y_for_existing;
    delete [] u_z_for_existing;

    delete [] tdiff;

    delete [] Jac_curr;
    delete [] Jac_curr_def;
    delete [] Jac_exis;
    delete [] Jac_exis_updated;
    //delete [] Je;

    delete [] arrayBasis0;
    delete [] arrayBasis1;
    delete [] arrayBasis2;

    delete [] xc_1o;
    delete [] yc_1o;
    delete [] zc_1o;
    delete [] xc_2o;
    delete [] yc_2o;
    delete [] zc_2o;

    delete [] dist;
    delete [] dist12;

    delete [] xc_1;
    delete [] yc_1;
    delete [] zc_1;
    delete [] xc_2;
    delete [] yc_2;
    delete [] zc_2;
    ///////////////////////////////

    return 0;
}


