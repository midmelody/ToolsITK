/*
 * =====================================================================================
 *
 *       Filename:  ConfigReader.h
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  08/30/2011 03:32:46 PM CDT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */

#ifndef __ConfigReader_h
#define __ConfigReader_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
    
struct LevelConfiguration
{
    int level_number;
    int iteration;
    int dilate_radius;
    int if_write_disp;

    float vm_sigma_min;
    float vm_sigma_max;
    float vm_iteration;

    std::vector< int > image_size_ratio;
    std::vector< int > spline_spacing;
    std::vector< float > epsilon;

    /*int image_size_ratio_x;
    int image_size_ratio_y;
    int image_size_ratio_z;

    int spline_spacing_x;
    int spline_spacing_y;
    int spline_spacing_z;

    float epsilon_x;
    float epsilon_y;
    float epsilon_z;*/
    
};

class Configuration
{
    public:
    //global setting
    int if_read_disp;
    std::string initial_disp_dir;
    std::string initial_disp_suffix;

    std::string moving_image_file;
    std::string fixed_image_file;
    std::string moving_mask_file;
    std::string fixed_mask_file;
    std::string moving_lmk_file;
    std::string fixed_lmk_file;
    std::string moving_vm_file;
    std::string fixed_vm_file;

    int readVM;

    std::string result_dir;
    std::string log_filename;
    std::string coeff_filename;
    std::string lmkerror_filename;

    int total_level_number;
    int if_save_fullsize_disp;
    int if_write_internal_image;

    float wt_tv;
    float wt_vm;
    float wt_lap;
    float wt_smooth;
    float alpha;
    float beta;
    float gamma;
    float wt_lmk;
    int blut_scale;

    int mask_lower_thr;
    int mask_upper_thr;

        
    float vm_sigma_min_global;
    float vm_sigma_max_global;
    float vm_iteration_global;


    std::vector< LevelConfiguration > level_setting;


    public:
    bool ReadConfiguration( char const* filename );
};




#endif

