/*
 * =====================================================================================
 *
 *       Filename:  ConfigReader.cxx
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  08/30/2011 09:00:01 PM CDT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */

#include "ConfigReader.h"

using namespace std;

bool Configuration::ReadConfiguration( char const* filename )
{
    ifstream infile( filename );
    if ( infile.good() == false ){
        return false;
    }

    string cmd_str;
    string in_str;
    string line;


    //////////////////////////////////////////////////////
    //reading global settings
    
    ////////////////////////////
    //set as default
    if_read_disp=0;
    mask_lower_thr=1;
    mask_upper_thr=255;
    wt_tv=1;
    wt_vm=0;
    wt_lap=0;
    wt_smooth=0;
    alpha=0;
    beta=0;
    gamma=0;
    wt_lmk=0;

    readVM=0;
    ////////////////////////////

    infile >> cmd_str;
    //while ( infile.eof() == false ){
    while ( cmd_str != "[Global_Setting_End]" ){

        std::cout<<cmd_str<<" ";
        if ( cmd_str == "[If_Read_Disp]" ){
            infile >> if_read_disp;
            std::cout<<if_read_disp<<std::endl;
        }
        if ( cmd_str == "[Initial_Disp_Dir]" ){
            infile >> initial_disp_dir;
            std::cout<<initial_disp_dir<<std::endl;
        }
        if ( cmd_str == "[Initial_Disp_Suffix]" ){
            infile >> initial_disp_suffix;
            std::cout<<initial_disp_suffix<<std::endl;
        }
        if ( cmd_str == "[Moving_Image_File]" ){
            infile >> moving_image_file;
            std::cout<<moving_image_file<<std::endl;
        }
        if ( cmd_str == "[Fixed_Image_File]" ){
            infile >> fixed_image_file;
            std::cout<<fixed_image_file<<std::endl;
        }
        if ( cmd_str == "[Moving_Mask_File]" ){
            infile >> moving_mask_file;
            std::cout<<moving_mask_file<<std::endl;
        }
        if ( cmd_str == "[Fixed_Mask_File]" ){
            infile >> fixed_mask_file;
            std::cout<<fixed_mask_file<<std::endl;
        }
        if ( cmd_str == "[Moving_Lmk_File]" ){
            infile >> moving_lmk_file;
            std::cout<<moving_lmk_file<<std::endl;
        }
        if ( cmd_str == "[Fixed_Lmk_File]" ){
            infile >> fixed_lmk_file;
            std::cout<<fixed_lmk_file<<std::endl;
        }

        if ( cmd_str == "[Moving_VM_File]" ){
            infile >> moving_vm_file;
            std::cout<<moving_vm_file<<std::endl;
        }
        if ( cmd_str == "[Fixed_VM_File]" ){
            infile >> fixed_vm_file;
            std::cout<<fixed_vm_file<<std::endl;
        }

        if ( cmd_str == "[Result_Dir]" ){
            infile >> result_dir;
            std::cout<<result_dir<<std::endl;
        }
        if ( cmd_str == "[Log_Filename]" ){
            infile >> log_filename;
            std::cout<<log_filename<<std::endl;
        }
        if ( cmd_str == "[Coeff_Filename]" ){
            infile >> coeff_filename;
            std::cout<<coeff_filename<<std::endl;
        }
        if ( cmd_str == "[Lmkerror_Filename]" ){
            infile >> lmkerror_filename;
            std::cout<<lmkerror_filename<<std::endl;
        }

        if ( cmd_str == "[Total_Level_Number]" ){
            infile >> total_level_number;
            std::cout<<total_level_number<<std::endl;
        }
        if ( cmd_str == "[If_Write_Internal_Image]" ){
            infile >> if_write_internal_image;
            std::cout<<if_write_internal_image<<std::endl;
        }
        if ( cmd_str == "[If_Save_Fullsize_Disp]" ){
            infile >> if_save_fullsize_disp;
            std::cout<<if_save_fullsize_disp<<std::endl;
        }

        if ( cmd_str == "[If_Read_VM]" ){
            infile >> readVM;
            std::cout<<readVM<<std::endl;
        }

        if ( cmd_str == "[Weight_SSTVD]" ){
            infile >> wt_tv;
            std::cout<<wt_tv<<std::endl;
        }
        if ( cmd_str == "[Weight_SSVMD]" ){
            infile >> wt_vm;
            std::cout<<wt_vm<<std::endl;
        }
        if ( cmd_str == "[Weight_Laplacian]" ){
            infile >> wt_lap;
            std::cout<<wt_lap<<std::endl;
        }
        if ( cmd_str == "[Weight_Smooth]" ){
            infile >> wt_smooth;
            std::cout<<wt_smooth<<std::endl;
        }
        if ( cmd_str == "[Smooth_Alpha]" ){
            infile >> alpha;
            std::cout<<alpha<<std::endl;
        }
        if ( cmd_str == "[Smooth_Beta]" ){
            infile >> beta;
            std::cout<<beta<<std::endl;
        }
        if ( cmd_str == "[Smooth_Gamma]" ){
            infile >> gamma;
            std::cout<<gamma<<std::endl;
        }
        if ( cmd_str == "[Weight_Lmk]" ){
            infile >> wt_lmk;
            std::cout<<wt_lmk<<std::endl;
        }

        if ( cmd_str == "[Blut_Scale]" ){
            infile >> blut_scale;
            std::cout<<blut_scale<<std::endl;
        }
        if ( cmd_str == "[Mask_Lower_Threshold]" ){
            infile >> mask_lower_thr;
            std::cout<<mask_lower_thr<<std::endl;
        }
        if ( cmd_str == "[Mask_Upper_Threshold]" ){
            infile >> mask_upper_thr;
            std::cout<<mask_upper_thr<<std::endl;
        }

        if ( cmd_str == "[VM_Sigma_Min_Global]" ){
            infile >> vm_sigma_min_global;
            std::cout<<vm_sigma_min_global<<std::endl;
        }
        if ( cmd_str == "[VM_Sigma_Max_Global]" ){
            infile >> vm_sigma_max_global;
            std::cout<<vm_sigma_max_global<<std::endl;
        }
        if ( cmd_str == "[VM_Iteration_Global]" ){
            infile >> vm_iteration_global;
            std::cout<<vm_iteration_global<<std::endl;
        }


        infile >> cmd_str;

    }


    //////////////////////////////////////////////////////
    //reading settings for each level
   
    LevelConfiguration example_level;

    ////////////////////////////
    //set as default
    example_level.dilate_radius=0;
    example_level.if_write_disp=0;
    for (int j=0; j<3; j++)
        example_level.epsilon.push_back(0.0);
    example_level.epsilon[1]=0.05;
    ////////////////////////////

    example_level.level_number=0;
    example_level.iteration=0;
        
    example_level.vm_sigma_min=0;
    example_level.vm_sigma_max=0;
    example_level.vm_iteration=0;

    for (int j=0; j<3; j++)
        example_level.image_size_ratio.push_back(1);
    for (int j=0; j<3; j++)
        example_level.spline_spacing.push_back(4);

    
    
    


    //vector< LevelConfiguration > level_setting( total_level_number, example_level );
    for ( int i=0; i<total_level_number; i++ )
    {
        level_setting.push_back( example_level );
        std::cout<<std::endl;

        infile >> cmd_str;
        while ( cmd_str != "[Level_Setting_End]" ){

            std::cout<<cmd_str<<" ";
            if ( cmd_str == "[Level_Number]" ){
                infile >> level_setting[i].level_number;
                std::cout<<level_setting[i].level_number<<std::endl;
            }
            if ( cmd_str == "[Iteration]" ){
                infile >> level_setting[i].iteration;
                std::cout<<level_setting[i].iteration<<std::endl;
            }
            if ( cmd_str == "[VM_Sigma_Min]" ){
                infile >> level_setting[i].vm_sigma_min;
                std::cout<<level_setting[i].vm_sigma_min<<std::endl;
            }
            if ( cmd_str == "[VM_Sigma_Max]" ){
                infile >> level_setting[i].vm_sigma_max;
                std::cout<<level_setting[i].vm_sigma_max<<std::endl;
            }
            if ( cmd_str == "[VM_Iteration]" ){
                infile >> level_setting[i].vm_iteration;
                std::cout<<level_setting[i].vm_iteration<<std::endl;
            }
            if ( cmd_str == "[Dilate_Radius]" ){
                infile >> level_setting[i].dilate_radius;
                std::cout<<level_setting[i].dilate_radius<<std::endl;
            }
            if ( cmd_str == "[If_Write_Disp]" ){
                infile >> level_setting[i].if_write_disp;
                std::cout<<level_setting[i].if_write_disp<<std::endl;
            }
            if ( cmd_str == "[Epsilon]" ){
                getline( infile, line );
                stringstream lineStream( line );
                string tmpString;
                getline(lineStream, tmpString, ' ');//skip the first blank space
                int k=0;
                while (getline(lineStream, tmpString, ' '))
                {
                    istringstream convertStream(tmpString);
                    float tmpVal;
                    convertStream >> tmpVal;
                    level_setting[i].epsilon[k]=tmpVal;
                    k++;
                }
                //std::cout<<k<<std::endl;

                for ( int j=0; j<3; j++ )
                    std::cout<<level_setting[i].epsilon[j]<<" ";
                std::cout<<std::endl;
            }
            if ( cmd_str == "[Image_Size_Ratio]" ){
                getline( infile, line );
                stringstream lineStream( line );
                string tmpString;
                getline(lineStream, tmpString, ' ');//skip the first blank space
                int k=0;
                while (getline(lineStream, tmpString, ' '))
                {
                    istringstream convertStream(tmpString);
                    int tmpVal;
                    convertStream >> tmpVal;
                    level_setting[i].image_size_ratio[k]=tmpVal;
                    k++;
                }
                //std::cout<<k<<std::endl;

                for ( int j=0; j<3; j++ )
                    std::cout<<level_setting[i].image_size_ratio[j]<<" ";
                std::cout<<std::endl;
            }
            if ( cmd_str == "[Spline_Spacing]" ){
                getline( infile, line );
                stringstream lineStream( line );
                string tmpString;
                getline(lineStream, tmpString, ' ');//skip the first blank space
                int k=0;
                while (getline(lineStream, tmpString, ' '))
                {
                    istringstream convertStream(tmpString);
                    int tmpVal;
                    convertStream >> tmpVal;
                    level_setting[i].spline_spacing[k]=tmpVal;
                    k++;
                }
                //std::cout<<k<<std::endl;

                for ( int j=0; j<3; j++ )
                    std::cout<<level_setting[i].spline_spacing[j]<<" ";
                std::cout<<std::endl;

                if ( level_setting[i].dilate_radius == 0 )
                {
                    int grid_space_x=level_setting[i].spline_spacing[0];
                    level_setting[i].dilate_radius=static_cast<int>(ceil(grid_space_x/2.479472335));
                }
            }
            

            infile >> cmd_str;
        }

        //level_setting.push_back( level );
    }


    return true;
}


