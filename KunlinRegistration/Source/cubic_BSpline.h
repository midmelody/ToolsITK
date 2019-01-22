/*
 * =====================================================================================
 *
 *       Filename:  cubic_BSpline.h
 *
 *    Description:  CURSOR>
 *
 *        Version:  1.0
 *        Created:  11/08/2009 07:32:25 PM CST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  first_name last_name (fl), fl@my-company.com
 *        Company:  my-company
 *
 * =====================================================================================
 */

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float quad_bsp ( float x )
{
    float y;
    float a;
    if ( x<-1.5 || x>1.5 )
    {
        y=0;
    }
    else
        if ( x<-0.5 )
        {
            a=x+1.5;
            y=a*a/2.0;
        }
        else
            if ( x<0.5 )
            {
                a=x+0.5;
                y=0.5+a-a*a;
            }
            else
            {
                a=x-0.5;
                y=0.5-a+a*a/2.0;
            }
    return y;
}

////////////////////////////////////////

//b(x)
float cubic_bsp( float x )
{
    float y;
    float a;
    if ( x<-2 || x>2 )
    {
        y=0;
    }
    else 
        if ( x<-1 )
        {
            a=x+2;
            y=a*a*a/6.0;
        }
        else 
            if ( x<0 )
            {
                a=x+1;
                y=1.0/6.0+a/2.0+a*a/2.0-a*a*a/2.0;
            }
            else
                if ( x<1 )
                {
                    y=2.0/3.0-x*x+x*x*x/2.0;
                }
                else
                {
                    a=x-1;
                    y=1.0/6.0-a/2.0+a*a/2.0-a*a*a/6.0;
                }
    return y;
}

double dcubic_bsp( double x )
{
    double y;
    double a;
    if ( x<-2 || x>2 )
    {
        y=0;
    }
    else 
        if ( x<-1 )
        {
            a=x+2;
            y=a*a*a/6.0;
        }
        else 
            if ( x<0 )
            {
                a=x+1;
                y=1.0/6.0+a/2.0+a*a/2.0-a*a*a/2.0;
            }
            else
                if ( x<1 )
                {
                    y=2.0/3.0-x*x+x*x*x/2.0;
                }
                else
                {
                    a=x-1;
                    y=1.0/6.0-a/2.0+a*a/2.0-a*a*a/6.0;
                }
    return y;
}

//first derivative of b(x)
float cubic_dbsp( float x )
{
    float y;
    float a;
    if ( x<-2 || x>2 )
    {
        y=0;
    }
    else 
        if ( x<-1 )
        {
            a=x+2;
            y=a*a/2.0;
        }
        else 
            if ( x<0 )
            {
                a=x+1;
                y=1.0/2.0+a-a*a*3.0/2.0;
            }
            else
                if ( x<1 )
                {
                    y=-2.0*x+x*x*3.0/2.0;
                }
                else
                {
                    a=x-1;
                    y=-1.0/2.0+a-a*a/2.0;                
                }
    return y;
}

double dcubic_dbsp( double x )
{
    double y;
    double a;
    if ( x<-2 || x>2 )
    {
        y=0;
    }
    else 
        if ( x<-1 )
        {
            a=x+2;
            y=a*a/2.0;
        }
        else 
            if ( x<0 )
            {
                a=x+1;
                y=1.0/2.0+a-a*a*3.0/2.0;
            }
            else
                if ( x<1 )
                {
                    y=-2.0*x+x*x*3.0/2.0;
                }
                else
                {
                    a=x-1;
                    y=-1.0/2.0+a-a*a/2.0;                
                }
    return y;
}

//second derivative of b(x)
float cubic_d2bsp( float x )
{
    float y;
    float a;
    if ( x<-2 || x>2 )
    {
        y=0;
    }
    else 
        if ( x<-1 )
        {
            y=x+2;
        }
        else 
            if ( x<0 )
            {
                a=x+1;
                y=1.0-a*3.0;
            }
            else
                if ( x<1 )
                {
                    y=-2.0+x*3.0;
                }
                else
                {
                    a=x-1;
                    y=1.0-a;                
                }
    return y;
}

double dcubic_d2bsp( double x )
{
    double y;
    double a;
    if ( x<-2 || x>2 )
    {
        y=0;
    }
    else 
        if ( x<-1 )
        {
            y=x+2;
        }
        else 
            if ( x<0 )
            {
                a=x+1;
                y=1.0-a*3.0;
            }
            else
                if ( x<1 )
                {
                    y=-2.0+x*3.0;
                }
                else
                {
                    a=x-1;
                    y=1.0-a;                
                }
    return y;
}
