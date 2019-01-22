/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    OptionList.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "OptionList.h"
#include <cstdlib>

OptionList::OptionList(int argc, char* argv[])
{
  std::string tag ;

  int index = 1 ;
  while (index < argc)
    {
      if (argv[index][0] == '-' && argv[index][1] == '-')
        {
          // tag
          tag = argv[index] ;
          tag = tag.erase(0, 2) ; // remove '--' 
        }
      else
        {
          // value 
          const std::string value = argv[index] ;
          m_Map.insert(std::make_pair<const std::string, std::string>(tag, value)) ;
        }
      index++ ;
    }
}

int OptionList::GetOption(const std::string option_tag, StringVector* values)
{
  values->clear() ;
  typedef OptionMap::const_iterator CI ;
  int count = 0 ;

  const OptionMap & Const_m_Map=m_Map; //Need to force use of const version of equal_range
  std::pair<CI, CI> bound ( Const_m_Map.equal_range(option_tag) ) ;
  for (CI i = bound.first ; i != bound.second ; ++i)
    {
      values->push_back(i->second) ;
      count++ ;
    }
  return count ;
}


int OptionList::DumpOption(const std::string option_tag, bool withTag,
                           bool withNewLine)
{
  typedef OptionMap::const_iterator CI ;
  const OptionMap & Const_m_Map=m_Map; //Need to force use of const version of equal_range
  std::pair<CI, CI> bound ( Const_m_Map.equal_range(option_tag) ) ;
  int count = 0 ;

  if ( bound.first !=  bound.second)
    {
      if (withTag)
        std::cout << "--" << option_tag << " " ;
      
      for (CI i =  bound.first ; i !=  bound.second ; ++i)
        {
          std::cout << i->second << " " ;
          count++ ;
        }
      
      if (withNewLine)
        std::cout << std::endl ;
      
      return count++ ;
    }

    return 0 ;
}


int OptionList::GetMultiDoubleOption(const std::string tag,
                                     std::vector<double>* args, 
                                     bool required)
{
  args->clear() ;
  
  StringVector temp_args ;
  int arg_no = this->GetOption(tag, &temp_args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return -1 ;
  
  if (temp_args[0] == "-")
    return -2 ;
  
  for (int i = 0 ; i < arg_no ; i++)
    {
      args->push_back( atof(temp_args[i].c_str()) ) ;
    }
  
  return arg_no ;
}



int OptionList::GetMultiDoubleOption(const std::string tag,
                                     itk::Array<double>* args, 
                                     bool required)
{
  // use a temporary std::vector<> because the size is not 
  // known in advance
  std::vector<double> tmp;
  int arg_no =  GetMultiDoubleOption( tag, &tmp, required);
  if( arg_no <= 0 )
    {
    return arg_no;
    }
 
  itk::Array<double> array( arg_no );
  
  for (int i = 0 ; i < arg_no ; i++)
    {
    array[i] = tmp[i];  
    }
  
  *args = array;

  return arg_no ;
}


double OptionList::GetDoubleOption(const std::string tag, double default_value,
                                   bool required)
{
  StringVector temp_args ;
  int arg_no = this->GetOption(tag, &temp_args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return default_value ;
  
  return atof(temp_args[0].c_str()) ;
}

bool OptionList::GetBooleanOption(const std::string tag, bool default_value, bool required)
{
  StringVector args ;
  int arg_no = this->GetOption(tag, &args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return default_value ;
  
  if (args[0] == "yes")
    {
      return true ;
    }
  else
    {
      return false ;
    }
}

int OptionList::GetMultiIntOption(const std::string tag,
                                  std::vector<int>* args, 
                                  bool required)
{
  args->clear() ;
  
  StringVector temp_args ;
  int arg_no = this->GetOption(tag, &temp_args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return -1 ;
  
  if (temp_args[0] == "-")
    return -2 ;
  
  for (int i = 0 ; i < arg_no ; i++)
    {
      args->push_back( atoi(temp_args[i].c_str()) ) ;
    }
  
  return arg_no ;
}

int OptionList::GetMultiUIntOption(const std::string tag,
                                  std::vector<unsigned int>* args, 
                                  bool required)
{
  args->clear() ;
  
  StringVector temp_args ;
  int arg_no = this->GetOption(tag, &temp_args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return -1 ;
  
  if (temp_args[0] == "-")
    return -2 ;
  
  for (int i = 0 ; i < arg_no ; i++)
    {
      args->push_back( (unsigned int) atoi(temp_args[i].c_str()) ) ;
    }
  
  return arg_no ;
}

int OptionList::GetMultiUCharOption(const std::string tag,
                                    std::vector<unsigned char>* args, 
                                    bool required)
{
  args->clear() ;
  
  StringVector temp_args ;
  int arg_no = this->GetOption(tag, &temp_args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return -1 ;
  
  if (temp_args[0] == "-")
    return -2 ;
  
  for (int i = 0 ; i < arg_no ; i++)
    {
      args->push_back( (unsigned char) atoi(temp_args[i].c_str()) ) ;
    }
  
  return arg_no ;
}

int OptionList::GetIntOption(const std::string tag, int default_value, bool required)
{
  StringVector args ;
  int arg_no = this->GetOption(tag, &args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return default_value ;
  
  return atoi(args[0].c_str()) ;
}

unsigned int OptionList::GetUIntOption(const std::string tag,
                                       unsigned int default_value, 
                                       bool required)
{
  StringVector args ;
  int arg_no = this->GetOption(tag, &args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return default_value ;
  
  return (unsigned int) atoi(args[0].c_str()) ;
}

unsigned char OptionList::GetUCharOption(const std::string tag, unsigned char default_value, bool required)
{
  StringVector args ;
  int arg_no = this->GetOption(tag, &args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return default_value ;
  
  return (unsigned char) atoi(args[0].c_str()) ;
}

int OptionList::GetStringOption(const std::string tag,
                                std::string* ret, 
                                bool required)
{
  StringVector args ;
  int arg_no = this->GetOption(tag, &args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return -1 ;
  
  *ret = args[0] ;
  return arg_no ;
}

int OptionList::GetMultiStringOption(const std::string tag,
                                     std::vector< std::string >* ret, 
                                     bool required)
{
  ret->clear() ;
  StringVector args ;
  int arg_no = this->GetOption(tag, &args) ;
  
  if (required && arg_no == 0)
    throw RequiredOptionMissing(tag) ;
  
  if (arg_no == 0)
    return -1 ;
  
  for (int i = 0 ; i < arg_no ; i++)
    {
      ret->push_back( args[i] ) ;
    }

  return arg_no ;
}
