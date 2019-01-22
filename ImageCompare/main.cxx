#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

//The program check if 2 images are the same and print the file name if so 

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile1 inputImageFile2 " << std::endl;
      return EXIT_FAILURE;
    }

  //ITK settings
  const unsigned int Dimension = 2;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::AbsoluteValueDifferenceImageFilter< ImageType, ImageType, ImageType > AbsDiffFilterType;
  typedef itk::StatisticsImageFilter< ImageType > StatFilterType;

  //Filters
  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  AbsDiffFilterType::Pointer absDiffFilter = AbsDiffFilterType::New();
  StatFilterType::Pointer statFilter = StatFilterType::New();

  //Parameters
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );

  //Pipeline
  try
    {
      reader1->Update();
      reader2->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  absDiffFilter->SetInput1(reader1->GetOutput());
  absDiffFilter->SetInput2(reader2->GetOutput());

  statFilter->SetInput(absDiffFilter->GetOutput());
  statFilter->Update();

  double total = statFilter->GetSum();
  if(total<1)
    {
      char whole1[50];
      strcpy(whole1, argv[1]);
      char folder1[10];
      char file1[10];
      int i = 0;
      while(whole1[i]!='/')
	{
	  folder1[i] = whole1[i];
	  i++;
	}
      folder1[i+1] = '\0';
      int j = i+1;
      while(whole1[j]!='.')
	{
	  file1[j-i-1] = whole1[j];
	  j++;
	}
      file1[j-i-1] = '\0';  
      //std::cout<<folder1<<"; "<<file1<<"."<<std::endl;

      char whole2[50];
      strcpy(whole2, argv[2]);
      char folder2[10];
      char file2[10];
      int a = 0;
      while(whole2[a]!='/')
	{
	  folder2[a] = whole2[a];
	  a++;
	}
      folder2[a+1] = '\0';
      int b = a+1;
      while(whole2[b]!='.')
	{
	  file2[b-a-1] = whole2[b];
	  b++;
	}
      file2[b-a-1] = '\0';  
      //std::cout<<folder2<<"; "<<file2<<"."<<std::endl;

      std::cout<<"mv "<<folder2<<"/"<<file2<<".png "<<folder2<<"/Ori/"<<file1<<".png"<<std::endl;
      std::cout<<"mv "<<folder2<<"/"<<file2<<"_bin.png "<<folder2<<"/Ori/"<<file1<<"_bin.png"<<std::endl;
    }
  return EXIT_SUCCESS;
}
