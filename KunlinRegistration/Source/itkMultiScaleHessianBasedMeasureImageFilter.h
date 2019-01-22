/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleHessianBasedMeasureImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007/04/01 23:13:46 $
  Version:   $Revision: 1.6 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiScaleHessianBasedMeasureImageFilter_h
#define __itkMultiScaleHessianBasedMeasureImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkHessianRecursiveGaussianImageFilter.h"

namespace itk
{
/**\class MultiScaleHessianBasedMeasureImageFilter
 * \brief A filter to enhance structures using Hessian eigensystem-based 
 * measures in a multiscale framework
 * 
 * The filter evaluates a Hessian-based enhancement measure, such as vesselness 
 * or objectness, at different scale levels. The Hessian-based measure is computed 
 * from the Hessian image at each scale level and the best response is selected. 
 *
 * Minimum and maximum sigma value can be set using SetMinSigma and SetMaxSigma
 * methods respectively. The number of scale levels is set using 
 * SetNumberOfSigmaSteps method. Exponentially distributed scale levels are 
 * computed within the bound set by the minimum and maximum sigma values 
 * 
 * The filter computes a second output image (accessed by the GetScalesOutput method)
 * containing the scales at which each pixel gave the best reponse. 
 *
 * \author Luca Antiga Ph.D.  Medical Imaging Unit,
 *                            Bioengineering Deparment, Mario Negri Institute, Italy.
 *
 * \sa HessianToObjectnessMeasureImageFilter 
 * \sa Hessian3DToVesselnessMeasureImageFilter 
 * \sa HessianSmoothed3DToVesselnessMeasureImageFilter 
 * \sa HessianRecursiveGaussianImageFilter 
 * \sa SymmetricEigenAnalysisImageFilter
 * \sa SymmetricSecondRankTensor
 * 
 * \ingroup IntensityImageFilters TensorObjects
 *
 */
template <class TInputImage,
          class THessianToMeasureFilter, 
          class TOutputImage=TInputImage >
class ITK_EXPORT MultiScaleHessianBasedMeasureImageFilter 
: public
ImageToImageFilter< TInputImage,TOutputImage > 
{
public:
  /** Standard class typedefs. */
  typedef MultiScaleHessianBasedMeasureImageFilter Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>              Superclass;

  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  typedef TInputImage                                    InputImageType;
  typedef TOutputImage                                   OutputImageType;

  typedef THessianToMeasureFilter                        HessianToMeasureFilterType;
  typedef typename HessianToMeasureFilterType::InputImageType    HessianImageType;

  typedef typename TInputImage::PixelType                InputPixelType;
  typedef typename TOutputImage::PixelType               OutputPixelType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);
 
  /** Hessian computation filter. */
  typedef HessianRecursiveGaussianImageFilter< InputImageType, HessianImageType> HessianFilterType;
 
  /** Update image buffer that holds the best objectness response */ 
  typedef Image< double, itkGetStaticConstMacro(ImageDimension) > UpdateBufferType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Set/Get macros for SigmaMin */
  itkSetMacro(SigmaMin, double);
  itkGetMacro(SigmaMin, double);
  
  /** Set/Get macros for SigmaMax */
  itkSetMacro(SigmaMax, double);
  itkGetMacro(SigmaMax, double);

  /** Set/Get macros for Number of Scales */
  itkSetMacro(NumberOfSigmaSteps, int);
  itkGetMacro(NumberOfSigmaSteps, int);

  typedef enum { EquispacedSigmaSteps = 0,
                 LogarithmicSigmaSteps = 1 } SigmaStepMethodType;

  /** Set/Get the method used to generate scale sequence (Equispaced or Logarithmic)*/
  itkSetMacro(SigmaStepMethod, SigmaStepMethodType);
  itkGetMacro(SigmaStepMethod, SigmaStepMethodType);
  void SetSigmaStepMethodToEquispaced()
  { this->SetSigmaStepMethod(Self::EquispacedSigmaSteps); }
  void SetSigmaStepMethodToLogarithmic()
  { this->SetSigmaStepMethod(Self::LogarithmicSigmaSteps); }

  /** Get the filter used to compute the Hessian based measure */
  HessianToMeasureFilterType* GetHessianToMeasureFilter()
  {
    return m_HessianToMeasureFilter;
  }

  /** Get the image containing the scales at which each pixel gave the best response*/
  OutputImageType* GetScalesOutput()
  { return  this->GetOutput(1); };

  /** Set the image containing the scales at which each pixel gave the best response*/
  void SetScalesOutput(OutputImageType *scalesImage)
  { this->SetNthOutput(1, scalesImage); };

  /** Get the image containing the Hessian computed at the best response scale*/
  HessianImageType* GetHessianOutput()
  { return static_cast<HessianImageType*>(this->ProcessObject::GetOutput(2)); };

  /** Set the image containing the Hessian computed at the best response scale*/
  void SetHessianOutput(HessianImageType *hessianImage)
  { this->SetNthOutput(2, hessianImage); };

  itkSetMacro(GenerateScalesOutput,bool);
  itkGetMacro(GenerateScalesOutput,bool);
  itkBooleanMacro(GenerateScalesOutput);

  itkSetMacro(GenerateHessianOutput,bool);
  itkGetMacro(GenerateHessianOutput,bool);
  itkBooleanMacro(GenerateHessianOutput);

protected:
  MultiScaleHessianBasedMeasureImageFilter();
  ~MultiScaleHessianBasedMeasureImageFilter() {};
public:
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

private:
  void UpdateMaximumResponse(double sigma);

  double ComputeSigmaValue(int scaleLevel);
  
  void AllocateUpdateBuffer();

  //purposely not implemented
  MultiScaleHessianBasedMeasureImageFilter(const Self&); 
  void operator=(const Self&); //purposely not implemented

  double m_SigmaMin;
  double m_SigmaMax;

  int m_NumberOfSigmaSteps;
  SigmaStepMethodType m_SigmaStepMethod;

  typename HessianToMeasureFilterType::Pointer m_HessianToMeasureFilter;
  typename HessianFilterType::Pointer m_HessianFilter;

  typename UpdateBufferType::Pointer m_UpdateBuffer;

  bool m_GenerateScalesOutput;
  bool m_GenerateHessianOutput;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiScaleHessianBasedMeasureImageFilter.cxx"
#endif
  
#endif
