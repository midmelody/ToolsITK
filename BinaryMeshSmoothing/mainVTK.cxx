#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataWriter.h>
#include <iostream>
#include "math.h"
#include "string.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile " << std::endl;
      return EXIT_FAILURE;
    }

  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();

  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

  //Pipeline
  reader->Update();
  
  smoother->SetInputConnection(reader->GetOutputPort());
  smoother->SetNumberOfIterations(20);
  smoother->BoundarySmoothingOn();
  smoother->FeatureEdgeSmoothingOn();
  smoother->SetFeatureAngle(120.0);
  smoother->SetPassBand(.001);
  smoother->NonManifoldSmoothingOn();
  smoother->NormalizeCoordinatesOn();
  smoother->Update();

  writer->SetInputConnection( smoother->GetOutputPort() );
  writer->Update();
 
  return EXIT_SUCCESS;
}
