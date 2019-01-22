#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#define PI 3.1415926
#include "Filtering3D.h"

int main(int argc, char *argv[])
{
  //Necessary parameters
  int i,j,k,a; //counter
  float maxScale = 0; //max scale value
  float sigmaNoise; //overall noise level
  int iteration; //iteration number of filtering
  //Get command line
  if( argc < 7 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage TScaleFolder outputImage noiseSigma iteration maxDT" << std::endl;
      return EXIT_FAILURE;
    }
  sigmaNoise = atof(argv[4]);
  iteration = atoi(argv[5]);
  float maxDT = atof(argv[6]);

  char fileName[100];

  //Read  shape tensor file and image
  ImageReaderType::Pointer ImageReader = ImageReaderType::New();
  ImageReader->SetFileName(argv[1]);
  ImageReader->Update();
  ImageReaderType::Pointer Size1Reader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Size1.img");
  Size1Reader->SetFileName(fileName);
  Size1Reader->Update();
  ImageReaderType::Pointer Ori1XReader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Ori1X.img");
  Ori1XReader->SetFileName(fileName); 
  Ori1XReader->Update(); 
  ImageReaderType::Pointer Ori1YReader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Ori1Y.img");
  Ori1YReader->SetFileName(fileName); 
  Ori1YReader->Update();
  ImageReaderType::Pointer Ori1ZReader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Ori1Z.img");
  Ori1ZReader->SetFileName(fileName); 
  Ori1ZReader->Update();
  ImageReaderType::Pointer Size2Reader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Size2.img");
  Size2Reader->SetFileName(fileName); 
  Size2Reader->Update();
  ImageReaderType::Pointer Ori2XReader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Ori2X.img");
  Ori2XReader->SetFileName(fileName); 
  Ori2XReader->Update(); 
  ImageReaderType::Pointer Ori2YReader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Ori2Y.img");
  Ori2YReader->SetFileName(fileName); 
  Ori2YReader->Update();
  ImageReaderType::Pointer Ori2ZReader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Ori2Z.img");
  Ori2ZReader->SetFileName(fileName); 
  Ori2ZReader->Update();
  ImageReaderType::Pointer Size3Reader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Size3.img");
  Size3Reader->SetFileName(fileName); 
  Size3Reader->Update();
  ImageReaderType::Pointer Ori3XReader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Ori3X.img");
  Ori3XReader->SetFileName(fileName); 
  Ori3XReader->Update(); 
  ImageReaderType::Pointer Ori3YReader = ImageReaderType::New();
  strcpy(fileName,argv[2]);
  strcat(fileName,"TensorScale/Ori3Y.img");
  Ori3YReader->SetFileName(fileName); 
  Ori3YReader->Update();
  ImageReaderType::Pointer Ori3ZReader = ImageReaderType::New();
  strcpy(fileName,argv[2]); 
  strcat(fileName,"TensorScale/Ori3Z.img");
  Ori3ZReader->SetFileName(fileName); 
  Ori3ZReader->Update();
  //Get image specs
  ImageType::SpacingType spacing = ImageReader->GetOutput()->GetSpacing();
  spaRow = spacing[0];
  spaCol = spacing[1];
  spaSli = spacing[2];

  ImageType::SizeType  size = ImageReader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];

  ImageType::IndexType start = ImageReader->GetOutput()->GetRequestedRegion().GetIndex();
  ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  
  //Read shape tensor file and image
  VECTOR * vector1  = (VECTOR *)malloc(sizeof(VECTOR) * pRow * pCol * pSli); 
  VECTOR * vector2  = (VECTOR *)malloc(sizeof(VECTOR) * pRow * pCol * pSli); 
  VECTOR * vector3  = (VECTOR *)malloc(sizeof(VECTOR) * pRow * pCol * pSli); 
  float * image = (float *)malloc(sizeof(float) * pRow * pCol * pSli); 
  ImageType::IndexType pixelIndex;
  
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  image(i,j,k) = ImageReader->GetOutput()->GetPixel(pixelIndex);	
	  vector1(i,j,k).x = Ori1XReader->GetOutput()->GetPixel(pixelIndex);
	  vector1(i,j,k).y = Ori1YReader->GetOutput()->GetPixel(pixelIndex);
	  vector1(i,j,k).z = Ori1ZReader->GetOutput()->GetPixel(pixelIndex);
	  vector1(i,j,k).length = Size1Reader->GetOutput()->GetPixel(pixelIndex);
	  vector2(i,j,k).x = Ori2XReader->GetOutput()->GetPixel(pixelIndex);
	  vector2(i,j,k).y = Ori2YReader->GetOutput()->GetPixel(pixelIndex);
	  vector2(i,j,k).z = Ori2ZReader->GetOutput()->GetPixel(pixelIndex);
	  vector2(i,j,k).length = Size2Reader->GetOutput()->GetPixel(pixelIndex);
	  vector3(i,j,k).x = Ori3XReader->GetOutput()->GetPixel(pixelIndex);
	  vector3(i,j,k).y = Ori3YReader->GetOutput()->GetPixel(pixelIndex);
	  vector3(i,j,k).z = Ori3ZReader->GetOutput()->GetPixel(pixelIndex);
	  vector3(i,j,k).length = Size3Reader->GetOutput()->GetPixel(pixelIndex);
	  if(maxScale<vector1(i,j,k).length) maxScale=vector1(i,j,k).length;
	  if(maxScale<vector2(i,j,k).length) maxScale=vector2(i,j,k).length;
	  if(maxScale<vector3(i,j,k).length) maxScale=vector3(i,j,k).length;
	}
  if(maxScale<maxDT) maxScale=maxDT;
  
  //----------------------------------------------------------------------------------
  //Compute diffusion magnitudes (ellipsoid length along 13 directions)
  DIFFUSIZE * diffuSize = (DIFFUSIZE *)malloc(sizeof(DIFFUSIZE) * pRow * pCol * pSli);
  VECTOR vec1, vec2, vec3; //ellipsoid data
  VECTOR t,v; //temp variable
  std::cout<<"Compute diffusion size: "<<std::endl;
  //Compute ellipsoid radius along 13 neighborhood directions (other 13 using symmetry to get)
  for(i=0;i<pRow;i++)
    {
      std::cout<<"\r";
      std::cout<<int((i+1)*100/pRow)<<"%"<<std::flush;
      for(j=0;j<pCol;j++)
	for(k=0;k<pSli;k++)
	  {
	    vec1.x = vector3(i,j,k).x;
	    vec1.y = vector3(i,j,k).y;
	    vec1.z = vector3(i,j,k).z;
	    vec1.length = vector3(i,j,k).length;
	    vec2.x = vector2(i,j,k).x;
	    vec2.y = vector2(i,j,k).y;
	    vec2.z = vector2(i,j,k).z;
	    vec2.length = vector2(i,j,k).length;
	    vec3.x = vector1(i,j,k).x;
	    vec3.y = vector1(i,j,k).y;
	    vec3.z = vector1(i,j,k).z;
	    vec3.length = vector1(i,j,k).length;
	    
	    diffuSize(i,j,k).upD  = elliRadius(vec1, vec2, vec3,  0,  0, -1);
	    diffuSize(i,j,k).upN  = elliRadius(vec1, vec2, vec3,  0, -1, -1);
	    diffuSize(i,j,k).upS  = elliRadius(vec1, vec2, vec3,  0,  1, -1);
	    diffuSize(i,j,k).upW  = elliRadius(vec1, vec2, vec3, -1,  0, -1);
	    diffuSize(i,j,k).upE  = elliRadius(vec1, vec2, vec3,  1,  0, -1);
	    diffuSize(i,j,k).upNW = elliRadius(vec1, vec2, vec3, -1, -1, -1);
	    diffuSize(i,j,k).upNE = elliRadius(vec1, vec2, vec3,  1, -1, -1);
	    diffuSize(i,j,k).upSW = elliRadius(vec1, vec2, vec3, -1,  1, -1);
	    diffuSize(i,j,k).upSE = elliRadius(vec1, vec2, vec3,  1,  1, -1);
	    diffuSize(i,j,k).cuN  = elliRadius(vec1, vec2, vec3,  0, -1,  0);
	    diffuSize(i,j,k).cuE  = elliRadius(vec1, vec2, vec3,  1,  0,  0);
	    diffuSize(i,j,k).cuNW = elliRadius(vec1, vec2, vec3, -1, -1,  0);
	    diffuSize(i,j,k).cuNE = elliRadius(vec1, vec2, vec3,  1, -1,  0);
	    /*
	    std::cout<<diffuSize(i,j,k).upD<<std::endl;
	    std::cout<<diffuSize(i,j,k).upN<<std::endl;
	    std::cout<<diffuSize(i,j,k).upS<<std::endl;
	    std::cout<<diffuSize(i,j,k).upW<<std::endl;
	    std::cout<<diffuSize(i,j,k).upE<<std::endl;
	    std::cout<<diffuSize(i,j,k).upNW<<std::endl;
	    std::cout<<diffuSize(i,j,k).upNE<<std::endl;
	    std::cout<<diffuSize(i,j,k).upSW<<std::endl;
	    std::cout<<diffuSize(i,j,k).upSE<<std::endl;
	    std::cout<<diffuSize(i,j,k).cuN<<std::endl;
	    std::cout<<diffuSize(i,j,k).cuE<<std::endl;
	    std::cout<<diffuSize(i,j,k).cuNW<<std::endl;
	    std::cout<<diffuSize(i,j,k).cuNE<<std::endl;
	    */
	  }
    }
  free(vector1);
  free(vector2);
  free(vector3);
  std::cout<<std::endl;
  
  //------------------------------------------------------------------------------------
  /*
  //Smooth the fields
  GaussFilterType::Pointer filterX = GaussFilterType::New();
  GaussFilterType::Pointer filterY = GaussFilterType::New();
  GaussFilterType::Pointer filterZ = GaussFilterType::New();
  filterX->SetDirection( 0 );   // 0 --> X direction
  filterY->SetDirection( 1 );   // 1 --> Y direction
  filterZ->SetDirection( 2 );   // 2 --> Z direction
  filterX->SetOrder( GaussFilterType::ZeroOrder );
  filterY->SetOrder( GaussFilterType::ZeroOrder );
  filterZ->SetOrder( GaussFilterType::ZeroOrder );
  filterX->SetSigma( spaRow*2 );
  filterY->SetSigma( spaCol*2 );
  filterZ->SetSigma( spaSli*2 );
  ImageType::Pointer sizeImage = ImageType::New();
  sizeImage->SetRegions( region );
  sizeImage->Allocate();
  //Link pipeline
  filterX->SetInput( sizeImage );
  filterY->SetInput( filterX->GetOutput() );
  filterZ->SetInput( filterY->GetOutput() );
  //Apply 13 data sets
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upD);
	  //if((i==12)&&(j==132)&&(k==32)) std::cout<<diffuSize(i,j,k).upD<<std::endl;
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upD = filterZ->GetOutput()->GetPixel(pixelIndex);\
	  std::cout<<diffuSize(i,j,k).upD<<std::endl;
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upN);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upN = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upS);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upS = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upE);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upE = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upW);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upW = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upNW);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upNW = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upNE);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upNE = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upSW);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upSW = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).upSE);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).upSE = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).cuN);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).cuN = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).cuE);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).cuE = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).cuNW);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).cuNW = filterZ->GetOutput()->GetPixel(pixelIndex);
	}

  sizeImage->DisconnectPipeline();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  sizeImage->SetPixel(pixelIndex,diffuSize(i,j,k).cuNE);
	}
  filterZ->Update();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  diffuSize(i,j,k).cuNE = filterZ->GetOutput()->GetPixel(pixelIndex);
	}  
  */
  //------------------------------------------------------------------------------------
  //Perform filtering for n iterations
  for (a=0; a<iteration; a++) 
    {
      std::cout<<"Loop: "<<a+1<<std::endl;
      shapeTensorFiltering(diffuSize, image, maxScale, sigmaNoise);
      std::cout<<std::endl;
    }
  //------------------------------------------------------------------------------------
  
  //Write output
  ImageType::Pointer OutImage = ImageType::New();
  OutImage->SetRegions( region );
  OutImage->SetSpacing(ImageReader->GetOutput()->GetSpacing());
  OutImage->Allocate();
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++)
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=k;
	  OutImage->SetPixel(pixelIndex,image(i,j,k));
	}
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( OutImage );
  writer->SetFileName( argv[3] );
  writer->Update();

  free(diffuSize);
  free(image);
}
