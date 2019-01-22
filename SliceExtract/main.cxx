#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
 
#define PI 3.1415926
#include "Compare.h"

int main(int argc, char *argv[])
{
  int i,j; //counters

  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImage sliceNo fileName orientation" << std::endl;
      return EXIT_FAILURE;
    }
 
  //Read image
  ImageReaderType::Pointer Reader = ImageReaderType::New();
  Reader->SetFileName( argv[1] );
  Reader->Update();
  
  //Get image specs
  ImageType::SizeType size = Reader->GetOutput()->GetRequestedRegion().GetSize();
  pRow = size[0];
  pCol = size[1];
  pSli = size[2];

  //Output
  FILE *fileW;
  fileW = fopen(argv[3],"w+");  
  float value;
  ImageType::IndexType pixelIndex;
  int row,col,sli;
  int ori = atoi(argv[4]);
  if(ori==1)
    {
      row = atoi(argv[2])-1;
      fprintf(fileW,"%d\n%d\n", pCol, pSli);
      for(i=0;i<pCol;i++)
	for(j=0;j<pSli;j++)
	  {
	    pixelIndex[0]=row;
	    pixelIndex[1]=i;
	    pixelIndex[2]=j;
	    value = Reader->GetOutput()->GetPixel(pixelIndex);
	    if(value!=value)
	      fprintf(fileW,"%f\n",0);
	    else
	      fprintf(fileW,"%f\n",value);
	  }     
    }
  else if(ori==2)
    {
      col = atoi(argv[2])-1;
      fprintf(fileW,"%d\n%d\n", pRow, pSli);
      for(i=0;i<pRow;i++)
	for(j=0;j<pSli;j++)
	  {
	    pixelIndex[0]=i;
	    pixelIndex[1]=col;
	    pixelIndex[2]=j;
	    value = Reader->GetOutput()->GetPixel(pixelIndex);
	    if(value!=value)
	      fprintf(fileW,"%f\n",0);
	    else
	      fprintf(fileW,"%f\n",value);
	  }     
    }
  else if(ori==3)
    {
      sli = atoi(argv[2])-1;
      fprintf(fileW,"%d\n%d\n", pRow, pCol);
      for(i=0;i<pRow;i++)
	for(j=0;j<pCol;j++)
	  {
	    pixelIndex[0]=i;
	    pixelIndex[1]=j;
	    pixelIndex[2]=sli;
	    value = Reader->GetOutput()->GetPixel(pixelIndex);
	    if(value!=value)
	      fprintf(fileW,"%f\n",0);
	    else
	      fprintf(fileW,"%f\n",value);
	  }     
    }
  else
    {
      std::cerr << "orientation have to be 1,2,3"<< std::endl;
      return EXIT_FAILURE;  
    }
  return EXIT_SUCCESS;
}


