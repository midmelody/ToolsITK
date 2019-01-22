//Function for computing and recording affinity relationship over whole image
float computeAffinityMap(AffinityMapType & affinityMap, float * inputImage, ImageType::SpacingType spacing, float sigma)
{
  int i, j, k;
  int a, b, c;
  int iNeighbor, jNeighbor, kNeighbor;
  unsigned int index, indexNeighbor;
  unsigned short kappa;
  float spacX, spacY, spacZ;
  float adjacency, affinity;
  float intensityP, intensityQ;
  IndexKappa record; 

  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)
      for(k=0;k<pSli;k++) 
	{
	  //index of candidate point
	  index = k*pRow*pCol+j*pRow+i;
	  //26 neighborhood
	  for(a=-1;a<=1;a++)
	    for(b=-1;b<=1;b++)
	      for(c=-1;c<=1;c++)	    
		{
		  //exclude candidate point itself
		  if((abs(a)+abs(b)+abs(c))!=0)
		    {
		      iNeighbor = i+a;
		      jNeighbor = j+b;
		      kNeighbor = k+c;
		      //check range
		      if((iNeighbor>=0)&&(iNeighbor<pRow)&&(jNeighbor>=0)&&(jNeighbor<pCol)&&(kNeighbor>=0)&&(kNeighbor<pSli))
			{

			  indexNeighbor = kNeighbor*pRow*pCol+jNeighbor*pRow+iNeighbor;
			  //Adjacency
			  spacX = spacing[0] * abs(a);
			  spacY = spacing[1] * abs(b);
			  spacZ = spacing[2] * abs(c);
			  adjacency = 1 / (1 + spacX*spacX+spacY*spacY+spacZ*spacZ);
			  //Affinity			  
			  intensityP = inputImage[index];
			  intensityQ = inputImage[indexNeighbor];
			  affinity =  adjacency * exp(-(intensityQ-intensityP)*(intensityQ-intensityP)/(2*sigma*sigma));
			  //Kappa
			  kappa = round(affinity*10000);
			  //Record in the map
			  record.index = indexNeighbor;
			  record.kappa = kappa;
			  affinityMap.insert(std::pair< unsigned int, IndexKappa >(index, record)); 			  
			}
		    }
		}	  
	}
}

/*

//Function for computing single seed set kappa-connectivity scene - Dijkstra kappa-FOEMS algorithm
void computeFC( float *inputImage, float *fuzzyConImage, Voxel *seed, int seedCount, float meanObj, float sigmaObj, float confidentBackground)
{ 
  //Compute Fuzzy Connectivity
  std::cout <<"computing FC"<<std::endl;
  Voxel cur,index,index1;
  float curConnectivity, nextConnectivity;
  int x, y, z,topIndex;
  int slice, row, column;
  char *fuzzyConFlag;
  int volumeSize = pRow * pCol * pSli;

  typedef std::list <Voxel> ListPointType;
  std::vector<ListPointType*> vwp;
  
  Voxel nbor[26] = { { -1, -1, -1 }, { 0, -1, -1 }, { 1, -1, -1},  
		     { -1, 0, -1 },  { 0, 0, -1 },  { 1, 0, -1 }, 
		     { -1, 1, -1 },  { 0, 1, -1 },  { 1, 1, -1 },
		     { -1, -1, 0 },  { 0, -1, 0 },  { 1, -1, 0},  
		     { -1, 0, 0 },                  { 1, 0, 0 }, 
		     { -1, 1, 0 },   { 0, 1, 0 },   { 1, 1, 0 },
		     { -1, -1,1 },   { 0, -1, 1},   { 1, -1, 1},  
		     { -1, 0, 1 },   { 0, 0, 1 },   { 1, 0, 1 }, 
		     { -1, 1, 1 },   { 0, 1, 1 },   { 1, 1, 1 }};
  
		     //----------------------Allocate distance transform array 
		     fuzzyConFlag = new char[volumeSize]; 
  if ( NULL == fuzzyConFlag )
    {
      printf("\r\nFailed to allocate the memory of flagFC!\n");
      return;
    }
  memset(fuzzyConFlag,0,volumeSize*sizeof(char));
  
  for(int i= 0;i<volumeSize;i++)  
    {
      if(inputImage[i] < confidentBackground)  //background, not compute FC for this part
	{
	  fuzzyConImage[i] = 0;
	  fuzzyConFlag[i] = 1;
	}
      else  
	{
	  fuzzyConImage[i] = 0;
	  fuzzyConFlag[i] = 0;
	}
    }
  
  //maxConnectivity 
  int maxConnectivity = 1000;
  assert(maxConnectivity <= maxHeapSize);   
  vwp.reserve(maxConnectivity+1);
  for(int i=0;i<=maxConnectivity;i++)
    vwp.push_back(new ListPointType);
  topIndex = maxConnectivity;
  
  // initialization of the seed points 
  for(int i= 0;i<seedCount;i++)   
    {
      if(inputImage(seed[i].x,seed[i].y,seed[i].z) < confidentBackground)
	{
	  std::cerr<<"seed point is not selected correctly"<<std::endl;
	  assert(0);
	}
      fuzzyConImage(seed[i].x,seed[i].y,seed[i].z) = maxConnectivity;
      fuzzyConFlag(seed[i].x,seed[i].y,seed[i].z) = 0;
      //put all seed points to the list as the start points for vessel
      cur.x = seed[i].x;
      cur.y = seed[i].y;
      cur.z = seed[i].z;
      vwp[maxConnectivity]->push_front(cur); 
    }
  
  float affinity;
  int listIndex;
  while((topIndex >= 0) && (!vwp[topIndex]->empty()))    
    {
      // return the first element in list 
      cur = vwp[topIndex]->front();
      // erase the first element in list 
      vwp[topIndex]->erase(vwp[topIndex]->begin());
      //Current connectivity 
      curConnectivity = fuzzyConImage(cur.x, cur.y, cur.z);
      //Haven't been visited
      if(fuzzyConFlag(cur.x, cur.y, cur.z) == 0) 
	{
	  //Mark as visited
	  //flagFC[cur.z*sliceSize+cur.y*pColumn+cur.x] = 1;
	  //Observe the 26-neighborhood
	  int start_element = 0;
	  int end_element = 26;
	  for(int ei = start_element;ei<end_element;ei++)    
	    {
	      x = cur.x + nbor[ei].x;
	      y = cur.y + nbor[ei].y;
	      z = cur.z + nbor[ei].z;
	      if (x >= 0 && x < pRow && y >= 0 && y < pCol && z >= 0 && z < pSli && fuzzyConFlag(x,y,z) == 0)   
		{
		  float intensityP = inputImage(cur.x, cur.y, cur.z);
		  float intensityQ = inputImage(x, y, z);
		  affinity = computeAffinity(intensityP,intensityQ,meanObj,sigmaObj);
		  listIndex = round(affinity*1000);
		  //Test affinity-based connectivity
		  if (affinity < curConnectivity) 
		    nextConnectivity = listIndex;
		  else 
		    nextConnectivity = curConnectivity;

		  //Update connectivity value
		  int oldValue = round(fuzzyConImage(x, y, z));
		  if( nextConnectivity > oldValue)  
		    {
		      index.x = x;	  
		      index.y = y;
		      index.z = z;
 		      vwp[nextConnectivity]->push_front(index);
		      fuzzyConImage(index.x, index.y, index.z) = nextConnectivity;
		    }
		}
	    }
	  //std::cout<<std::endl;
	}
      while((topIndex >= 0) && (vwp[topIndex]->empty()))  
	topIndex--;
    }

}
*/
//Function for computing Iterative Multi-Seed Multi-Object Fuzzy Connectivity
