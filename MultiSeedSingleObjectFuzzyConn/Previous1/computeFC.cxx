//Function for computing and recording affinity relationship over whole image
float computeAffinityMap(AffinityType * affinityMap, float * inputImage, ImageType::SpacingType spacing, float sigma)
{
  int i, j, k;
  int a, b, c;
  int iNeighbor, jNeighbor, kNeighbor;
  unsigned int index, indexNeighbor, indexRecord;
  unsigned short kappa;
  float spacX, spacY, spacZ;
  float normalize;
  float adjacency, affinity;
  float intensityP, intensityQ;
 
  normalize = spacing[0]*spacing[0]+spacing[1]*spacing[1]+spacing[2]*spacing[2];

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
		  iNeighbor = i+a;
		  jNeighbor = j+b;
		  kNeighbor = k+c;
		  //record location in 26 neighborhood
		  indexRecord = (c+1)*9 + (b+1)*3 +(a+1);
		  //check range
		  if((iNeighbor>=0)&&(iNeighbor<pRow)&&(jNeighbor>=0)&&(jNeighbor<pCol)&&(kNeighbor>=0)&&(kNeighbor<pSli))
		    {
		      indexNeighbor = kNeighbor*pRow*pCol+jNeighbor*pRow+iNeighbor;
		      //Adjacency
		      spacX = spacing[0] * abs(a);
		      spacY = spacing[1] * abs(b);
		      spacZ = spacing[2] * abs(c);
		      adjacency = 1 / (1 + 0.05*(spacX*spacX+spacY*spacY+spacZ*spacZ)/normalize);
		      //Affinity			  
		      intensityP = inputImage[index];
		      intensityQ = inputImage[indexNeighbor];
		      affinity =  adjacency * exp(-(intensityQ-intensityP)*(intensityQ-intensityP)/(sigma*sigma));
		      //Kappa
		      kappa = round(affinity*maxAffinity);
		      //Record in the map
		      affinityMap[index].affinity[indexRecord] = kappa;  
		    }
		  else
		    affinityMap[index].affinity[indexRecord] = 0;  
		}	  
	}
}
 

//Function for computing single object kappa-connectivity scene - Dijkstra kFOEMS algorithm
void kFOEMS( AffinityType * affinityMap, SeedListType & seedList, unsigned short *fuzzyConImage)
{ 
  //Compute Fuzzy Connectivity
  std::cout <<"Computing kFOEMS"<<std::endl;
  int volumeSize = pRow * pCol * pSli;
  typedef std::list <unsigned int> ListIndexType;
  std::vector<ListIndexType*> candidateQueue;
  candidateQueue.reserve(maxAffinity+1);
  for(int i=0;i<=maxAffinity;i++)
    candidateQueue.push_back(new ListIndexType);
  int topIndex = maxAffinity;

  // initialization of the seed points
  SeedListType::iterator seedIt;
  for(seedIt=seedList.begin(); seedIt!=seedList.end(); ++seedIt)
    {
      fuzzyConImage[*seedIt] = maxAffinity;
      candidateQueue[maxAffinity]->push_front(*seedIt);
    }

  //go over list
  unsigned int curIndex, recIndex, neiIndex;
  unsigned short curConnectivity, nexConnectivity, oriConnectivity;
  int a,b,c;
  int i,j,k;
  int iNeighbor,jNeighbor,kNeighbor;
  while((topIndex >= 0) && (!candidateQueue[topIndex]->empty()))    
    {
      //Pick strongest FC index in queue
      curIndex = candidateQueue[topIndex]->front();
      //Get 3D index
      int i = curIndex % pRow;
      int j = ((curIndex-i) / pRow) % pCol;
      int k = (curIndex-i-j*pRow) / (pRow*pCol);
      //remove it from queue
      candidateQueue[topIndex]->erase(candidateQueue[topIndex]->begin());
      //current connectivity
      curConnectivity = fuzzyConImage[curIndex];
      //examine the affinity map
      //26 neighborhood
      for(a=-1;a<=1;a++)
	for(b=-1;b<=1;b++)
	  for(c=-1;c<=1;c++)	    
	    {
	      iNeighbor = i+a;
	      jNeighbor = j+b;
	      kNeighbor = k+c;
	      //record location in 26 neighborhood
	      recIndex = (c+1)*9 + (b+1)*3 +(a+1);
	      int kappa = affinityMap[curIndex].affinity[recIndex];
	      //check neighbor with positive affinity
	      if(kappa > 0)
		{
		  neiIndex = kNeighbor*pRow*pCol+jNeighbor*pRow+iNeighbor;
		  oriConnectivity = fuzzyConImage[neiIndex];
		  if(kappa < curConnectivity)
		    nexConnectivity = kappa;
		  else
		    nexConnectivity = curConnectivity;
		  //update connectivity for neighbors
		  if( nexConnectivity > oriConnectivity)
		    {
		      fuzzyConImage[neiIndex] = nexConnectivity;
		      //push to the candidate queue
		      candidateQueue[nexConnectivity]->push_front(neiIndex);
		    }  		  
		}
	    }
      //current strongest FC all processes, move to next level under current max
      while((topIndex >= 0) && (candidateQueue[topIndex]->empty()))  
	topIndex--;     
    }
}
