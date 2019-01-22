//Function for computing and recording affinity relationship over whole image
float computeAffinityMap(AffinityMapType & affinityMap, float * inputImage, ImageType::SpacingType spacing, float sigma)
{
  int i, j, k;
  int a, b, c;
  int iNeighbor, jNeighbor, kNeighbor;
  unsigned int index, indexNeighbor;
  unsigned short kappa;
  float spacX, spacY, spacZ;
  float normalize;
  float adjacency, affinity;
  float intensityP, intensityQ;
  IndexKappa record; 

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
			  adjacency = 1 / (1 + (spacX*spacX+spacY*spacY+spacZ*spacZ)/normalize);
			  //Affinity			  
			  intensityP = inputImage[index];
			  intensityQ = inputImage[indexNeighbor];
			  affinity =  adjacency * exp(-(intensityQ-intensityP)*(intensityQ-intensityP)/(sigma*sigma));
			  //Kappa
			  kappa = round(affinity*maxAffinity);
			  //Record in the map
			  record.index = indexNeighbor;
			  record.kappa = kappa;
			  affinityMap.insert(std::pair< unsigned int, IndexKappa >(index, record)); 			  
			}
		    }
		}	  
	}
}
 

//Function for computing single object kappa-connectivity scene - Dijkstra kFOEMS algorithm
void kFOEMS( AffinityMapType & affinityMap, SeedListType & seedList, unsigned short *fuzzyConImage)
{ 
  //Compute Fuzzy Connectivity
  std::cout <<"computing kFOEMS"<<std::endl;
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
  unsigned int curIndex;
  unsigned short curConnectivity, nexConnectivity, oriConnectivity;
  while((topIndex >= 0) && (!candidateQueue[topIndex]->empty()))    
    {
      //Pick strongest FC index in queue
      curIndex = candidateQueue[topIndex]->front();
      //remove it from queue
      candidateQueue[topIndex]->erase(candidateQueue[topIndex]->begin());
      //current connectivity
      curConnectivity = fuzzyConImage[curIndex];
      //examine the affinity map
      std::pair< AffinityMapType::iterator, AffinityMapType::iterator > affinityRange;
      //locate range
      affinityRange = affinityMap.equal_range(curIndex);
      //iterate connected voxels to the current candidate point
      AffinityMapType::iterator affinityIt;
      for (affinityIt=affinityRange.first; affinityIt!=affinityRange.second; ++affinityIt)
	{
	  int neighborIndex = (*affinityIt).second.index;
	  oriConnectivity = fuzzyConImage[neighborIndex];
	  int kappa = (*affinityIt).second.kappa;
	  //disabled link has kappa value 0
	  if(kappa > 0)
	    {
	      if(kappa < curConnectivity)
		nexConnectivity = kappa;
	      else
		nexConnectivity = curConnectivity;
	      //update connectivity for neighbors
	      if( nexConnectivity > oriConnectivity)
		{
		  fuzzyConImage[neighborIndex] = nexConnectivity;
		  //push to the candidate queue
		  candidateQueue[nexConnectivity]->push_front(neighborIndex);
		}  
	    }
	}
      //current strongest FC all processes, move to next level under current max
      while((topIndex >= 0) && (candidateQueue[topIndex]->empty()))  
	topIndex--;     
    }
}

//Function for computing Iterative Multi-Seed Multi-Object Fuzzy Connectivity - kIRMOFC algorithm
void kIRMOFC( AffinityMapType & affinityMap, SeedMapType & seedMap, SeedLabelListType & seedLabelList, unsigned short *labelSegImage)
{
  int volumeSize = pRow * pCol * pSli;
  std::cout<<"******************************************"<<std::endl;
  std::cout <<"computing kIRMOFC"<<std::endl<<std::endl;
  //iterate through all groups (labels)
  SeedLabelListType::iterator seedLabelIt;
  SeedListType seedListS;
  SeedListType seedListW;
  for( seedLabelIt=seedLabelList.begin(); seedLabelIt!=seedLabelList.end(); ++seedLabelIt)
    {
      int label = *seedLabelIt;
      std::cout<<"Current Label: "<<label<<std::endl;

      //iterate through seedMap, divide to two groups, S and W
      seedListS.clear();
      seedListW.clear();
      SeedMapType::iterator seedMapIt;
      int seedIndex;
      for(seedMapIt=seedMap.begin(); seedMapIt!=seedMap.end(); ++seedMapIt)
	{
	  seedIndex = (*seedMapIt).second;
	  if((*seedMapIt).first==label)
	    seedListS.push_front(seedIndex);
	  else
	    seedListW.push_front(seedIndex);	    
	}

      //Print seed list
      SeedListType::iterator seedListIt;
      std::cout<<"Candidate seeds:"<<std::endl;
      for(seedListIt=seedListS.begin(); seedListIt!=seedListS.end(); ++seedListIt)
	{
	  seedIndex = *seedListIt;
	  int xN = seedIndex % pRow;
	  int yN = ((seedIndex-xN) / pRow) % pCol;
	  int zN = (seedIndex-xN-yN*pRow) / (pRow*pCol);
	  std::cout<<"("<<xN<<" "<<yN<<" "<<zN<<"); ";
	}
      std::cout<<std::endl;
      std::cout<<"Background seeds:"<<std::endl;
      for(seedListIt=seedListW.begin(); seedListIt!=seedListW.end(); ++seedListIt)
	{
	  seedIndex = *seedListIt;
	  int xN = seedIndex % pRow;
	  int yN = ((seedIndex-xN) / pRow) % pCol;
	  int zN = (seedIndex-xN-yN*pRow) / (pRow*pCol);
	  std::cout<<"("<<xN<<" "<<yN<<" "<<zN<<"); ";
	}
      std::cout<<std::endl;

      //marker and FC image over candidate seed (K and S)
      bool *markerImage;
      markerImage = (bool *)malloc(sizeof(bool) * pRow * pCol * pSli); 
      memset(markerImage,0,pRow*pCol*pSli*sizeof(bool));
      unsigned short *KSImage;
      KSImage = (unsigned short *)malloc(sizeof(unsigned short) * pRow * pCol * pSli); 
      memset(KSImage,0,pRow*pCol*pSli*sizeof(unsigned short));   
      //compute FC over affinityMap and candidate seed
      std::cout<<"Compute FC_KS"<<std::endl;
      kFOEMS( affinityMap, seedListS, KSImage);

      //duplicate affinityMap for constrained affinityMap
      AffinityMapType affinityMapS = affinityMap;
      //FC Image over constrained affinityMap and background seed
      unsigned short *KsWImage;
      KsWImage = (unsigned short *)malloc(sizeof(unsigned short) * pRow * pCol * pSli);
      memset(KsWImage,0,pRow*pCol*pSli*sizeof(unsigned short)); 
      
      //Iteration
      bool flag = true;
      int count = 0;
      while(flag)
	{
	  flag = false;
	  //Compute FC over Ks and W
	  std::cout<<"Compute FC_KsW time "<<count+1<<std::endl;	  
	  kFOEMS( affinityMapS, seedListW, KsWImage);

	  //go through image 
	  for(int i=0; i<volumeSize; i++)
	    {
	      //hasn't been marked and KS stronger than KsW, mark
	      if((!markerImage[i])&&(KSImage[i]>=KsWImage[i]))
		{
		  markerImage[i] = true;
		  flag = true;
		  //constrain affinityMap by setting affinity value associated with i to 0
		  std::pair< AffinityMapType::iterator, AffinityMapType::iterator > affinityRange;
		  //locate range containing i
		  affinityRange = affinityMapS.equal_range(i);
		  //set kappa to 0
		  AffinityMapType::iterator affinityIt;
		  for (affinityIt=affinityRange.first; affinityIt!=affinityRange.second; ++affinityIt)
		    {
		      int neighborIndex = (*affinityIt).second.index;
		      //set kappa to 0
		      (*affinityIt).second.kappa = 0;
		      //erase the other side
		      std::pair< AffinityMapType::iterator, AffinityMapType::iterator > affinityRangeConuter;
		      //locate range containing neighborIndex
		      affinityRangeConuter = affinityMapS.equal_range(neighborIndex);
		      //look for link neighborIndex-i
		      AffinityMapType::iterator affinityItConuter;
		      for (affinityItConuter=affinityRangeConuter.first; affinityItConuter!=affinityRangeConuter.second; ++affinityItConuter)
			{
			  int counterIndex =  (*affinityItConuter).second.index;
			  if(counterIndex == i)
			    //set kappa to 0
			    (*affinityItConuter).second.kappa = 0;			    
			}
		    }
		}
	    } 
	  count++; 
	}
      std::cout<<count<<" iterations"<<std::endl;
      std::cout<<"-----------------------------------------"<<std::endl;

      for(int i=0; i<volumeSize; i++)
	{
	  if(markerImage[i])
	    labelSegImage[i] = label;
	}
    }
}
