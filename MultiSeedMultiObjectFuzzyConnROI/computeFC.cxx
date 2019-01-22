//Function for computing and recording affinity relationship over whole image
float computeAffinityMap(AffinityType * affinityMap, float * inputImage, ImageType::SpacingType spacing, float sigma, bool * ROIImage)
{
  int i, j, k;
  int a, b, c;
  int iNeighbor, jNeighbor, kNeighbor;
  unsigned int curIndex, neiIndex, recIndex;
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
	  curIndex = k*pRow*pCol+j*pRow+i;
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
		  //check range
		  if((iNeighbor>=0)&&(iNeighbor<pRow)&&(jNeighbor>=0)&&(jNeighbor<pCol)&&(kNeighbor>=0)&&(kNeighbor<pSli)&&(ROIImage(i,j,k)))
		    {
		      neiIndex = kNeighbor*pRow*pCol+jNeighbor*pRow+iNeighbor;
		      //Adjacency
		      spacX = spacing[0] * abs(a);
		      spacY = spacing[1] * abs(b);
		      spacZ = spacing[2] * abs(c);
		      adjacency = 1 / (1 + 0.05*(spacX*spacX+spacY*spacY+spacZ*spacZ)/normalize);
		      //Affinity			  
		      intensityP = inputImage[curIndex];
		      intensityQ = inputImage[neiIndex];
		      affinity =  adjacency * exp(-(intensityQ-intensityP)*(intensityQ-intensityP)/(sigma*sigma));
		      //Kappa
		      kappa = round(affinity*maxAffinity);
		      //Record in the map
		      affinityMap[curIndex].affinity[recIndex] = kappa;  
		    }
		  else
		    affinityMap[curIndex].affinity[recIndex] = 0;  
		}	  
	}
}
 

//Function for computing single object kappa-connectivity scene - Dijkstra kFOEMS algorithm
void kFOEMS( AffinityType * affinityMap, SeedListType & seedList, unsigned short *fuzzyConImage, bool * ROIImage)
{ 
  //Compute Fuzzy Connectivity
  std::cout <<"Computing kFOEMS.......";
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
	      if(ROIImage(iNeighbor,jNeighbor,kNeighbor))
		{
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
	    }
      //current strongest FC all processes, move to next level under current max
      while((topIndex >= 0) && (candidateQueue[topIndex]->empty()))  
	topIndex--;     
    }
  std::cout <<"finished"<<std::endl;
}

//Function for computing Iterative Multi-Seed Multi-Object Fuzzy Connectivity - kIRMOFC algorithm
void kIRMOFC( AffinityType * affinityMap, SeedMapType & seedMap, SeedLabelListType & seedLabelList, unsigned short *labelSegImage, bool * ROIImage)
{
  int volumeSize = pRow * pCol * pSli;
  unsigned int curIndex, recIndex, neiIndex;
  int i,j,k;
  int a,b,c;
  int iNeighbor,jNeighbor,kNeighbor;
  std::cout <<"Computing kIRMOFC"<<std::endl;
  std::cout<<"-----------------------------------------"<<std::endl; 
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
      /*
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
      */
      //marker and FC image over candidate seed (K and S)
      bool *markerImage;
      markerImage = (bool *)malloc(sizeof(bool) * pRow * pCol * pSli); 
      memset(markerImage,0,pRow*pCol*pSli*sizeof(bool));
      unsigned short *KSImage;
      KSImage = (unsigned short *)malloc(sizeof(unsigned short) * pRow * pCol * pSli); 
      memset(KSImage,0,pRow*pCol*pSli*sizeof(unsigned short));   
      //compute FC over affinityMap and candidate seed
      std::cout<<"Compute FC_KS"<<std::endl;
      kFOEMS( affinityMap, seedListS, KSImage, ROIImage);

      //duplicate affinityMap for constrained affinityMap
      AffinityType *affinityMapS;
      affinityMapS = (AffinityType *)malloc(sizeof(AffinityType) * pRow * pCol * pSli); 
      for(i=0;i<volumeSize;i++)
	for(j=0;j<27;j++)
	  {
	    affinityMapS[i].affinity[j] = affinityMap[i].affinity[j];
	  }

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
	  kFOEMS( affinityMapS, seedListW, KsWImage, ROIImage);

	  //go through image 
	  for(curIndex=0; curIndex<volumeSize; curIndex++)
	    {
	      //hasn't been marked and KS stronger than KsW, mark
	      if((!markerImage[curIndex])&&(KSImage[curIndex]>=KsWImage[curIndex]))
		{
		  markerImage[curIndex] = true;
		  flag = true;
		  //Get 3D index
		  int i = curIndex % pRow;
		  int j = ((curIndex-i) / pRow) % pCol;
		  int k = (curIndex-i-j*pRow) / (pRow*pCol);
		  //constrain affinityMap by setting affinity value associated with i to 0
		  for(a=-1;a<=1;a++)
		    for(b=-1;b<=1;b++)
		      for(c=-1;c<=1;c++)
			{
			  if((abs(a)+abs(b)+abs(c))>0)//exclude candidate itself
			    {
			      iNeighbor = i+a;
			      jNeighbor = j+b;
			      kNeighbor = k+c;		

			      //disable i-neighbor link
			      recIndex = (c+1)*9 + (b+1)*3 +(a+1);
			      affinityMapS[curIndex].affinity[recIndex] = 0;

			      //disable neighbor-i link
			      //check range
			      if((iNeighbor>=0)&&(iNeighbor<pRow)&&(jNeighbor>=0)&&(jNeighbor<pCol)&&(kNeighbor>=0)&&(kNeighbor<pSli))
				{
				  recIndex = (-c+1)*9 + (-b+1)*3 +(-a+1);
				  neiIndex = kNeighbor*pRow*pCol+jNeighbor*pRow+iNeighbor;	    
				  affinityMapS[neiIndex].affinity[recIndex] = 0;   
				}
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
  std::cout<<"finished"<<std::endl;
  std::cout<<std::endl;
}
