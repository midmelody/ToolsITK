//Function for computing and recording affinity relationship over whole image
float computeAffinityMap(AffinityType * affinityMap, VectorType * inputImage, ImageType::SpacingType spacing, VectorType sigma, VectorType mean, VectorType darkObject)
{
  int i, j, k; 
  int a, b, c;
  int iNeighbor, jNeighbor, kNeighbor;
  unsigned int curIndex, neiIndex, recIndex;
  unsigned short kappa;
  float spacX, spacY, spacZ;
  float normalize;
  float adjacency, affinity;
  //stats for intensity, vesselness and scale
  float intAffinity, intP, intQ, intSigma, intMean;
  float vesAffinity, vesP, vesQ, vesSigma, vesMean;
  float scaAffinity, scaP, scaQ, scaSigma, scaMean;
  bool intFlag, vesFlag, scaFlag;

  normalize = spacing[0]*spacing[0]+spacing[1]*spacing[1]+spacing[2]*spacing[2];

  intSigma = float(sigma.intensity);
  vesSigma = float(sigma.vesselness);
  scaSigma = float(sigma.scale);

  intMean = float(mean.intensity);
  vesMean = float(mean.vesselness);
  scaMean = float(mean.scale);

  intFlag = darkObject.intensity;
  vesFlag = darkObject.vesselness;
  scaFlag = darkObject.scale;

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
		  if((iNeighbor>=0)&&(iNeighbor<pRow)&&(jNeighbor>=0)&&(jNeighbor<pCol)&&(kNeighbor>=0)&&(kNeighbor<pSli))
		    {
		      neiIndex = kNeighbor*pRow*pCol+jNeighbor*pRow+iNeighbor;
		      //Adjacency
		      spacX = spacing[0] * abs(a);
		      spacY = spacing[1] * abs(b);
		      spacZ = spacing[2] * abs(c);
		      adjacency = 1 / (1 + 0.05*(spacX*spacX+spacY*spacY+spacZ*spacZ)/normalize);
		      //Affinity			  
		      intP = float(inputImage[curIndex].intensity);
		      intQ = float(inputImage[neiIndex].intensity);
		      intAffinity = exp(-(intQ-intP)*(intQ-intP)/(intSigma*intSigma));
		      if(intFlag)
			{
			  if(intQ < intMean)
			    intAffinity = intAffinity;
			  else
			    intAffinity = intAffinity * exp(-(intQ-intMean)*(intQ-intMean)/(intSigma*intSigma));
			}
		      else
			{
			  if(intQ > intMean)
			    intAffinity = intAffinity;
			  else
			    intAffinity = intAffinity * exp(-(intQ-intMean)*(intQ-intMean)/(intSigma*intSigma));
			}
		      intAffinity = sqrt(intAffinity);
		      intAffinity = intAffinity * adjacency;

		      /*
		      vesP = float(inputImage[curIndex].vesselness);
		      vesQ = float(inputImage[neiIndex].vesselness);
		      vesAffinity = exp(-(vesQ-vesP)*(vesQ-vesP)/(vesSigma*vesSigma));
		      if(vesFlag)
			{
			  if(vesQ < vesMean)
			    vesAffinity = vesAffinity;
			  else
			    vesAffinity = vesAffinity * exp(-(vesQ-vesMean)*(vesQ-vesMean)/(vesSigma*vesSigma));
			}
		      else
			{
			  if(vesQ > vesMean)
			    vesAffinity = vesAffinity;
			  else
			    vesAffinity = vesAffinity * exp(-(vesQ-vesMean)*(vesQ-vesMean)/(vesSigma*vesSigma));
			}
		      vesAffinity = sqrt(vesAffinity);
		      vesAffinity = vesAffinity * adjacency;

		      /*
		      scaP = float(inputImage[curIndex].scale);
		      scaQ = float(inputImage[neiIndex].scale);
		      scaAffinity = exp(-(scaQ-scaP)*(scaQ-scaP)/(scaSigma*scaSigma));
		      if(scaFlag)
			{
			  if(scaQ < scaMean)
			    scaAffinity = scaAffinity;
			  else
			    scaAffinity = scaAffinity * exp(-(scaQ-scaMean)*(scaQ-scaMean)/(scaSigma*scaSigma));
			}
		      else
			{
			  if(scaQ > scaMean)
			    scaAffinity = scaAffinity;
			  else
			    scaAffinity = scaAffinity * exp(-(scaQ-scaMean)*(scaQ-scaMean)/(scaSigma*scaSigma));
			}
		      scaAffinity = sqrt(scaAffinity);
		      scaAffinity = scaAffinity * adjacency;

		      if(scaP>scaMean)
			affinity = intAffinity;
		      else
		      */
		      affinity = intAffinity;
		      //affinity = vesAffinity*(1-scaP/scaMean) + intAffinity*scaP/scaMean;    
		      //affinity = (intAffinity + vesAffinity + scaAffinity)/3;
		      //affinity = (intAffinity + vesAffinity)/2;
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
void kFOEMS( AffinityType * affinityMap, SeedListType & seedList, unsigned short *fuzzyConImage)
{ 
  //Compute Fuzzy Connectivity
  std::cout <<"Computing kFOEMS......."<<std::flush;
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
  std::cout <<"finished"<<std::endl;
}
