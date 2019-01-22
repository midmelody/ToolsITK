//Function for computing single object kappa-connectivity scene - Dijkstra kFOEMS algorithm
void kFOEMS(VectorType * inputImage, ImageType::SpacingType spacing, VectorType sigmaHigh, VectorType sigmaLow, VectorType mean, unsigned int * seedList, unsigned int seedCount, unsigned int *fuzzyConImage)
{ 
  //Compute Fuzzy Connectivity
  int volumeSize = pRow * pCol * pSli;
  typedef std::list <unsigned int> ListIndexType;
  std::vector<ListIndexType*> candidateQueue;
  candidateQueue.reserve(maxAffinity+1);
  for(int i=0;i<=maxAffinity;i++)
    candidateQueue.push_back(new ListIndexType);
  int topIndex = maxAffinity;

  // initialization of the seed points 
  int seedIt;
  unsigned int seedIndex;
  for(seedIt=0; seedIt<seedCount; seedIt++)
    {
      seedIndex = seedList[seedIt];
      fuzzyConImage[seedIndex] = maxAffinity;
      candidateQueue[maxAffinity]->push_front(seedIndex);
    }

  //Index of current candidate and neighbor voxel
  unsigned int curIndex, neiIndex;
  //connectivities
  unsigned short curConnectivity, nexConnectivity, oriConnectivity;
  //Variables
  int a,b,c;
  int i,j,k;
  int iNeighbor,jNeighbor,kNeighbor;
  //affinity computation 
  int kappa;
  float spacX, spacY, spacZ;
  float normalize;
  float adjacency, affinity;
  //stats for intensity, vesselness and scale
  float CTAffinity, CTP, CTQ, CTPQ, CTSigmaHigh, CTSigmaLow, CTMean;
  float PETAffinity, PETP, PETQ, PETPQ, PETSigma, PETMean;
  //Get value
  normalize = spacing[0]*spacing[0]+spacing[1]*spacing[1]+spacing[2]*spacing[2];
  CTSigmaHigh = float(sigmaHigh.CT);
  CTSigmaLow = float(sigmaLow.CT);
  PETSigma = float(sigmaHigh.PET);
  CTMean = float(mean.CT);
  PETMean = float(mean.PET);

  std::cout<<std::endl<<"Parameters: "<<std::endl;
  std::cout<<"CTSigmaHigh "<<CTSigmaHigh<<"CTSigmaLow "<<CTSigmaLow<<" CTMean "<<CTMean<<" PETSigma "<<PETSigma<<" PETMean "<<PETMean<<std::endl;

  //go over list
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
      //examine the affinity for 26 neighborhood
      for(a=-1;a<=1;a++)
	for(b=-1;b<=1;b++)
	  for(c=-1;c<=1;c++)
	    {
	      iNeighbor = i+a;
	      jNeighbor = j+b;
	      kNeighbor = k+c;
	      //check range and mask
	      if((iNeighbor>=0)&&(iNeighbor<pRow)&&(jNeighbor>=0)&&(jNeighbor<pCol)&&(kNeighbor>=0)&&(kNeighbor<pSli))
		{
		  neiIndex = kNeighbor*pRow*pCol+jNeighbor*pRow+iNeighbor;
		  //Adjacency
		  spacX = spacing[0] * abs(a);
		  spacY = spacing[1] * abs(b);
		  spacZ = spacing[2] * abs(c);
		  adjacency = 1 / (1 + 0.05*(spacX*spacX+spacY*spacY+spacZ*spacZ)/normalize);
		  //Affinity - CT - attenuate if away from expected value
		  CTP = float(inputImage[curIndex].CT);
		  CTQ = float(inputImage[neiIndex].CT);
		  CTPQ = (CTP+CTQ)/2;
		  if(CTPQ<CTMean)
		    CTAffinity = exp(-(CTQ-CTP)*(CTQ-CTP)/(CTSigmaLow*CTSigmaLow));
		  else
		    CTAffinity = exp(-(CTQ-CTP)*(CTQ-CTP)/(CTSigmaHigh*CTSigmaHigh));
		  //choose smaller one in calculation ( greater attenuation )
		  CTPQ = (CTP+CTQ)/2;
		  if(CTPQ<CTMean)
		    CTAffinity = CTAffinity * exp(-(CTPQ-CTMean)*(CTPQ-CTMean)/(CTSigmaLow*CTSigmaLow));
		  else
		    CTAffinity = CTAffinity * exp(-(CTPQ-CTMean)*(CTPQ-CTMean)/(CTSigmaHigh*CTSigmaHigh));
		  CTAffinity = sqrt(CTAffinity);
		  CTAffinity = CTAffinity * adjacency;

		  //Affinity - PET - attenuate if less than expected value
		  PETP = float(inputImage[curIndex].PET);
		  PETQ = float(inputImage[neiIndex].PET);
		  PETAffinity = exp(-(PETQ-PETP)*(PETQ-PETP)/(PETSigma*PETSigma));
		  if(PETP < PETQ)
		    PETPQ = PETP;
		  else
		    PETPQ = PETQ;
		  if(PETPQ < PETMean)
		    PETAffinity = PETAffinity * exp(-(PETPQ-PETMean)*(PETPQ-PETMean)/(PETSigma*PETSigma));
		  PETAffinity = sqrt(PETAffinity);
		  PETAffinity = PETAffinity * adjacency;
		  
		  //weighted
		  float weiCT = 0.7;
		  float weiPET = 0.3;
		  affinity = PETAffinity*weiPET + CTAffinity*weiCT;
		  float setZero = maxAffinity*sqrt(PETAffinity*CTAffinity);
		  if(setZero<100)
		    affinity = 0;
		  
		  //Kappa
		  kappa = round(affinity*maxAffinity);	
		}
	      else
		kappa = 0;
		
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
