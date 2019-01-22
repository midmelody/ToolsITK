//Function for computing single object kappa-connectivity scene - Dijkstra kFOEMS algorithm
void kFOEMS(VectorType * inputImage, ImageType::SpacingType spacing, VectorType sigma, VectorType mean, unsigned int *seedList, unsigned int seedCount, unsigned int *fuzzyConImage, bool *maskImage, bool *lungImage)
{ 
  /*
  int obsX = 251;
  int obsY = 201;
  int obsZ = 234;
  int obsIndex = obsZ*pRow*pCol+obsY*pRow+obsX;
  */
 
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
  float intAffinity, intP, intQ, intPQ, intSigma, intMean;
  float vesAffinity, vesP, vesQ, vesPQ, vesSigma, vesMean;
  float gryAffinity, gryP, gryQ, gryPQ, grySigma, gryMean;
  float scaAffinity, scaP, scaQ, scaPQ, scaThre;
  //Get value
  normalize = spacing[0]*spacing[0]+spacing[1]*spacing[1]+spacing[2]*spacing[2];
  intSigma = float(sigma.intensity);
  vesSigma = float(sigma.vesselness);
  grySigma = float(sigma.grayReconst);
  intMean = float(mean.intensity);
  vesMean = float(mean.vesselness);
  scaThre = float(mean.scale);
  gryMean = float(mean.grayReconst);

  std::cout<<"Parameters: "<<std::endl;
  std::cout<<"IntSigma: "<<intSigma<<"; intMean: "<<intMean<<"; vesSigma: "<<vesSigma<<"; vesMean: "<<vesMean<<"; grySigma: "<<grySigma<<"; gryMean: "<<gryMean<<"; scaThre: "<<scaThre<<std::endl;

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
	      if((iNeighbor>=0)&&(iNeighbor<pRow)&&(jNeighbor>=0)&&(jNeighbor<pCol)&&(kNeighbor>=0)&&(kNeighbor<pSli)&&(maskImage(iNeighbor,jNeighbor,kNeighbor)))
		{
		  neiIndex = kNeighbor*pRow*pCol+jNeighbor*pRow+iNeighbor;
		  //Adjacency
		  spacX = spacing[0] * std::abs(a);
		  spacY = spacing[1] * std::abs(b);
		  spacZ = spacing[2] * std::abs(c);
		  adjacency = 1 / (1 + 0.05*(spacX*spacX+spacY*spacY+spacZ*spacZ)/normalize);
		  
		  //Affinity - Intensity			  
		  intP = float(inputImage[curIndex].intensity);
		  intQ = float(inputImage[neiIndex].intensity);
		  intAffinity = exp(-(intQ-intP)*(intQ-intP)/(intSigma*intSigma));
		  //bright object, attenuate if value smaller than mean
		  //choose larger one in calculation ( less attenuation )
		  if(intP > intQ)
		    intPQ = intP;
		  else
		    intPQ = intQ;
		  if(intPQ > intMean)
		    intAffinity = intAffinity;
		  else
		    intAffinity = intAffinity * exp(-(intPQ-intMean)*(intPQ-intMean)/(intSigma*intSigma));
		  intAffinity = sqrt(intAffinity);
		  intAffinity = intAffinity*adjacency;

		  //Affinity - Vesselness  
		  vesP = float(inputImage[curIndex].vesselness);
		  vesQ = float(inputImage[neiIndex].vesselness);
		  vesAffinity = exp(-(vesQ-vesP)*(vesQ-vesP)/(vesSigma*vesSigma));
		  if(vesP > vesQ)
		    vesPQ = vesP;
		  else
		    vesPQ = vesQ;
		  if(vesPQ > vesMean)
		    vesAffinity = vesAffinity;
		  else
		    vesAffinity = vesAffinity * exp(-(vesPQ-vesMean)*(vesPQ-vesMean)/(vesSigma*vesSigma));
		  vesAffinity = sqrt(vesAffinity);
		  vesAffinity = vesAffinity * adjacency;
	
		  //Affinity - Grayscale reconstructed image  
		  gryP = float(inputImage[curIndex].grayReconst);
		  gryQ = float(inputImage[neiIndex].grayReconst);
		  gryAffinity = exp(-(gryQ-gryP)*(gryQ-gryP)/(grySigma*grySigma));
		  if(gryP > gryQ)
		    gryPQ = gryP;
		  else
		    gryPQ = gryQ;
		  if(gryPQ > gryMean)
		    gryAffinity = gryAffinity;
		  else
		    gryAffinity = gryAffinity * exp(-(gryPQ-gryMean)*(gryPQ-gryMean)/(grySigma*grySigma));
		  gryAffinity = sqrt(gryAffinity);
		  gryAffinity = gryAffinity * adjacency;
	  
		  //weighted by scale
		  scaP = float(inputImage[curIndex].scale);
		  scaQ = float(inputImage[neiIndex].scale);
		  if(scaP<scaQ)
		    scaPQ = scaP;
		  else
		    scaPQ = scaQ;
		  
		  if(scaPQ>scaThre)
		    affinity = intAffinity;
		  else
		    //affinity = gryAffinity*(1-scaPQ/scaThre) + intAffinity*scaPQ/scaThre;		    
		    affinity = sqrt(gryAffinity*vesAffinity)*(1-scaPQ/scaThre) + intAffinity*scaPQ/scaThre;		
	
		  if(lungImage[curIndex]||lungImage[neiIndex])
	  	    affinity = sqrt(gryAffinity*vesAffinity)*(1-scaPQ/500) + intAffinity*scaPQ/500;		
	
		  //Kappa
		  kappa = round(affinity*maxAffinity);	
		  
		  /*
		  if(curIndex == obsIndex)
		    {
		      std::cout<<intP<<" "<<intQ<<" "<<std::endl;
		      std::cout<<vesP<<" "<<vesQ<<" "<<std::endl;
		      std::cout<<scaP<<" "<<scaQ<<" "<<std::endl;
		      std::cout<<intMean<<" "<<intSigma<<" "<<std::endl;
		      std::cout<<vesMean<<" "<<vesSigma<<" "<<std::endl;
		      std::cout<<intAffinity<<" "<<vesAffinity<<" "<<kappa<<std::endl;
		      std::cout<<scaThre<<std::endl;
		      std::cout<<"------------------------------------------------"<<std::endl;
		    }
		  */
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
