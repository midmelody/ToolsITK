COMPILE_OPTIONS		= -D_DEBUG
## COMPILE_OPTIONS		= 

PotFieldDIR	 =	PotField
HCBDIR		 = 	HCB
StreamLNDIR	 =	StreamLn
CritPtsDIR	 =	CritPts
ExpandVolDIR	 =	ExpandVol
MakeSolidVol_DIR = 	MakeSolidVol
HighDiverg_DIR	 = 	HighDiverg
pfSkelDIR	 =	.

TNT_HOME = TNT
JAMA_HOME = JAMA

all: pfSkel

common.o : common.h common.cpp
	g++ -c common.cpp -o common.o $(COMPILE_OPTIONS)

potVect.o : $(PotFieldDIR)/potVect.h $(PotFieldDIR)/potVect.cpp common.h
	g++ -c $(PotFieldDIR)/potVect.cpp -o potVect.o $(COMPILE_OPTIONS)


StreamLn.o:	$(StreamLNDIR)/StreamLn.cpp $(StreamLNDIR)/StreamLn.h common.h
	g++ -c $(StreamLNDIR)/StreamLn.cpp -o StreamLn.o $(COMPILE_OPTIONS)


CritPts.o:	$(CritPtsDIR)/CritPts.cpp $(CritPtsDIR)/CritPts.h common.h
	g++ -c $(CritPtsDIR)/CritPts.cpp -o CritPts.o $(COMPILE_OPTIONS) -I$(TNT_HOME) -I$(JAMA_HOME)

expandVol.o:	$(ExpandVolDIR)/expandVol.cpp  $(ExpandVolDIR)/expandVol.h common.h
	g++ -c $(ExpandVolDIR)/expandVol.cpp -o expandVol.o $(COMPILE_OPTIONS)

makeSolidVol.o : $(MakeSolidVol_DIR)/makeSolidVol.cpp $(MakeSolidVol_DIR)/makeSolidVol.h
	g++ -c $(MakeSolidVol_DIR)/makeSolidVol.cpp -o makeSolidVol.o $(COMPILE_OPTIONS)

HighDiverg.o : $(HighDiverg_DIR)/HighDiverg.cpp $(HighDiverg_DIR)/HighDiverg.h common.h
	g++ -c $(HighDiverg_DIR)/HighDiverg.cpp -o HighDiverg.o $(COMPILE_OPTIONS)

pfSkel.o : $(pfSkelDIR)/pfSkel.cpp common.h
	g++ -c $(pfSkelDIR)/pfSkel.cpp -o pfSkel.o $(COMPILE_OPTIONS) -I$(TNT_HOME) -I$(JAMA_HOME)

pfSkel: driver.cpp pfSkel.o common.o potVect.o StreamLn.o CritPts.o expandVol.o makeSolidVol.o HighDiverg.o
	g++ driver.cpp pfSkel.o common.o expandVol.o potVect.o CritPts.o StreamLn.o makeSolidVol.o HighDiverg.o $(COMPILE_OPTIONS) -I$(TNT_HOME) -I$(JAMA_HOME) -o pfSkel -lm

clean:
	rm -f *.o pfSkel *~

################################
## Not needed for now because the module is not working:
##
##
## BoundPoints.o : $(HCBDIR)/BoundPoints.h $(HCBDIR)/BoundPoints.cpp $(HCBDIR)/def.h
##	g++ -c $(HCBDIR)/BoundPoints.cpp -o BoundPoints.o $(COMPILE_OPTIONS)
##
## hcBoundHaimeiri3.o : $(HCBDIR)/hcBoundHaimeiri3.cpp $(HCBDIR)/hcBoundHaimeiri3.h $(HCBDIR)/def.h
##	g++ -c $(HCBDIR)/hcBoundHaimeiri3.cpp -o hcBoundHaimeiri3.o $(COMPILE_OPTIONS) -I$(TNT_HOME) -I$(JAMA_HOME)
##
## hcBound.o : $(HCBDIR)/hcBound.h $(HCBDIR)/hcBound.cpp $(HCBDIR)/def.h
##	g++ -c $(HCBDIR)/hcBound.cpp -o hcBound.o $(COMPILE_OPTIONS)
##
###########################
