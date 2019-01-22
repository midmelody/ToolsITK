UTILS_CODE	= ../../../../Utils/code
UTILS_INCLUDE	= ../../../../Utils/include

TNT_HOME = /home/cornea/Software/TNT
JAMA_HOME = /home/cornea/Software/JAMA

COMPILE_OPTIONS		= -D_DEBUG

common.o : ../common.h ../common.cpp
	g++ -c ../common.cpp -o common.o $(COMPILE_OPTIONS)

BoundPoints.o : common.o BoundPoints.h BoundPoints.cpp def.h
	g++ -c BoundPoints.cpp -o BoundPoints.o $(COMPILE_OPTIONS)

hcBoundHaimeiri3.o : hcBoundHaimeiri3.cpp hcBoundHaimeiri3.h def.h
	g++ -c hcBoundHaimeiri3.cpp -o hcBoundHaimeiri3.o $(COMPILE_OPTIONS) -I$(TNT_HOME) -I$(JAMA_HOME)

hcBound.o : hcBound.h hcBound.cpp def.h
	g++ -c hcBound.cpp -o hcBound.o $(COMPILE_OPTIONS)
	
	

driver : driver.cpp hcBound.o hcBoundHaimeiri3.o BoundPoints.o common.o 
	g++ driver.cpp hcBound.o hcBoundHaimeiri3.o BoundPoints.o common.o -o driver $(COMPILE_OPTIONS) -lm



clean:
	rm -f *.o driver
