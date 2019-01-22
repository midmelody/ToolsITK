COMPILE_OPTIONS		= -D_DEBUG

TNT_HOME = ../TNT
JAMA_HOME = ../JAMA

all: driver

common.o:	../common.h ../common.cpp
	g++ -c ../common.cpp -o common.o $(COMPILE_OPTIONS)

CritPts.o:	CritPts.cpp CritPts.h
	g++ -c CritPts.cpp -o CritPts.o $(COMPILE_OPTIONS) -I$(TNT_HOME) -I$(JAMA_HOME)

driver:		driver.cpp CritPts.o common.o
	g++ driver.cpp CritPts.o common.o -o driver $(COMPILE_OPTIONS) -I$(TNT_HOME) -I$(JAMA_HOME)

clean:	
	rm -f driver *.o *~
