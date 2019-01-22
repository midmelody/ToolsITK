COMPILE_OPTIONS		= -D_DEBUG

all: driver
common.o:	../common.h ../common.cpp
	g++ -c ../common.cpp -o common.o $(COMPILE_OPTIONS)

StreamLn.o:	StreamLn.cpp StreamLn.h
	g++ -c StreamLn.cpp -o StreamLn.o $(COMPILE_OPTIONS)

driver:		driver.cpp StreamLn.o common.o
	g++ driver.cpp StreamLn.o common.o -o driver $(COMPILE_OPTIONS)

clean:	
	rm -f driver *.o *~
