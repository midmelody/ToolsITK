COMPILE_OPTIONS		= -D_DEBUG

all: driver

common.o : ../common.h ../common.cpp
	g++ -c ../common.cpp -o common.o $(COMPILE_OPTIONS)

expandVol.o : expandVol.h expandVol.cpp
	g++ -c expandVol.cpp -o expandVol.o $(COMPILE_OPTIONS)

driver: driver.cpp common.o expandVol.o
	g++	driver.cpp common.o expandVol.o	-o driver $(COMPILE_OPTIONS)


clean:
	rm -f *.o driver *~
