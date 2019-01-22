COMPILE_OPTIONS		= -D_DEBUG

all: driver

common.o : ../common.h ../common.cpp
	g++ -c ../common.cpp -o common.o $(COMPILE_OPTIONS)
potVect.o : potVect.h potVect.cpp ../common.h
	g++ -c potVect.cpp -o potVect.o $(COMPILE_OPTIONS)

driver: driver.cpp common.o potVect.o
	g++	driver.cpp common.o potVect.o	-o driver $(COMPILE_OPTIONS)


clean:
	rm -f *.o driver *~
