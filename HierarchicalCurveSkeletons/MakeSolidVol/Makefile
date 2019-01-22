all: driver

makeSolidVol.o : makeSolidVol.cpp makeSolidVol.h
	g++ -c makeSolidVol.cpp -o makeSolidVol.o

driver: makeSolidVol.o driver.cpp
	g++ driver.cpp makeSolidVol.o -o driver

clean:
	rm -f driver *.o *~

