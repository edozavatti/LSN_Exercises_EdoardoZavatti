CC = g++
CFLAGS = -Wall -O3 --std=c++11

all : esercizio_8.2.exe

esercizio_8.2.exe : esercizio_8.2.o random.o Posizione.o integralMC.o
	$(CC) IntegralMC.o Posizione.o random.o esercizio_8.2.o -o esercizio_8.2.exe
esercizio_8.2.o : esercizio_8.2.cpp Posizione.h random.h funzioni.h
	$(CC) -c esercizio_8.2.cpp -o esercizio_8.2.o $(CFLAGS)
sercizio_8.2.o : esercizio_8.2.cpp Posizione.h random.h funzioni.h
	$(CC) -c esercizio_8.2.cpp -o esercizio_8.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Posizione.o : Posizione.cpp Posizione.h
	$(CC) -c Posizione.cpp -o Posizione.o $(CFLAGS)	
integralMC.o : integralMc.cpp integralMC.h FunzioneBase.h
	$(CC) -c integralMC.cpp -o integralMC.o $(CFLAGS)

clean:
	rm *.o  
	rm *.dat

