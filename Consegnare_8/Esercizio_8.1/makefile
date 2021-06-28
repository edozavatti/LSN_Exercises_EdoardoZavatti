CC = g++
CFLAGS = -Wall -O3 --std=c++11

all : esercizio_8.1.exe

esercizio_8.1.exe : esercizio_8.1.o random.o Posizione.o integralMC.o
	$(CC) IntegralMC.o Posizione.o random.o esercizio_8.1.o -o esercizio_8.1.exe
esercizio_8.1.o : esercizio_8.1.cpp Posizione.h random.h funzioni.h
	$(CC) -c esercizio_8.1.cpp -o esercizio_8.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Posizione.o : Posizione.cpp Posizione.h
	$(CC) -c Posizione.cpp -o Posizione.o $(CFLAGS)	
integralMC.o : integralMc.cpp integralMC.h FunzioneBase.h
	$(CC) -c integralMC.cpp -o integralMC.o $(CFLAGS)

clean:
	rm *.o  
	rm *.dat

