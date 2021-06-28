CC = g++
CFLAGS = -Wall -O3 --std=c++11

all : esercizio_5.1.exe esercizio_5.2.exe

esercizio_5.1.exe : esercizio_5.1.o random.o Posizione.o 
	$(CC) Posizione.o random.o esercizio_5.1.o -o esercizio_5.1.exe
esercizio_5.1.o : esercizio_5.1.cpp Posizione.h random.h funzioni.h
	$(CC) -c esercizio_5.1.cpp -o esercizio_5.1.o $(CFLAGS)
esercizio_5.2.exe : esercizio_5.2.o random.o Posizione.o 
	$(CC) Posizione.o random.o esercizio_5.2.o -o esercizio_5.2.exe
esercizio_5.2.o : esercizio_5.2.cpp Posizione.h random.h funzioni.h
	$(CC) -c esercizio_5.2.cpp -o esercizio_5.2.o $(CFLAGS)	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Posizione.o : Posizione.cpp Posizione.h
	$(CC) -c Posizione.cpp -o Posizione.o $(CFLAGS)	
integralMC.o : integralMc.cpp integralMC.h FunzioneBase.h
	$(CC) -c integralMC.cpp -o integralMC.o $(CFLAGS)

clean:
	rm *.o esercizio_5.2.exe 
	rm *.dat