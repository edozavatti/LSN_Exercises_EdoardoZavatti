CC = g++
CFLAGS = -Wall -O3 --std=c++11

make: esercizio_2.2.exe esercizio_2.2.bis.exe

esercizio_2.2.exe : esercizio_2.2.o random.o Posizione.o 
	$(CC) Posizione.o random.o esercizio_2.2.o -o esercizio_2.2.exe
esercizio_2.2.bis.exe : esercizio_2.2.bis.o random.o Posizione.o
	$(CC) Posizione.o random.o esercizio_2.2.bis.o -o esercizio_2.2.bis.exe
esercizio_2.2.bis.o : esercizio_2.2.bis.cpp funzioni.h Posizione.h random.h
	$(CC) -c esercizio_2.2.bis.cpp -o esercizio_2.2.bis.o $(CFLAGS)
esercizio_2.2.o : esercizio_2.2.cpp Posizione.h random.h funzioni.h
	$(CC) -c esercizio_2.2.cpp -o esercizio_2.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Posizione.o : Posizione.cpp Posizione.h
	$(CC) -c Posizione.cpp -o Posizione.o $(CFLAGS)	

clean:
	rm *.o esercizio_2.2.exe 
	rm *.dat