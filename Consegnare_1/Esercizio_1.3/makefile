CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio_1.3.exe : esercizio_1.3.o random.o
	$(CC) random.o esercizio_1.3.o -o esercizio_1.3.exe
esercizio_1.3.o : esercizio_1.3.cpp funzioni.h random.h
	$(CC) -c esercizio_1.3.cpp -o esercizio_1.3.o 
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o 
	
clean:
	rm *.o esercizio_1.3.exe seed.out
	rm -f *.dat