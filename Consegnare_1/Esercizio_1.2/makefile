CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio_1.2.exe : esercizio_1.2.o random.o
	g++ -std=c++11 random.o esercizio_1.2.o -o esercizio_1.2.exe
esercizio_1.2.o : esercizio_1.2.cpp funzioni.h
	g++ -std=c++11 -c esercizio_1.2.cpp -o esercizio_1.2.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -std=c++11 -c random.cpp -o random.o $(CFLAGS)
	
clean:
	rm *.o esercizio_1.2.exe seed.out