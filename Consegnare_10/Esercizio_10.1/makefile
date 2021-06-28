CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio_quad_10.1.exe : esercizio_quad_10.1.o random.o
	$(CC) random.o esercizio_quad_10.1.o -o esercizio_quad_10.1.exe
esercizio_10.1.exe : esercizio_10.1.o random.o
	$(CC) random.o esercizio_10.1.o -o esercizio_10.1.exe
esercizio_10.1.o : esercizio_10.1.cpp random.h funzioni.h
	$(CC) -c esercizio_10.1.cpp -o esercizio_10.1.o $(CFLAGS)
esercizio_quad_10.1.o : esercizio_quad_10.1.cpp random.h funzioni.h
	$(CC) -c esercizio_quad_10.1.cpp -o esercizio_quad_10.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o esercizio_10.1.exe 
	rm *.dat

