CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: main.exe main_bis.exe

main_bis.exe : main_bis.o random.o
	$(CC) random.o main_bis.o -o main_bis.exe
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp funzioni.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
main_bis.o : main_bis.cpp funzioni.h
	$(CC) -c main_bis.cpp -o main_bis.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
	rm *.dat
