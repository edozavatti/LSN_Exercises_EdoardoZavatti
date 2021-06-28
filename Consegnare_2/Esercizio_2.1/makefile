CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio_2.1.exe : esercizio_2.1.o random.o integralMC.o
	$(CC) integralMC.o random.o esercizio_2.1.o -o esercizio_2.1.exe
esercizio_2.1.o : esercizio_2.1.cpp funzioni.h integralMC.h FunzioneBase.h
	$(CC) -c esercizio_2.1.cpp -o esercizio_2.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
integralMC.o : integralMc.cpp integralMC.h FunzioneBase.h
	$(CC) -c integralMC.cpp -o integralMC.o $(CFLAGS)	

clean:
	rm *.o esercizio_2.1.exe 
	rm media_I_par.dat
	rm errore_I_par.dat