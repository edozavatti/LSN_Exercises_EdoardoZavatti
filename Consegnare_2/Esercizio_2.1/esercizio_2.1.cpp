#include "integralMC.h"
#include "random.h"
#include "FunzioneBase.h"
#include "funzioni.h"
#include <cmath>
#include <fstream>
#include <vector>

double error(double av, double av_2, int n) {
	
		if (n == 0) {
			return 0;
		}
		else {
			return sqrt((av_2 - pow(av, 2))/n);
		}
	}; 

int main() {

	//Funzione_2_1* f = new Funzione_2_1();
	Funzione_2_1_par* f_par = new Funzione_2_1_par();
	
	integralMC myMC;
	
	int M = 10000;
	int N = 100; //numero di blocchi
	int L = (int)M/N; //lanci per blocco
	
//	vector<double> media_medie;
//	vector<double> Err;
	
//	double sum_prog = 0.;
//	double sum_prog_2 = 0.;
	
	vector<double> media_medie_par;
	vector<double> Err_par;
	
	double sum_prog_par = 0.;
	double sum_prog_2_par = 0.;
	
	
	for (int k=0; k<N ; k++){
	
		cout << k<< endl;

		double I, I_par, errore_media, errore_media_par;		
		
//		double sum = 0.;
//		double A_i;
		
		double sum_par = 0.;
		double A_i_par;
		
		for (int i=0; i<L; i++){   //entro nel blocco k_i
		
			int punti = 1000;		
		//	I = myMC.integralAVE(0, 1, f, punti);
			I_par = myMC.integralAVE_par(0, 1, f_par, punti);
		//	cout << I << endl;
		//	sum += I;
			sum_par += M_PI/2. * I_par;
			
		} //ho finito i lanci nel blocco			

//		A_i = sum/L; //media del blocco i
		A_i_par = sum_par/L; //media del blocco i
		cout << A_i_par << endl;
//		sum_prog += A_i; //somma progressiva delle medie
//		sum_prog_2 += pow(A_i, 2); //somma progressiva delle medie al quadrato
//		errore_media = error(sum_prog/(k+1), sum_prog_2/(k+1), (k+1));  //errore statistico sulla media

		sum_prog_par += A_i_par; //somma progressiva delle medie
		sum_prog_2_par += pow(A_i_par, 2); //somma progressiva delle medie al quadrato
		errore_media_par = error(sum_prog_par/(k+1), sum_prog_2_par/(k+1), (k+1));  //errore statistico sulla media	
		
//		media_medie.push_back(sum_prog/(k+1) - 1);              //salvo la media di I su k+1 blocchi
//		Err.push_back(errore_media);         //salvo l'errore di I su k+1 blocchi		
	
		media_medie_par.push_back(sum_prog_par/(k+1) - 1);              //salvo la media di I su k+1 blocchi
		Err_par.push_back(errore_media_par);         //salvo l'errore di I su k+1 blocchi		
	
	
	}
			
//		Print("media_I.dat", media_medie);
//		Print("errore_I.dat", Err);
		Print("media_I_par.dat", media_medie_par);
		Print("errore_I_par.dat", Err_par);
	
	return 0;
	
}