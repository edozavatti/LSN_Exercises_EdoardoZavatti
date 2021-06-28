/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include "City.h"
#include "random.h"
#include <algorithm>

using namespace std;

Random my_rand;
int N_pop = 600;
int N = 32;
vector<vector<double> > Matrix;
//vector<vector<double> > Matrix_2;

/*void Print(vector< vector< int > > a){
    int N = a.size();
    for(int i = 0; i < N; i++){
        int M = a[i].size();
        for(int j = 0; j < M; j++){
            cout << a[i][j] << " " ;
        }
        //cout << L2(a[i]) << endl;
        cout << endl;
    }
}*/

template <typename T> void Print(vector< vector< T > > a){
    int N = a.size();
    for(int i = 0; i < N; i++){
        int M = a[i].size();
        for(int j = 0; j < M; j++){
            cout << a[i][j] << " " ;
        }
        cout << endl;
    }
}

template <typename T> void Print(vector< T > a){
    int N = a.size();
    for(int i = 0; i < N; i++){    
        cout << a[i] << " " ;
    }
    cout << endl;
}

vector<double> Generate_Pos_circ(int N){
	vector<double> theta;
	for(int i=0; i<N; i++){
		theta.push_back(my_rand.Rannyu(0., 2*M_PI));
        //theta.push_back(2*M_PI/double(N) * i);
	}
    return theta;
}

vector<vector<double>> Generate_Pos_square(int N){
	vector <double> x;
	vector <double> y;
	vector<vector<double>> result;
	for(int i=0; i<N; i++){
		x.push_back(my_rand.Rannyu());
		y.push_back(my_rand.Rannyu());
	}
	result.push_back(x);
	result.push_back(y);
    return result;
}

vector<vector<double>> Generate_Dist_square(vector<vector<double>> pos, int N){
	
	vector<double> v;
	vector<vector<double>> result;
	double Mij;
		
	for(int i=0; i<N; i++){
		v.clear();
		for(int j=0; j<N; j++){
			Mij = pow((pos[0][i] - pos[0][j]), 2) + pow((pos[1][i] - pos[1][j]), 2);
			v.push_back(Mij);
		}
		
		result.push_back(v);
	}
	
    return result;
}

vector<vector<double>> Generate_Dist_circ(vector<double> theta, int N){
    
    vector<vector<double>> result;
    //GENERAZIONE MATRICE DISTANZE
	for(int i = 0; i < N; i++){
        vector<double> my_vector;
        double x_i = cos(theta[i]);
        double y_i = sin(theta[i]);
        for(int j = 0; j < N; j++){
            double x_j = cos(theta[j]);
            double y_j = sin(theta[j]);
            double M_ij = (pow((x_i - x_j), 2) + pow((y_i - y_j), 2));
            my_vector.push_back(M_ij);
        }

        result.push_back(my_vector);
    
    }
    //Print(Matrix);
    return result;
}

void Check(vector<int> a){
	int N = a.size();
	int colonna;
	vector<vector<int> > matrix;
	for(int i = 0; i<N; i++){
		//creo matrice di zeri
    	vector<int> myvector;
    	for(int j = 0; j<N; j++){
        	int tempVal = 0;
        	myvector.push_back(tempVal);
    	}
    	matrix.push_back(myvector);
	}
    if(a[0] != 1){
        cout << "Check fallito: la prima città non è 1" << endl;
    }
	for(int i = 0; i < N; i++){
        colonna = a[i]-1;
		for(int j = 0; j < N; j++){
			if(j != colonna){
				matrix[i][j] = 0;
			}
			else{
				matrix[i][j] = 1;
			}
		}
	}

    //Print(matrix);

	int sum_col, sum_rig;

	for(int i = 0; i < N; i++){
		sum_col = 0;
		for(int j = 0; j < N; j++){
			sum_col += matrix[j][i];	
		}
		if(sum_col != 1){
			cout << "check fallito" << endl;
		}
	}
	for(int i = 0; i < N; i++){
		sum_rig = 0;
		for(int j = 0; j < N; j++){
			sum_rig += matrix[i][j];	
		}
		if(sum_rig != 1){
			cout << "check fallito" << endl;
		}
	}

}

double L2(vector<int> C){
    double sum = 0.;
    int N = C.size();
    int a, b;
    for(int i = 0; i < N-1; i++){  
        a = C[i] - 1;
		b = C[i+1] -1;  
        sum += Matrix[a][b];    
	}
    sum += Matrix[a][0];
    
    return sum;
}

void Disorder(vector<int> &a){
    int N = a.size();
	int i, j;
	int n=0;
	while(n<16){
		do{
			i = int(my_rand.Rannyu(1, N));
			j = int(my_rand.Rannyu(1, N));
		} while(i == j);
		n++;
		//cout << "eseguo swap fra " << i << " e " << j << endl;
        iter_swap(a.begin() + i,a.begin() + j); 
    }
}

vector<int> Generate_pop(int N){
    vector<int> a;
	for(int i = 0; i < N; i++){
        a.push_back(i + 1);
    }
    Disorder(a);
    return a;
}

template <typename T> void Fitness_sort (vector<vector<T> >& v ) {
	//sorting array
	for (int i=0; i<v.size(); i++) {
		for (int j=i+1; j<v.size(); j++) {
			if (L2(v[j]) < L2(v[i])) {
				v[i].swap(v[j]);
			}
		}
	}
}

int Selection(vector<vector<int> > P){
    int N = P.size();
    double p = 3.;
    double r = my_rand.Rannyu(0., 1.);
    int a = int(N * pow(r, p));
    return a;
}

int Selection(vector<vector<int> > P, int b){
    int N = P.size();
    double p = 3.;
    int a;
    do{
        double r = my_rand.Rannyu();
        a = int(N * pow(r, p));    
    }while(a == b);
    return a;
}

/*vector<int> Zeros(vector<int> a, int n){
    for(int i = 0; i < n; i++){
        a.push_back(0);
    }
}*/

void Pair_Permutation(vector<int> &a){
    int N = a.size();
    int i = int(my_rand.Rannyu(1, N));
	int j;
	do{
		j = int(my_rand.Rannyu(1, N));
	} while(i == j);
    //cout << "eseguo permutazione di coppia fra " << a[i] << " e " << a[j] << endl;
    iter_swap(a.begin() + i,a.begin() + j);     
}

void M_Permutation(vector<int> &a){
    int N = a.size();
	int i, j;
	int n=0;
    int M = my_rand.Rannyu(2, int(N)/int(2));
    //cout << "eseguo " << M << " permutazioni" << endl;
	while(n<M){
		do{
			i = int(my_rand.Rannyu(1, N));
			j = int(my_rand.Rannyu(1, N));
		} while(i == j);
		n++;
		//cout << "eseguo permutazione di coppia fra " << a[i] << " e " << a[j] << endl;
        iter_swap(a.begin() + i,a.begin() + j); 
    }
    Check(a);
}

vector<int> Shift(vector<int> a){
    int N = a.size();
    int start = my_rand.Rannyu(1, int(N)/int(2));
    int n_shift = my_rand.Rannyu(1, N-start);
    int shift = my_rand.Rannyu(1, N-n_shift-start);
    //cout << endl;
    //cout << "eseguo shift di " << shift << " posizioni, a partire da " << a[start] << " di " << n_shift << " elementi" << endl;
    //Print(a);
    vector<int> my_vector;
    for(int i = 0; i < start; i++){//riempio il vettore temporaneo con elementi che non shiftano
        my_vector.push_back(a[i]);
    }
    for(int i = 0; i < shift; i++){
        my_vector.push_back(a[N - shift + i]);
    }
    for(int i = 0; i < N - start - shift; i++){
        my_vector.push_back(a[start + i]);
    }
    Check(my_vector);
    return my_vector;
}

void Inversion(vector<int> &v){
	int N = v.size();
	int m = my_rand.Rannyu(1, N/int(2));
    m = 2*m;
    int a;
	//cout << "m " << m << endl;
    do{
        a = my_rand.Rannyu(1, N-1);
    } while(a > m);
	//cout << "eseguo inversione da " << v[a] << " a " << v[m] << endl;
    int intervallo = m - a;
	for(int i=0; i<intervallo/int(2); i++){
		iter_swap(v.begin() + a + i, v.begin() + m - i);
	}
}

vector<vector<int> > Populate(vector<vector<int> > P, int k){
    vector<vector<int> > New_population;
    int N = P[0].size();
    int N_pop = int(P.size())/int(2);

    for(int j = 0; j < N_pop; j++){
        vector<int> Figlio_1, Figlio_2;
        //Seleziono i genitori
        int a = Selection(P);
        int b = Selection(P, a);
        vector<int> Genitore_1, Genitore_2;
        for(int i = 0; i < N; i++){
            Genitore_1.push_back(P[a][i]);
            Genitore_2.push_back(P[b][i]);
        }
       
        if(my_rand.Rannyu() > 0.25){ //chiamo un crossover con probabilità del 50%
            //Seleziono il punto in cui fare crossover
            int c = int(my_rand.Rannyu(1, N));
            //Copio i cromosomi fino al punto selezionato
            for(int i = 0; i < c; i++){
                Figlio_1.push_back(Genitore_1[i]);
                Figlio_2.push_back(Genitore_2[i]);
            }
            //Eseguo il crossover
            //cout << "crossover a " << c + 1 << endl;
            //cout << "Figli : " << endl;
            for(int i = 0; i < N; i++){
                int check = Genitore_2[i];
                //cout << "elemento " << i << " da scambiare = " << check << endl;;
                bool cross = true;
                for(int j = 0; j < c; j++){
                    if (check == Genitore_1[j]){
                        cross = false;
                    }
                }
                //cout << cross << endl;
                if(cross != false){
                    Figlio_1.push_back(check);
                    //cout << " scambiato " << check << endl;
                }
            }
            
            for(int i = 0; i < N; i++){
                int check = Genitore_1[i];
                //cout << "elemento " << i << " da scambiare = " << check << endl;;
                bool cross = true;
                for(int j = 0; j < c; j++){
                    if (check == Genitore_2[j]){
                        cross = false;
                    }
                }   
                //cout << cross << endl;
                if(cross != false){
                    Figlio_2.push_back(check);
                    //cout << " scambiato " << check << endl;
                }
            }      
        }
        else { //non avviene crossover, i figli sono uguali ai genitori
            Figlio_1 = Genitore_1;
            Figlio_2 = Genitore_2;
        }
        
        //Check(Figlio_1);
        //Check(Figlio_2);

        /*cout << "Figli prima delle mutazioni: " << endl;
        Print(Figlio_1);
        Print(Figlio_2);
        cout << endl;*/
        
        if(k != 600){
        
        //MUTAZIONI
        //permutazioni singole
        if(my_rand.Rannyu() > 0.7){
            Pair_Permutation(Figlio_1);
            //Check(Figlio_1);
        }
        if(my_rand.Rannyu() > 0.7){
            Pair_Permutation(Figlio_2);
            //Check(Figlio_2);
        }
        //permutazioni multiple
        if(my_rand.Rannyu() > 0.7){
            M_Permutation(Figlio_1);
            //Check(Figlio_1);
        }
        if(my_rand.Rannyu() > 0.7){
            M_Permutation(Figlio_2);
            //Check(Figlio_2);
        }
        
        //shift
        if(my_rand.Rannyu() > 0.7){
            //Print(Figlio_1);
            Genitore_1 = Shift(Figlio_1);
            //Print(Figlio_1);
        }
        if(my_rand.Rannyu() > 0.7){
            Genitore_2 = Shift(Figlio_2);
        }
        
        //inversion
        if(my_rand.Rannyu() > 0.7){
            //Print(Figlio_1);
            Inversion(Figlio_1);
            //Print(Figlio_1);
            //Check(Figlio_1);
        }
        if(my_rand.Rannyu() > 0.7){
            Inversion(Figlio_2);
            //Check(Figlio_2);
        }
        
        }
        
        else cout << "ultima generazione: non faccio mutazioni" << endl;
        
        Check(Figlio_1);
        Check(Figlio_2); 

        New_population.push_back(Figlio_1);
        New_population.push_back(Figlio_2);
        
    }
    return New_population;
}

int main (int argc, char *argv[]){

	ofstream stampa, media, media_best;

   int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
  		Primes >> p1 >> p2 ;
	} else std::cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
  		while ( !input.eof() ){
     		input >> property;
     		if( property == "RANDOMSEED" ){
        		input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        		my_rand.SetRandom(seed,p1,p2);
     		}
  		}
  		input.close();
	} else std::cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	my_rand.SaveSeed();

	
    vector< vector<int> > Population;

    //Genero posizioni dentro un quadrato
    vector<vector<double> > square = Generate_Pos_square(N);
    Print(square);
    cout << endl;
    //Genero matrice distanze
    Matrix = Generate_Dist_square(square, N);
    Print(Matrix);

    for(int i = 0; i < N_pop; i++){
        vector<int> my_vector;
        my_vector = Generate_pop(N);
        Check(my_vector);
        Population.push_back(my_vector);
        //cout << "calcolo percorso " << i << endl;
        //cout << "route = " << L2(my_vector) << endl;
    }
	
    stampa.open("Route_quad_0", ios::app);
    for(int i = 0; i < N; i++){
        double x_i = square[0][i];
        double y_i = square[1][i];
        stampa << x_i << " " << y_i << " " << endl;
    }
    stampa.close();

    /*cout << "popolazione dopo il fitness sort" << endl;
    cout << endl;
    Fitness_sort(Population);
    Print(Population);
    cout << endl;
    int sel = Selection(Population);
    int sel_2 = Selection(Population, sel);
    cout << "Primo cromosoma scelto: " << sel << endl;
    cout << "Secondo cromosoma scelto: " << sel_2 << endl;
    */
    cout << "Popolazione di partenza" << endl;
    Print(Population);
    cout << endl;
        
   vector<vector<int> > New_Population;
   int n=1;
   while(n < 501){
       string filename = "Route_quad_" + to_string(n);
        New_Population = Populate(Population, n);
        Fitness_sort(New_Population);
        double route = 0.;
        int num = 0;
        for(int i = 0; i < int(N_pop)/int(2); i++) {
            route += L2(New_Population[i]);
            num ++;
        }
        media.open("Media_quad",ios::app);
        media_best.open("Media_quad_migliore", ios::app);
        route = route/double(num);
        media << route << endl;
        media_best << L2(New_Population[0]) << endl;
        media.close();
        media_best.close();
        if(n%100 == 0){
            stampa.open(filename, ios::app);
            for(int i = 0; i < N; i++){
                double x_i = square[0][Population[0][i]-1];
                double y_i = square[1][Population[0][i]-1];
                stampa << x_i << " " << y_i << " " << endl;
            }
            //cout << Population[0][0] << endl;
            stampa << square[0][0] << " " << square[1][0] << endl;
            stampa.close();
        }
        n++;
        cout << "Generazione " << n << endl;
        stampa.close();
        Population = New_Population;

    }
\
    cout << "Popolazione di arrivo" << endl;
    Print(New_Population);
    cout << endl;
    


    return 0;
} 







/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/