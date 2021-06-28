#include <iostream>
#include <fstream>
#include <cmath>
#include <vector> //container
#include <algorithm> //functions

using namespace std;

template <typename T> T computeMedia (const vector<T>& v ) {
	
	//computing media	
	T sum = 0;
	for (int i=0; i < v.size(); i++){
		sum += v[i];		
	}
	T media = sum/T(v.size());
	return media;	

};

template <typename T> vector<T> computeMediavett (const vector <vector <T> >& v ) {
	
	int dim = v[0].size();

	vector<T> mean(dim);
 	vector<T> sum(dim);
 	
 	for (int i = 0; i < dim; i++){
 	
 		sum[i] = 0.;
 	
 	}
 	
	for (int i = 0; i < dim; i++) {
	
		for (int j = 0; j < v.size(); j++){
		
		sum[i] += v[j][i] ;
		
		}
		 
	mean[i] = sum[i]/(v.size());
	
	}
	
	return mean;

};

template <typename T> T computeVarianza (const vector<T>& v ) {
	
	//computing media	
	T media = computeMedia<double> (v);

	//computing varianza
	T sum2 = 0;
	for (int i=0; i < v.size(); i++){
		sum2 += pow(v[i]-media, 2);
	}
	T varianza = sum2/T((v.size()-1));
	return varianza;

};

template <typename T> T computeVarianza (const vector<T>& v, double media ) {
	
	//computing varianza
	T sum2 = 0;
	for (int i=0; i < v.size(); i++){
		sum2 += pow(v[i]-media, 2);
	}
	T varianza = sum2/T((v.size()-1));
	return varianza;

};

template <typename T>  vector<T> computeDevstdvett (const vector< vector<T> >& v ) {

	int dim = v[0].size();

	vector<T> mean(dim);
	mean = computeMediavett(v);
	
	vector<T> sum2(dim);
	vector<T> var(dim);
	
	for (int i = 0; i < dim; i++){
 	
 		sum2[i] = 0.;
 	
 	}

	for (int i = 0; i < dim; i++){
		for (int j = 0; j < v.size(); j++){
			
			sum2[i] += pow(v[j][i] - mean[i], 2); 
	
		}

		var[i] = sqrt(sum2[i]/T((v.size() - 1)));

	}

	return var;

};

double computeCorrelation (const vector<double>& v, const vector<double>& w) {

	double dev_v = sqrt(computeVarianza(v));
	double dev_w = sqrt(computeVarianza(w));
	
	double med_v = computeMedia(v);
	double med_w = computeMedia(w);
	
	double sum = 0;
	
	for (int i = 0; i < v.size(); i++){

		sum += (v[i] - med_v)*(w[i] - med_w);
		
	}
	
	return sum/((v.size()) * dev_v * dev_w);

};

template <typename T> T computeCorrelationvett (const vector<vector<T> >& v) {

	vector <T> dev(2);
	dev = computeDevstdvett (v);
	
	vector <T> med(2);
	med = computeMediavett(v);
	
	double sum = 0;
	
	for (int j = 0; j < v.size(); j++ ){
		
		sum += (v[j][0] - med[0])*(v[j][1] - med[1]);
	
	}
	
	return sum/((v.size()) * dev[0] * dev[1]);

};

template <typename T> void sortArray (vector<T>& v ) {

	//sorting array
		for (int i=0; i<v.size(); i++) {
		for (int j=i+1; j<v.size(); j++) {
			if (v[j]<v[i]) {
				double tmp = v[j];
				v[j]=v[i];
				v[i]=tmp;
			}
		}
	}

};

template <typename T> T Mediana (vector<T> v) {

	//sorting array
	//creating copy
	vector<T> vcopy = v;	
	sort (vcopy.begin(), vcopy.end());
	
//finding mediana
	T mediana = 0;
	if (vcopy.size()%2==0){
		mediana = (vcopy[vcopy.size()/2]+vcopy[vcopy.size()/2-1])/2;
}
	else {
		mediana = vcopy[vcopy.size()/2];
}
	return mediana;

};

template <typename T> vector<T> ReadfromFile (char* filename, unsigned int size ) {

	vector<T> data;

	ifstream fin(filename);
		
	if (!fin){
	cout << "Error: cannot open file " << filename << endl;
	exit(2);
	} else {
	for (int i=0; i<size; i++){
		T temp ;
		fin >> temp; 
		data.push_back( temp );

	}
	fin.close();

	return data;
}
};

template <typename T> vector<T> ReadAll (string filename) {

	vector<T> data;

	ifstream fin(filename);
		
	if (!fin){
	cout << "Error: cannot open file " << filename << endl;
	exit(2);
	} else {
	while (!fin.eof()) {
		T temp ;
		fin >> temp; 
		data.push_back( temp );
		
	}
	fin.close();

	return data;
}
};

template <typename T> void Print (const char* filename ,const vector<T>& v) {
	
	ofstream outputFile(filename);
	if (!outputFile) {
		cerr << "cannot create file" << endl;
	} else {
	for (int z=0; z<v.size(); z++) {
		outputFile << v[z] << endl;
	}
	}
	outputFile.close();

};

template <typename T> void Print (const vector<T>& v ) {

	for (int i=0; i<v.size(); i++) {
		cout << "element number " << i+1 << " = " << v[i] << endl;
	}
};
