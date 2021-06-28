#ifndef __FunzioneBase_h__
#define __FunzioneBase_h__

#include <iostream>
#include <cmath>
#include <math.h>

using namespace std;

class FunzioneBase {
	
	public:
	
	virtual double Eval(double x) const = 0;
 
};

class FunzioneVettoriale {
	
	public:
	
	virtual double Eval(double x, double y, double z) const = 0;
 
};

class Parabola : public FunzioneBase {

	public:
	
	Parabola();	
	Parabola(double a, double b, double c);
	~Parabola();

	virtual double Eval(double x) const { return m_a*x*x + m_b*x + m_c; }
	void SetA (double a) { m_a = a; }
	double GetA() const { return m_a; }	
	void SetB (double b) { m_b = b; }
	double GetB() const { return m_b; }
	void SetC (double c) { m_c = c; }
	double GetC() const { return m_c; }

	private:

	double m_a, m_b, m_c;

};
class Seno : public FunzioneBase {

	public :
	
	Seno(){};	
	~Seno(){};

	double Eval(double x) const { return sin(x); }
};

class Id : public FunzioneBase {

	public :
	
	Id(){};	
	~Id(){};

	double Eval(double x) const { return 1; }
};

class Funzione_2_1 : public FunzioneBase {

	public :
	
	Funzione_2_1(){};	
	~Funzione_2_1(){};

	double Eval(double x) const { return 0.5 * M_PI * cos(0.5 * M_PI * x); }
};

class Funzione_2_1_par : public FunzioneBase {

	public :
	
	Funzione_2_1_par(){};	
	~Funzione_2_1_par(){};

	double Eval(double x) const { return cos(0.5 * M_PI * x)/(2*(1-x)); }
};

class Psi_100 : public FunzioneVettoriale {

	public :
	
	Psi_100(){};	
	~Psi_100(){};

	double Eval(double x, double y, double z) const { return pow(1./sqrt(M_PI) * exp(- sqrt((x*x + y*y + z*z))) , 2); }
};

class Psi_210 : public FunzioneVettoriale {

	public :
	
	Psi_210(){};	
	~Psi_210(){};

	double Eval(double x, double y, double z) const { return pow(1./8. * sqrt(2.)/sqrt(M_PI) * z * exp(- 1./2. * sqrt((x*x + y*y + z*z))) , 2); }
};

class Psi_sigma_mu : public FunzioneVettoriale {

	public :
	
	Psi_sigma_mu(){};	
	~Psi_sigma_mu(){};

	double Eval(double x, double mu, double sigma) const { 
		
		double f = exp(-pow((x-mu), 2)/(2. * pow(sigma, 2))) + exp(-pow((x+mu), 2)/(2. * pow(sigma, 2)));
		
		return pow(f, 2);
	
	}
};

class H_Psi : public FunzioneVettoriale {

	public :
	
	H_Psi(FunzioneBase *potential, double m, double h){
		//m_V = potential;
		m_m = m;
		m_hbar = h;
	}	
	~H_Psi(){};

	double Eval(double x, double mu, double sigma) const { 
		
		double V = pow(x, 4) - 2.5 * pow(x, 2);
		double kinetic = exp(-pow(x-mu, 2)/(2.*pow(sigma, 2))) * (- pow(sigma, 2) + pow(x-mu, 2))/(pow(sigma, 4)) + exp(-pow(x+mu, 2)/(2.*pow(sigma, 2))) * (pow(x+mu, 2) - pow(sigma, 2))/(pow(sigma, 4));
		double f = exp(-pow(x-mu, 2)/(2.*pow(sigma, 2))) + exp(-pow(x+mu, 2)/(2.*pow(sigma, 2)));

		return (-pow(m_hbar, 2)/(2. * m_m) * kinetic + V * f )/f ;
	
	}

	private :

	FunzioneBase* m_V;
	double m_m, m_hbar;

};

class Potenziale_doppia_buca : public FunzioneBase {

	public :
	
	Potenziale_doppia_buca(){};	
	~Potenziale_doppia_buca(){};

	double Eval(double x) const { return pow(x, 4) - 2.5 * pow(x, 2); }
};

#endif //__FunzioneBase_h__
