/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie, iw, igofr;
double ptail, vtail;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
double err_pot, err_kin, err_etot, err_temp, err_gdir, err_pres;

// averages
double acc,att;
double blk_av[m_props];
double glob_av[m_props], glob_av2[m_props];
double blk_norm, accepted, attempted;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double walker[m_props];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblk, n_bins;
double delta, bin_size;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double, double, int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
