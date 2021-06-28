/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

double error(double av, double av_2, int n) {
	
	if (n == 0) {
		return 0;
	}
	else {
		return sqrt((av_2 - pow(av, 2))/n);
	}
}

int main(){ 

  ofstream Epot_ave, Ekin_ave, Etot_ave, Temp_ave;
  ofstream Epot_error, Ekin_error, Etot_error, Temp_error;

  Input();             //Inizialization
  int nconf = 1;
  
  double sum_pot=0, sum_kin=0, sum_temp=0, sum_etot=0;
  double sum2_pot=0, sum2_kin=0, sum2_temp=0, sum2_etot=0;

  Epot_ave.open("Esercizio_4.3/output_epot_avegas_new.dat",ios::app);
  Ekin_ave.open("Esercizio_4.3/output_ekin_avegas_new.dat",ios::app);
  Temp_ave.open("Esercizio_4.3/output_temp_avegas_new.dat",ios::app);
  Etot_ave.open("Esercizio_4.3/output_etot_avegas_new.dat",ios::app);

  Epot_error.open("Esercizio_4.3/output_epot_errorgas_new.dat",ios::app);
  Ekin_error.open("Esercizio_4.3/output_ekin_errorgas_new.dat",ios::app);
  Temp_error.open("Esercizio_4.3/output_temp_errorgas_new.dat",ios::app);
  Etot_error.open("Esercizio_4.3/output_etot_errorgas_new.dat",ios::app);

  for(int iblock=1; iblock <= nblocks; ++iblock){ // ciclo sui blocchi
    double pot_ave = 0;
    double kin_ave = 0;
    double etot_ave = 0;
    double temp_ave = 0;
    if(iblock%iprint == 0) cout << "Number of blocks: " << iblock << endl; 
    double measures = 0;
    for(int istep=1; istep<=nstep_block; ++istep){ //entro nel blocco i
      Move();           //Move particles with Verlet algorithm
      if(istep%10 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        
        pot_ave += stima_pot;  //Potential energy per particle
        kin_ave += stima_kin; //Kinetic energy per particle
        temp_ave += stima_temp; //Temperature
        etot_ave += stima_etot; //Total energy per particle
        
        measures++;
        nconf += 1;
      }
    }  
    double A_i = pot_ave/double(measures); //energia potenziale media del blocco i
    double B_i = kin_ave/double(measures);
    double C_i = temp_ave/double(measures);
    double D_i = etot_ave/double(measures);
    
    sum_pot += A_i; //sommo le energie di ogni blocco
    sum_kin += B_i;
    sum_temp += C_i;
    sum_etot += D_i;
    
    sum2_pot += pow(A_i, 2);
    sum2_kin += pow(B_i, 2);
    sum2_temp += pow(C_i, 2);
    sum2_etot += pow(D_i, 2);
     
    Epot_ave << sum_pot/double(iblock)  << endl;
    Ekin_ave << sum_kin/double(iblock)  << endl;
    Temp_ave << sum_temp/double(iblock) << endl;
    Etot_ave << sum_etot/double(iblock) << endl;
    
    Epot_error << error(sum_pot/double(iblock), sum2_pot/double(iblock), iblock) << endl;
    Ekin_error << error(sum_kin/double(iblock), sum2_kin/double(iblock), iblock) << endl;
    Temp_error << error(sum_temp/double(iblock), sum2_temp/double(iblock), iblock) << endl;
    Etot_error << error(sum_etot/double(iblock), sum2_etot/double(iblock), iblock) << endl;

  }

  Epot_ave.close();
  Ekin_ave.close();
  Temp_ave.close();
  Etot_ave.close(); 

  Epot_error.close();
  Ekin_error.close();
  Temp_error.close();
  Etot_error.close();

  ConfFinal();         //Write final configuration to restart
  
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation

  //char b;
  bool a;
  /*cout << "Is it a restart? [y/n]" << endl;
  cin >> b ;
  if (b == 'y'){
    a = true;
    cout << "Restart" << endl;
  }
  else if (b == 'n') {
    a = false;
    cout << "Start from scratch" << endl;
  }*/

  ifstream ReadInput,ReadConf, ReadOldConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl; //densità numerica non di massa
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nblocks;
  ReadInput >> a;

  nstep_block = int(nstep)/int(nblocks);

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "number of blocks = " << nblocks << endl << endl;
  cout << a << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  if (a == true){ //Restart case
    //Read configuration a time t - dt
    cout << "Read pre-initial configuration from file old_config.0 " << endl << endl;
    ReadOldConf.open("old_config.0");
    for (int i=0; i<npart; ++i){
      ReadOldConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadOldConf.close();

    Move(); //move the particles via Verlet algorithm
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = (x[i] - xold[i])/delta;
      vy[i] = (y[i] - yold[i])/delta;
      vz[i] = (z[i] - zold[i])/delta;
      
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }

    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i) sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    sumv2 /= (double)npart;

  //  double T = 2.0/3.0 * sumv2;
  //  fs = et/T;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor, fissa la temperatura iniziale a quella target
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }

  }

  else {
    //Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];  // sottraggo il sumv alla velocità di ogni particella per avere che il centro di massa 
     vy[i] = vy[i] - sumv[1];  // fermo, senza velocità di deriva
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor, fissa la temperatura iniziale a quella target
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }

  }
   return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  /*ofstream Epot, Ekin, Etot, Temp;
  
  Epot.open("Esercizio_4.3/output_epotgas_new.dat",ios::app);
  Ekin.open("Esercizio_4.3/output_ekingas_new.dat",ios::app);
  Temp.open("Esercizio_4.3/output_tempgas_new.dat",ios::app);
  Etot.open("Esercizio_4.3/output_etotgas_new.dat",ios::app);*/

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle

  /*Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();*/

  return;
}


void ConfFinal(void){ //Write final configuration

  ofstream WriteConf, WriteOldConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  
  cout << "Print pre-final configuration to file config.final " << endl << endl;
  WriteOldConf.open("old_config.final");

  for (int i=0; i<npart; ++i){
    WriteOldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteOldConf.close();
  
  return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}


double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
