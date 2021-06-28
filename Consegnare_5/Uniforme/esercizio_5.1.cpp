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
#include "random.h"
#include "FunzioneBase.h"
#include "Posizione.h"

using namespace std;

double error(double av, double av_2, int n) {
	
	if (n == 0) {
		return 0;
	}
	else {
		return sqrt((av_2 - pow(av, 2))/n);
	}
}
 
int main (int argc, char *argv[]){

   ofstream x_coord, y_coord, z_coord, r_ave, r_ave_error;

   Psi_100* f_1 = new Psi_100();
 //  Psi_210* f_1 = new Psi_210();

   double x, y, z;
   int Nstep = 1000000;
   int Nblocks = 100;
   int Nstep_block = int(Nstep/Nblocks);
   double passo = 2.45;

   x=3./2.;
   y=0.;
   z=0.;
   
   double ave = 0.;
   double ave2 = 0.;
   double sum_alpha = 0.;

   Posizione p(x, y, z);

   x_coord.open("output_x.dat",ios::app);
   y_coord.open("output_y.dat",ios::app);
   z_coord.open("output_z.dat",ios::app);
   r_ave.open("output_r_ave.dat",ios::app);
   r_ave_error.open("output_r_ave_error.dat",ios::app);

   for (int i = 1; i <= Nblocks; i++){
      double sum = 0.;
      //double sum2 = 0.;
      for(int step = 1; step <= Nstep_block; step++){
         p.Metropolis(passo, *f_1);
         if (step%500 == 0){
            x_coord << p.getX() << endl;
            y_coord << p.getY() << endl;
            z_coord << p.getZ() << endl;
         }

         double r = p.getR();
         double a = p.getAlpha();
         sum += r; 
         sum_alpha += a;
         
      }
      ave += sum/Nstep_block;
      ave2 += pow(sum/Nstep_block, 2);

      r_ave << ave/double(i) << endl;
      r_ave_error << error(ave/(double(i)), ave2/double(i), i) << endl;

   }

   x_coord.close();
   y_coord.close();
   z_coord.close();
   r_ave.close();
   r_ave_error.close();

   cout << "Acceptance = " << sum_alpha/double(Nstep) << endl;


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
