#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "City.h"

using namespace std;

City :: City(double x, double y, int n){
    m_x = x;
    m_y = y;
    m_num = n;
}

City :: ~City(){}

double City :: GetX(){
   return m_x;
}

double City :: GetY(){
   return m_y;
}

double City :: Getnum(){
   return m_num;
}

void City :: SetX(double x){
   m_x = x;
}

void City :: SetY(double y){
   m_y = y;
}

void City :: Setnum(int num){
   m_num = num;
}