/*
 * Liear oscillator
 * 2.12.2015
 * A.M.Hamada
 */

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void Analaytic(const double dt1, const double dt2, const double N1, const double N2, double* a1, double* a2);
void Write (const double N1, const double N2, double* a1, double* a2, double* f1, double* f2, double* h1, double* h2, const double dt1, const double dt2);
void Forward (const double dt1, const double dt2, const double N1, const double N2, double* f1, double* f2);
void Backward (const double dt1, const double dt2, const double N1, const double N2, double* h1, double* h2);

int main (void){
    const double dt1 = M_PI / 10 ;
    const double dt2 = M_PI / 100 ;
    const int N1 = (20 * M_PI) / dt1;
    const int N2 = (20 * M_PI) / dt2;
    double a1[N1], f1[N1], h1[N1];
    double a2[N2], f2[N2], h2[N2];
    
      Analaytic (dt1, dt2, N1, N2, a1, a2);
      Forward (dt1, dt2, N1, N2, f1, f2);
      Backward (dt1, dt2, N1, N2, h1, h2); 
      Write (N1, N2, a1, a2, f1, f2, h1, h2, dt1, dt2);
      
  return 0;
}

void Analaytic(const double dt1, const double dt2, const double N1, const double N2, double* a1, double* a2){
      for(int i=0; i<N1; i++) a1[i] = cos (i*dt1);
      for(int i=0; i<N2; i++) a2[i] = cos (i*dt2);
}
void Write (const double N1, const double N2, double* a1, double* a2, double* f1, double* f2, double* b1, double* b2, const double dt1, const double dt2){
      ofstream print_1("Euler_Os_dt1.txt");
      for (int i=0; i<N1; i++) print_1 << dt1 * i << "\t" << a1[i] << "\t" << f1[i] << "\t" << b1[i] << endl;
      
      ofstream print_2("Euler_Os_dt2.txt");
      for (int i=0; i<N2; i++) print_2 << dt2 * i << "\t" << a2[i] << "\t" << f2[i] << "\t" << b2[i] << endl;
}
void Forward (const double dt1, const double dt2, const double N1, const double N2, double* f1, double* f2){
  double x[2];
  x[0] = 1;
  x[1] = 0;
  f1[0] = x[0];
  f2[0] = x[0];
  for (int i = 0 ; i< N1; i++){
    double temp = x[0];
    x[0] = x[0] + x[1] * dt1;
    x[1] = x[1] - temp * dt1;
    f1[i] = x[0];
  }
  x[0] = 1;
  x[1] = 0;
  for (int i = 0; i<N2; i++){
    double temp = x[0];
    x[0] = x[0] + x[1] * dt2;
    x[1] = x[1] - temp * dt2;
    f2[i] = x[0];
  }
} 
void Backward (const double dt1, const double dt2, const double N1, const double N2, double* h1, double* h2){
  double x[2];
  x[0] = 1;
  x[1] = 0;
  h1[0] = x[0];
  h2[0] = x[0];
  for (int i = 0 ; i< N1; i++){
    double temp = x[0];
    x[0] = x[0] + x[1] * dt1;
    x[0] = x[0] / (1+ pow(dt1,2));
    x[1] = x[1] - temp * dt1;
    x[1] = x[1] / (1+ pow(dt1,2));
    h1[i] = x[0];
  }
  x[0] = 1;
  x[1] = 0;
  for (int i = 0; i<N2; i++){
    double temp = x[0];
    x[0] = x[0] + x[1] * dt2;
    x[0] = x[0] / (1+ pow(dt2,2));
    x[1] = x[1] - temp * dt2;
    x[1] = x[1] / (1+ pow(dt2,2));
    h2[i] = x[0];
  }
}
/*
 *gnuplot> set term png
 * gnuplot> set out "plot.png"
 * gnuplot> plot "noisy.txt", "filtered.txt" 
 */