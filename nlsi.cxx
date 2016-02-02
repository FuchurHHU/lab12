#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);
void Lstep(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N);
void Nstep(cmplx* psi0, cmplx* psiS, const double dt, const int N);
//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psiS = new cmplx[Nx]; // Psi tilde als Zwischenergebnis
	cmplx* h;
	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-1; j++) {
		  
		  Lstep(psiS, psi0, dt, dx, Nx); // 1. Schritt der linearen Gleichung (vgl Mitschrift); Jetzt ist das alte Psi in psi0 und das neue in psiS gespeichert
		  
		  Nstep(psi0, psiS, dt,Nx); // 2. Schritt ""0; wir errechnen jetzt das neue psi und nutzen dabei psiS; dazu überschreiben wir direkt psi0, was dann im nächsten Schritt direkt wieder genutzt wird
// 		  h = psi0;
// 		  psi0 = psiS;
// 		  psiS = h; 
		  
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}
	delete[] psi0;
	delete[] psiS;
	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] =2* f/cosh(eta * x); // Anfangshöhe ändert Verlauf der Fkt (Versuche mit * 2 weg)
	}
}
void Lstep(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N) //aus lab 11 geklaut
{

  cmplx* d=new cmplx[N];
  cmplx* u=new cmplx[N];
  cmplx* l=new cmplx[N];

  for(int i=0;i<N;i++) d[i] = cmplx(1.0,-2.0*dt/(dx*dx));
  for(int i=0;i<N;i++) u[i] = cmplx(0.0,dt/(dx*dx));
  for(int i=0;i<N;i++) l[i] = cmplx(0.0,dt/(dx*dx));
  
  for(int i=1;i<N;i++){
  
    d[i]  -= u[i-1]*l[i]/d[i-1];
    f0[i] -= f0[i-1]*l[i]/d[i-1];
  }
  
  f1[N-1] = f0[N-1]/d[N-1];  //f0 und d sind schon die geschlaengelten
  
  for(int i=N-2;i>=0;i--){
  
    f1[i] = (f0[i] - u[i]*f1[i+1]) / d[i];    
  }
}
  void Nstep(cmplx* psi0, cmplx* psiS, const double dt, const int N){
  
    for (int i=0; i<N;i++){
    cmplx v = cmplx(0,-norm(psiS[i])*dt); 
    psi0[i] = psiS[i] * exp(v);
    }
  }
  