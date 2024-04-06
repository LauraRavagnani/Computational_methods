#include<iostream>
#include<cmath>
#include<complex>
#include<fstream>

using namespace std;

void risolvi_matrice (int N, complex<double> *psi, complex <double> *F, complex <double> **M){		//psi = a  //F = b
	complex<double>* alpha = (complex<double>*) malloc (N * sizeof(complex<double>));
	complex<double>* beta = (complex<double>*) malloc (N * sizeof(complex<double>));

	alpha[0] = - M[0][0] / M[0][1];
	beta[0] = F[0] / M[0][1];
	alpha[N-1] = 0.;
	beta[N-1] = 0.;

	for (int i = 1; i < N; i++){
		alpha[i] = - (M[i][i-1] / (M[i][i+1] * alpha[i-1]) + M[i][i] / M[i][i+1]);
		beta[i] = F[i] / M[i][i+1] + M[i][i-1] / M[i][i+1] * beta[i-1] / alpha[i-1];
	}

	psi[N-1] = (F[N-1] / M[N-1][N-2] + beta[N-2] / alpha[N-2]) / (1. / alpha[N-2] + M[N-1][N-1] / M[N-1][N-2]);

	for (int i = N - 2; i >= 0; i--){
		psi[i] = (psi[i+1] - beta[i]) / alpha[i];
	}

	return;
}

double V (double x){		//funzione potenziale
	double V0 = 1.7;
	double a = 250., b = 260.;
	double V;
	if (x > a && x < b)
		V = V0;
	else
		V = 0;
	
	return V;
}

void norma (complex<double>* psi, double h, double* x){
	double norm = 0.;
	double Nx = 1000;
	double N = Nx - 2;

	for (int i = 0; i < N - 1; i++){
		norm += abs(psi[i] * conj(psi[i]) * h);
	}

	for (int i = 0; i < N -1; i++){
		psi[i] /= sqrt(norm);
	}

	return;
}

int main(){
	double L = 500, x0 = 200, q = 2, sigma = 20, a = 250, b = 260, Nx = 1000;
	int tstep = 10000;
	double V0 = 1.7, dt = 0.1;
	int N = Nx - 2;
	double h = L / (Nx-1);
	double* x = (double*) malloc(Nx * sizeof(double));
	complex<double>* F = (complex<double>*) malloc(N * sizeof(complex<double>));
	complex<double>* psi = (complex<double>*) malloc(N * sizeof(complex<double>));
	complex<double>* m = (complex<double>*) malloc (N * N * sizeof(complex<double>));
	complex<double>** M = (complex<double>**) malloc(N * sizeof(complex<double>));

	ofstream oFile("simulazione.dat");

	for (int i = 0; i < N; i++)		//GRIGLIA ASSE x SPAZIATA DI h
		x[i] = i * h;

	for (int i = 0; i < N; i++){
		M[i] =& m[i * N];
		psi[i] = exp(polar(q * x[i], M_PI/2.)) * exp(-pow(x[i] - x0, 2.) / (2. * sigma * sigma));
	}

	norma (psi, h, x);

	for (int i = 0; i < N; i++){						//RIEMPIO MATRICE M
		for (int j = 0; j < N; j++){
			if (j == i)
				M[i][j] = polar(4. * h * h / dt, M_PI/2) - 2. - 2. * h * h * V(x[i]);
			else if (j == i + 1 || j == i - 1)
				M[i][j] = 1.;
			else
				M[i][j] = 0.;
		}
	}

	F[0] = - psi[1] + 2. * psi[0] + polar(4. * h * h / dt, M_PI/2) * psi[0] + 2. * h * h * V(x[0]) * psi[0];

	for (int j = 0; j < N-1; j++){
	F[j] = - psi[j+1] + 2. * psi[j] - psi[j-1] + polar(4. * h * h / dt, M_PI/2) * psi[j] + 2. * h * h * V(x[j]) * psi[j];	
	}

	F[N-1] = 2. * psi[N-1] - psi[N-2] + polar(4. * h * h / dt, M_PI/2) * psi[N-1] + 2. * h * h * V(x[N-1]) * psi[N-1];

	for(int i = 0; i < tstep; i++){
		risolvi_matrice(N, psi, F, M);

		F[0] = - psi[1] + 2. * psi[0] + polar(4. * h * h / dt, M_PI/2) * psi[0] + 2. * h * h * V(x[0]) * psi[0];

		for (int j = 0; j < N-1; j++){
			F[j] = - psi[j+1] + 2. * psi[j] - psi[j-1] + polar(4. * h * h / dt, M_PI/2) * psi[j] + 2. * h * h * V(x[j]) * psi[j];	
		}

		F[N-1] = 2. * psi[N-1] - psi[N-2] + polar(4. * h * h / dt, M_PI/2) * psi[N-1] + 2. * h * h * V(x[N-1]) * psi[N-1];

		if (i % 100 == 0){
			for (int k = 0; k < N; k++){
				oFile << x[k] << "\t" << abs(psi[k] * conj(psi[k])) << "\n";
			}		

			oFile << "\n";
		}	
	}

	return 0;
	
}












