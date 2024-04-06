#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main () {

	int k = 0, n = 0, m = 0, p = 0;
	double theta0, omega0, Ts, T, y1teta, y1omega, y2teta, y2omega, y3teta, y3omega, y4teta, y4omega, yrkteta, yrkomega;
	double y1tetas, y1omegas, y2tetas, y2omegas, y3tetas, y3omegas, y4tetas, y4omegas, yrktetas, yrkomegas; 
	double dt;
	double l = 1;
	double g = 9.81;
	int rows = 2;
	int cols = 30000;
	ofstream oFile("dati.txt");

	cout << "angolo di partenza: " << "\n";
	cin >> theta0;
	cout << "omega di partenza: " << "\n";
	cin >> omega0;
	cout << "passo temporale: " << "\n";
	cin >> dt;

	/*double **yes = new double *[rows];			//allocazione di un puntatore per ogni riga	
	for (int i = 0; i < rows; i++)
		yes[i] = new double [cols];				//allocazione dell'array che contiene ogni riga

	yes[0][0] = theta0;
	yes[1][0] = omega0;

	do{
		for (int i = 0; i < rows; i++){			//eulero esplicito
			for (int j = 1; j < cols; j++){
				if (i == 0){
					yes[i][j] = yes[i][j-1] + yes[i+1][j-1] * dt;
				}
				if (i == 1){
					yes[i][j] = yes[i][j-1] -g/l * sin(yes[i-1][j-1]) * dt; 
				}
			}
		}
		k++;
	} while (yes[0][k] * yes[0][k+1] > 0 || yes[0][k] * yes[0][k+1] == 0);

	//cout << k * dt * 4 << "\n";

	for (int i = 0; i < rows; i++)				//deallocazione delle singole righe
		delete[] yes[i];

	delete[] yes;									//deallocazione dei puntatori alle righe

	double **ye = new double *[rows];			//allocazione di un puntatore per ogni riga	
	for (int i = 0; i < rows; i++)
		ye[i] = new double [cols];				//allocazione dell'array che contiene ogni riga


	ye[0][0] = theta0;
	ye[1][0] = omega0;

	do{
		for (int i = 0; i < rows; i++){			//eulero esplicito approssimato
			for (int j = 1; j < cols; j++){
				if (i == 0){
					ye[i][j] = ye[i][j-1] + ye[i+1][j-1] * dt;
				}
				if (i == 1){
					ye[i][j] = ye[i][j-1] -g/l * (ye[i-1][j-1]) * dt; 
				}
			}
		}
		n++;
	} while (ye[0][n] * ye[0][n+1] > 0 || ye[0][n] * ye[0][n+1] == 0);

	//cout << n * dt * 4 << "\n";

	for (int i = 0; i < rows; i++)				//deallocazione delle singole righe
		delete[] ye[i];

	delete[] ye;								//deallocazione dei puntatori alle righe*/

	
	//for (theta0 = 0.001; theta0 < 3.14/2.; theta0 = theta0 + 0.01, omega0 = omega0){			//runge kutta


		do {

	y1tetas = theta0;
	y1omegas = omega0;
	m = 0;

do
{ 
		if (m > 0)
		{
			y1tetas = yrktetas;
			y1omegas = yrkomegas;
		}
		y2tetas = y1tetas + (y1omegas * dt/2.);
		y2omegas = y1omegas - (g/l * sin(y1tetas) * dt/2.);

		y3tetas = y1tetas + (y2omegas * dt/2.);
		y3omegas = y1omegas - (g/l * sin(y2tetas) * dt/2.);

		y4tetas = y1tetas + (y3omegas * dt);
		y4omegas = y1omegas - (g/l * sin(y3tetas) * dt);

		yrktetas = y1tetas + (y1omegas + 2 * y2omegas + 2 * y3omegas + y4omegas) * dt/6.;
		yrkomegas = y1omegas + (-g/l * sin(y1tetas) + 2 * (-g/l * sin(y2tetas)) + 2 * (-g/l * sin(y3tetas)) + (-g/l * sin(y4tetas))) * dt/6.;
		m++;

} while (yrktetas * y1tetas >= 0);

	Ts = dt * m * 4;	

y1teta = theta0;
y1omega = omega0;
p = 0;

	do
{
		if (p > 0)
		{
			y1teta = yrkteta;
			y1omega = yrkomega;
		}
		y2teta = y1teta + (y1omega * dt/2.);
		y2omega = y1omega - (g/l * (y1teta) * (dt/2.));

		y3teta = y1teta + (y2omega * dt/2.);
		y3omega = y1omega - (g/l * (y2teta) * (dt/2.));

		y4teta = y1teta + (y3omega * dt);
		y4omega = y1omega - (g/l * (y3teta) * dt);

		yrkteta = y1teta + (y1omega + 2 * y2omega + 2 * y3omega + y4omega) * dt/6.;
		yrkomega = y1omega + (-g/l * (y1teta) + 2 * (-g/l * (y2teta)) + 2 * (-g/l * (y3teta)) + (-g/l * (y4teta))) * dt/6.;
		p++;

} while (yrkteta * y1teta >= 0);

  T = 4 * dt * p;

  oFile << theta0 << "\t" << Ts/T <<"\n";

  theta0 = theta0 + 0.01;

} while (theta0 < 3.14/2.);

	return 0;
}
