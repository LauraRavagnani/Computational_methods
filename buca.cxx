#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

int main() {
	int L = 10; //Bohr
	int N = 1000; //numero di grid step lungo x
	int m = 1, h = 1, e = 1;
	int E_min = 0; //Hartree
	int E_max = 10; //Hartree
	double DE = 0.001; //Hartree
	double s = 1e-8; //Hartree
	int a = 3; //Bohr
	int n;
	int b = 7; //Bohr
	int V_0 = 1; //Hartree
	int k = 1; //Hartree/Bohr^2
	double dx = 0.001; 
	double x;
	double yb = 0, y0 = 0, y1 = 0, y2 = 0, Ez = 0, yz = 0, ybz = 0, En = 0;
	vector <double> Yb;						//vettore con tutti i valori della funzione corrispondenti a E
	vector <double> Energie;				//vettore con tutti i valori delle energie da 0 a 10
	vector <double> avl;					//vettore con autovalori E
	vector <double> ast;					//vettore con autostati corrispondenti a avl
	vector <double> Esol;					//vettore con soluzioni analitiche En
	
	ofstream autovalori ("autovalori.txt");
	ofstream bucapiatta ("bucapiatta.txt");
	ofstream scalino ("scalino.txt");
	ofstream armonico ("armonico.txt");

	cout << "Seleziona il potenziale desiderato: " << "\n" << "1 - Buca piatta" << "\n" << "2 - Potenziale a scalino" << "\n" << "3 - Potenziale armonico" << "\n";
	cin >> n;

if (n == 1){
	for (double E = E_min; E < E_max; E += DE){
		x = 0;
		y0 = 0;
		y1 = 1;
		double M = m * dx * dx/(6*h*h);
		do {
			y2 = (1/(1 + M * E)) * ((- M * 10 * E + 2) * y1 + (- M * E - 1) * y0);

			y0 = y1; 
			y1 = y2;

			x += dx;

		} while (x <= L);

		yb = y1;

		autovalori << E << "\t" << yb <<"\n";
		Yb.push_back(yb);
		Energie.push_back(E);
	}

	for (int i = 0; i < Yb.size(); i++){
		if (Yb[i] * Yb[i+1] < 0){
			do {
				Ez = (Energie[i] + Energie[i+1])/2.;

				x = 0;
				y0 = 0;
				y1 = 1;
				double M = m * dx * dx/(6*h*h);
				do {
					yz = (1/(1 + M * Ez)) * ((- M * 10 * Ez + 2) * y1 + (- M * Ez - 1) * y0);

					y0 = y1; 
					y1 = yz;

					x += dx;

				} while (x <= L);

				ybz = y1;

				if (Yb[i] * ybz < 0){
					Yb[i+1] = ybz;
					Energie[i+1] = Ez;
				}

				else if (ybz * Yb[i+1] < 0){
					Yb[i] = ybz;
					Energie[i] = Ez;
				}

			} while (abs(Energie[i] - Energie[i+1]) > s);

			avl.push_back(Ez);
			ast.push_back(ybz);
		}
	}
 
 	for (int i = 0; i < 14; i++){							
 		En = h * h/(2. * m) * pow(M_PI/L * (i+1), 2);
 		Esol.push_back(En);
 		cout << avl[i] <<"\t" << Esol[i] <<"\n";				//verifica valore degli avl con soluzione analitica

 		x = 0;
		y0 = 0;
		y1 = 1;
		double M = m * dx * dx/(6*h*h);
		do {
			y2 = (1/(1 + M * avl[i])) * ((- M * 10 * avl[i] + 2) * y1 + (- M * avl[i] - 1) * y0);

			y0 = y1; 
			y1 = y2;

			bucapiatta << x << "\t" << y1 <<"\n";

			x += dx;

		} while (x <= L);

 	}
}

else if (n == 2){
	for (double E = E_min; E < E_max; E += DE){
		x = 0;
		y0 = 0;
		y1 = 1;
		double M = m * dx * dx/(6*h*h);
		do {
			y2 = (1/(1 + M * E)) * ((- M * 10 * E + 2) * y1 + (- M * E - 1) * y0);

			y0 = y1; 
			y1 = y2;

			x += dx;

		} while (x <= a);

		do {
			y2 = (1/(1 + M * (E-V_0))) * ((- M * 10 * (E-V_0) + 2) * y1 + (- M * (E-V_0) - 1) * y0);

			y0 = y1; 
			y1 = y2;

			x += dx;

		} while (x <= b);

		do {
			y2 = (1/(1 + M * E)) * ((- M * 10 * E + 2) * y1 + (- M * E - 1) * y0);

			y0 = y1; 
			y1 = y2;

			x += dx;

		} while (x <= L);

		yb = y1;

		autovalori << E << "\t" << yb <<"\n";
		Yb.push_back(yb);
		Energie.push_back(E);
	
	}

	for (int i = 0; i < Yb.size(); i++){
		if (Yb[i] * Yb[i+1] < 0){
			do {
				Ez = (Energie[i] + Energie[i+1])/2.;

				x = 0;
				y0 = 0;
				y1 = 1;
				double M = m * dx * dx/(6*h*h);
				if (x <= a) {
					yz = (1/(1 + M * Ez)) * ((- M * 10 * Ez + 2) * y1 + (- M * Ez - 1) * y0);

					y0 = y1; 
					y1 = yz;

					x += dx;
				}

				else if (x > a && x <= b) {
					yz = (1/(1 + M * (Ez-V_0))) * ((- M * 10 * (Ez-V_0) + 2) * y1 + (- M * (Ez-V_0) - 1) * y0);

					y0 = y1; 
					y1 = yz;

					x += dx;
				}

				else if (x > b && x <= L) {
					yz = (1/(1 + M * Ez)) * ((- M * 10 * Ez + 2) * y1 + (- M * Ez - 1) * y0);

					y0 = y1; 
					y1 = yz;

					x += dx;
				}

				ybz = y1;

				if (Yb[i] * ybz < 0){
					Yb[i+1] = ybz;
					Energie[i+1] = Ez;
				}

				else if (ybz * Yb[i+1] < 0){
					Yb[i] = ybz;
					Energie[i] = Ez;
				}

			} while (abs(Energie[i] - Energie[i+1]) > s);

			avl.push_back(Ez);
			ast.push_back(ybz);
		}
	}
 
 	for (int i = 0; i < 13; i++){							
 		x = 0;
		y0 = 0;
		y1 = 1;
		double M = m * dx * dx/(6*h*h);
		do {
			y2 = (1/(1 + M * avl[i])) * ((- M * 10 * avl[i] + 2) * y1 + (- M * avl[i] - 1) * y0);

			y0 = y1; 
			y1 = y2;

			scalino << x << "\t" << y1 <<"\n";

			x += dx;

		} while (x <= a);

		do {
			y2 = (1/(1 + M * (avl[i]-V_0))) * ((- M * 10 * (avl[i]-V_0) + 2) * y1 + (- M * (avl[i]-V_0) - 1) * y0);

			y0 = y1; 
			y1 = y2;

			scalino << x << "\t" << y1 <<"\n";

			x += dx;

		} while (x <= b);

		do {
			y2 = (1/(1 + M * avl[i])) * ((- M * 10 * avl[i] + 2) * y1 + (- M * avl[i] - 1) * y0);

			y0 = y1; 
			y1 = y2;

			scalino << x << "\t" << y1 <<"\n";

			x += dx;

		} while (x <= L);

 	}
}

if (n == 3){
	for (double E = E_min; E < E_max; E += DE){
		x = 0;
		y0 = 0;
		y1 = 1;
		double M = m * dx * dx/(6*h*h);
		do {
			y2 = (1/(1 + M * (E-k*pow((x-L/2.), 2)))) * ((- M * 10 * (E-k*pow((x-L/2.), 2)) + 2) * y1 + (- M * (E-k*pow((x-L/2.), 2)) - 1) * y0);

			y0 = y1; 
			y1 = y2;

			x += dx;

		} while (x <= L);

		yb = y1;

		autovalori << E << "\t" << yb <<"\n";
		Yb.push_back(yb);
		Energie.push_back(E);
	}

	for (int i = 0; i < Yb.size(); i++){
		if (Yb[i] * Yb[i+1] < 0){
			do {
				Ez = (Energie[i] + Energie[i+1])/2.;

				x = 0;
				y0 = 0;
				y1 = 1;
				double M = m * dx * dx/(6*h*h);
				do {
					yz = (1/(1 + M * (Ez-k*pow((x-L/2.), 2)))) * ((- M * 10 * (Ez-k*pow((x-L/2.), 2)) + 2) * y1 + (- M * (Ez-k*pow((x-L/2.), 2)) - 1) * y0);

					y0 = y1; 
					y1 = yz;

					x += dx;

				} while (x <= L);

				ybz = y1;

				if (Yb[i] * ybz < 0){
					Yb[i+1] = ybz;
					Energie[i+1] = Ez;
				}

				else if (ybz * Yb[i+1] < 0){
					Yb[i] = ybz;
					Energie[i] = Ez;
				}

			} while (abs(Energie[i] - Energie[i+1]) > s);

			avl.push_back(Ez);
			ast.push_back(ybz);
		}
	}
 
 	for (int i = 0; i < 5; i++){							
 		x = 0;
		y0 = 0;
		y1 = 1;
		double M = m * dx * dx/(6*h*h);
		do {
			y2 = (1/(1 + M * (avl[i]-k*pow((x-L/2.), 2)))) * ((- M * 10 * (avl[i]-k*pow((x-L/2.), 2)) + 2) * y1 + (- M * (avl[i]-k*pow((x-L/2.), 2)) - 1) * y0);

			y0 = y1; 
			y1 = y2;

			armonico << x << "\t" << y1 <<"\n";

			x += dx;

		} while (x <= L);

 	}
}



	return 0;
}