#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double fgrav (double r, double R, double m1, double m2){
	double G = 6.67e-20;							//[km^3/(kg s^2)]
	double F;
	F = G * (m1 * m2 / pow(R, 3)) * r;

	return F;
}

void calcola_forze (double *rx, double *ry, double *rz, double *vx, double *vy, double *vz, double *masse, int N){
	int dt = 86400;								//dati sulle posizioni sono presi ogni 86400 s (1 giorno)
	double R;									//distanza tra i due corpi
	double G = 6.67e-20;						//costante di gravitazione universale

	double *rx1; double *ry1; double *rz1;		//componenti posizioni al tempo n+1	
	double *vx1; double *vy1; double *vz1;		//componenti velocità al tempo n+1
	double *fx; double *fy; double *fz;			//componenti forza al tempo n	
	double *fx1; double *fy1; double *fz1;		//componenti forza al tempo n+1

	rx1 = (double*) malloc (N * sizeof(double));
	ry1 = (double*) malloc (N * sizeof(double));
	rz1 = (double*) malloc (N * sizeof(double));
	vx1 = (double*) malloc (N * sizeof(double));
	vy1 = (double*) malloc (N * sizeof(double));
	vz1 = (double*) malloc (N * sizeof(double)); 
	fx = (double*) malloc (N * sizeof(double));
	fy = (double*) malloc (N * sizeof(double));
	fz = (double*) malloc (N * sizeof(double));
	fx1 = (double*) malloc (N * sizeof(double));
	fy1 = (double*) malloc (N * sizeof(double));
	fz1 = (double*) malloc (N * sizeof(double));

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (j != i){
				R = sqrt(pow(rx[i] - rx[j],2) + pow(ry[i] - ry[j],2) + pow(rz[i] - rz[j],2));
				fx[i] += fgrav(rx[j]-rx[i], R, masse[i], masse[j]);
				fy[i] += fgrav(ry[j]-ry[i], R, masse[i], masse[j]);
				fz[i] += fgrav(rz[j]-rz[i], R, masse[i], masse[j]);
			}
		}
	}

	for (int i = 0; i < N; i++){
		rx1[i] = rx[i] + vx[i] * dt + 1./(2. * masse[i]) * fx[i] * dt * dt;
		ry1[i] = ry[i] + vy[i] * dt + 1./(2. * masse[i]) * fy[i] * dt * dt;
		rz1[i] = rz[i] + vz[i] * dt + 1./(2. * masse[i]) * fz[i] * dt * dt;
	}

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (i != j){
				R = sqrt(pow(rx1[i] - rx1[j],2) + pow(ry1[i] - ry1[j],2) + pow(rz1[i] - rz1[j],2));
				fx1[i] += fgrav(rx1[j]-rx1[i], R, masse[i], masse[j]);
				fy1[i] += fgrav(ry1[j]-ry1[i], R, masse[i], masse[j]);
				fz1[i] += fgrav(rz1[j]-rz1[i], R, masse[i], masse[j]);
			}
		}
	}

	for (int i = 0; i < N; i++){
		vx1[i] = vx[i] + 1./(2 * masse[i]) * (fx[i] + fx1[i]) * dt;
		vy1[i] = vy[i] + 1./(2 * masse[i]) * (fy[i] + fy1[i]) * dt;
		vz1[i] = vz[i] + 1./(2 * masse[i]) * (fz[i] + fz1[i]) * dt;
	}

	for (int i = 0; i < N; i++){
		rx[i] = rx1[i]; ry[i] = ry1[i]; rz[i] = rz1[i];
		vx[i] = vx1[i]; vy[i] = vy1[i]; vz[i] = vz1[i];
	}

	free(fx); free(fy); free(fz);

	return;
}

double calcola_energia (double *rx, double *ry, double *rz, double *vx, double *vy, double *vz, double *masse, int N){
	double energia = 0;
	double *potenziale; double *cinetica;
	double G = 6.67e-20;
	double R;

	potenziale = (double*) malloc (N * sizeof(double));
	cinetica = (double*) malloc (N * sizeof(double));

	for (int i = 0; i < N; i++)					//inizializzo a 0 le componenti
		potenziale[i] = 0;

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (i != j){
				R = sqrt(pow(rx[i] - rx[j],2) + pow(ry[i] - ry[j], 2) + pow(rz[i] - rz[j], 2));
				potenziale[i] -= 0.5 * (G * masse[i] * masse[j] / R) ;
				cinetica[i] = 0.5 * masse[i] * (pow(vx[i], 2) + pow(vy[i], 2) + pow(vz[i], 2));

				energia += (potenziale[i] + cinetica[i]);
			}
		}
	}	

	return energia;
}

int main(){
	int N;											//numero corpi
	int Ndati;										//numero dati									
	double R;										//distanza tra i due corpi
	int T = 0;										//indice di tempo
	double vcmx, vcmy, vcmz;						//componenti velocità del centro di massa																		
	double *dati;									//dati: #corpi; massa (kg); rx (km), ry, rz; vx (km/s), vy, vz
	double *masse;
	double *rx; double *ry; double *rz;				//componenti posizioni al tempo n
	double *vx; double *vy; double *vz;				//componenti velocità al tempo n

	cout << "Quanti corpi:" << "\n";
	cin >> N;

	Ndati = N * 7 + 1;

	dati = (double*) malloc (Ndati * sizeof(double));
	masse = (double*) malloc (N * sizeof(double));
	rx = (double*) malloc (N * sizeof(double));
	ry = (double*) malloc (N * sizeof(double));
	rz = (double*) malloc (N * sizeof(double));
	vx = (double*) malloc (N * sizeof(double));
	vy = (double*) malloc (N * sizeof(double));
	vz = (double*) malloc (N * sizeof(double));

	ofstream Sole ("sole.dat");
	ofstream Terra ("terra.dat");
	ofstream Giove ("giove.dat");
	ofstream Saturno ("saturno.dat");
	ofstream Voyager ("voyager.dat");
	ofstream Energia ("energia.txt");

	ifstream iFile ("datiesame.txt");

	double x;
	int j = 0;
	while (iFile >> x){								//riempio array con tutti i dati							
		dati[j] = x;
		j++;
	}

	for (int i = 0; i < N; i++){					//divido i dati
		masse[i] = dati[i*7+1];
		rx[i] = dati[i*7+2];
		ry[i] = dati[i*7+3];
		rz[i] = dati[i*7+4];
		vx[i] = dati[i*7+5];
		vy[i] = dati[i*7+6];
		vz[i] = dati[i*7+7];
	}

	double *nasa;									//dati sonda presi dal database della NASA giorno per giorno 
	double *rxnasa;double *rynasa; double *rznasa;	//componenti posizioni sonda NASA							

	nasa = (double*) malloc (7968 * sizeof(double));
	rxnasa = (double*) malloc (1992 * sizeof(double));
	rynasa = (double*) malloc (1992 * sizeof(double));
	rznasa = (double*) malloc (1992 * sizeof(double));

	ifstream iFilenasa ("nasa.txt");

	int s = 0;
	while (iFilenasa >> x){							//riempio array con tutti i dati
		nasa[s] = x;
		s++;
	}

	for(int i = 0; i < 1992; i++){					//divido i dati
		rxnasa[i] = nasa[i*4+1];
		rynasa[i] = nasa[i*4+2];
		rznasa[i] = nasa[i*4+3];
	}

	do{
		calcola_forze(rx, ry, rz, vx, vy, vz, masse, N);

		if (T % 10 == 0){							//stampo dati ogni 10 giorni
		Sole << T+1 << "\t" << rx[0] << "\t" << ry[0] << "\t" << rz[0] << "\n";
		Terra << T+1 << "\t" << rx[3] << "\t" << ry[3] << "\t" << rz[3] << "\n";
		Giove << T+1 << "\t" << rx[5] << "\t" << ry[5] << "\t" << rz[5] << "\n";
		Saturno << T+1 << "\t" << rx[6] << "\t" << ry[6] << "\t" << rz[6] << "\n";
		Voyager << T+1 << "\t" << rx[9] << "\t" << ry[9] << "\t" << rz[9] << "\n";

		Energia << T << "\t" << calcola_energia(rx, ry, rz, vx, vy, vz, masse, N) <<"\n";
		}

		T += 1;

	} while (T < 580);

	ifstream iFile1 ("9apr1979.txt");						//dati prima correzione

	int k = 0;
	while (iFile1 >> x){								
		dati[k] = x;
		k++;
	}

	rx[9] = dati[0]; ry[9] = dati[1]; rz[9] = dati[2];		//cambio posizione sonda
	vx[9] = dati[3]; vy[9] = dati[4]; vz[9] = dati[5];		//cambio velocità sonda

	do{
		calcola_forze(rx, ry, rz, vx, vy, vz, masse, N);
		
		if (T % 10 == 0){
		Sole << T+1 << "\t" << rx[0] << "\t" << ry[0] << "\t" << rz[0] << "\n";
		Terra << T+1 << "\t" << rx[3] << "\t" << ry[3] << "\t" << rz[3] << "\n";
		Giove << T+1 << "\t" << rx[5] << "\t" << ry[5] << "\t" << rz[5] << "\n";
		Saturno << T+1 << "\t" << rx[6] << "\t" << ry[6] << "\t" << rz[6] << "\n";
		Voyager << T+1 << "\t" << rx[9] << "\t" << ry[9] << "\t" << rz[9] << "\n";

		Energia << T << "\t" << calcola_energia(rx, ry, rz, vx, vy, vz, masse, N) <<"\n";
		}

		T += 1;

	} while (T >= 580 && T < 760);

	ifstream iFile2 ("10ott1979.txt");						//dati seconda correzione

	int q = 0;
	while (iFile2 >> x){								
		dati[q] = x;
		q++;
	}

	vx[9] = dati[3]; vy[9] = dati[4]; vz[9] = dati[5];		//cambio velocità sonda

	do{
		calcola_forze(rx, ry, rz, vx, vy, vz, masse, N);

		if (T % 10 == 0){
		Sole << T+1 << "\t" << rx[0] << "\t" << ry[0] << "\t" << rz[0] << "\n";
		Terra << T+1 << "\t" << rx[3] << "\t" << ry[3] << "\t" << rz[3] << "\n";
		Giove << T+1 << "\t" << rx[5] << "\t" << ry[5] << "\t" << rz[5] << "\n";
		Saturno << T+1 << "\t" << rx[6] << "\t" << ry[6] << "\t" << rz[6] << "\n";
		Voyager << T+1 << "\t" << rx[9] << "\t" << ry[9] << "\t" << rz[9] << "\n";
		Energia << T << "\t" << calcola_energia(rx, ry, rz, vx, vy, vz, masse, N) <<"\n";
		}
	
		T += 1;

	} while (T >= 760 && T < 1170);

	ifstream iFile3 ("12nov1980.txt");						//dati terza correzione

	int p = 0;
	while (iFile3 >> x){								
		dati[p] = x;
		p++;
	}

	rx[9] = dati[0]; ry[9] = dati[1]; rz[9] = dati[2];		//cambio posizione sonda
	vx[9] = dati[3]; vy[9] = dati[4]; vz[9] = dati[5];		//cambio velocità sonda

	do{
		calcola_forze(rx, ry, rz, vx, vy, vz, masse, N);
		
		if (T % 10 == 0){
		Sole << T+1 << "\t" << rx[0] << "\t" << ry[0] << "\t" << rz[0] << "\n";
		Terra << T+1 << "\t" << rx[3] << "\t" << ry[3] << "\t" << rz[3] << "\n";
		Giove << T+1 << "\t" << rx[5] << "\t" << ry[5] << "\t" << rz[5] << "\n";
		Saturno << T+1 << "\t" << rx[6] << "\t" << ry[6] << "\t" << rz[6] << "\n";
		Voyager << T+1 << "\t" << rx[9] << "\t" << ry[9] << "\t" << rz[9] << "\n";

		Energia << T << "\t" << calcola_energia(rx, ry, rz, vx, vy, vz, masse, N) <<"\n";
		}

		T += 1;

	} while (T >= 1170 && T < 2000);

	return 0;
}