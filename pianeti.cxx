#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

int main () {
	double x;
	vector <double> dati;
	vector <double> massa;
	vector <double> Rxv;			//componenti ai tempi i
	vector <double> Ryv;
	vector <double> Rzv;
	vector <double> Vx;
	vector <double> Vy;
	vector <double> Vz;
	vector <double> Vx0v;		//velocità corrette
	vector <double> Vy0v;		
	vector <double> Vz0v;
	double Rx[10] = {0};		//componenti ai tempi i+1 (array)
	double Ry[10] = {0};
	double Rz[10] = {0};
	double Vx0[10] = {0};
	double Vy0[10] = {0};
	double Vz0[10] = {0};
	double vcmx, vcmy, vcmz, M = 0, VMx = 0, VMy = 0, VMz = 0, dt = 86400, G = 6.67e-11;
	int T = 0;

	ifstream iFile ("dati.txt");
	if (!iFile)
		cout<<"errore apertura file" <<endl;

	ofstream oFileSo ("object1.dat");	
	ofstream oFileMe ("object2.dat");	
	ofstream oFileV ("object3.dat");	
	ofstream oFileT ("object4.dat");	
	ofstream oFileMa ("object5.dat");	
	ofstream oFileG ("object6.dat");	
	ofstream oFileSa ("object7.dat");	
	ofstream oFileU ("object8.dat");	
	ofstream oFileN ("object9.dat");	
	ofstream oFileH ("object10.dat");
	ofstream oFileEnergia ("energia.txt");	

	while (iFile >> x)
		dati.push_back(x);

	for (int i = 0; i < dati.size(); i = i + 7){
		massa.push_back(dati[i+1]);
		Rxv.push_back(dati[i+2]*1000);
		Ryv.push_back(dati[i+3]*1000);
		Rzv.push_back(dati[i+4]*1000);
		Vx.push_back(dati[i+5]*1000);
		Vy.push_back(dati[i+6]*1000);
		Vz.push_back(dati[i+7]*1000);
	}
	
	for (int j = 0; j < 10; j++){				//per calcolare velocità rispetto al cm
		M = M + massa[j];
		VMx = VMx + massa[j] * Vx[j]; 
		VMy = VMy + massa[j] * Vy[j]; 
		VMz = VMz + massa[j] * Vz[j]; 
	}

	vcmx = VMx/M;			//velocità del cm
	vcmy = VMy/M;
	vcmz = VMz/M;

	for (int k = 0; k < 10; k++){
		Vx0v.push_back(Vx[k] - vcmx);
		Vy0v.push_back(Vy[k] - vcmy);
		Vz0v.push_back(Vz[k] - vcmz);
	}

	for (int i = 0; i < 10; i++){						//scrivo velocità e posizioni come array
		Rx[i] = Rxv[i]; Ry[i] = Ryv[i]; Rz[i] = Rzv[i];
		Vx0[i] = Vx0v[i]; Vy0[i] = Vy0v[i]; Vz0[i] = Vz0v[i];
	}


do{		
	double Rx1[10] = {0};		//componenti ai tempi i+1
	double Ry1[10] = {0};
	double Rz1[10] = {0};
	double Vx01[10] = {0};
	double Vy01[10] = {0};
	double Vz01[10] = {0};
	double fx[10] = {0};
	double fy[10] = {0};
	double fz[10] = {0};
	double fx1[10] = {0};
	double fy1[10] = {0};
	double fz1[10] = {0};
	double energiap[10] = {0}; double energiac[10] = {0};
	double energ = 0;
	for (int i = 0; i < 10; i++){
		for (int j = 0; j < 10; j++){
			if (j != i){
				fx[i] = fx[i] + (G * massa[i] * massa[j] * (Rx[j] - Rx[i])/(pow(pow(Rx[i] - Rx[j],2) + pow(Ry[i] - Ry[j],2) + pow(Rz[i] - Rz[j],2),3./2.)));
				fy[i] = fy[i] + (G * massa[i] * massa[j] * (Ry[j] - Ry[i])/(pow(pow(Rx[i] - Rx[j],2) + pow(Ry[i] - Ry[j],2) + pow(Rz[i] - Rz[j],2),3./2.)));
				fz[i] = fz[i] + (G * massa[i] * massa[j] * (Rz[j] - Rz[i])/(pow(pow(Rx[i] - Rx[j],2) + pow(Ry[i] - Ry[j],2) + pow(Rz[i] - Rz[j],2),3./2.)));
			}
		}
	}

	for (int i = 0; i < 10; i++){
		Rx1[i] = Rx[i] + Vx0[i] * dt + 1./(2. * massa[i]) * fx[i] * dt * dt;
		Ry1[i] = Ry[i] + Vy0[i] * dt + 1./(2. * massa[i]) * fy[i] * dt * dt;
		Rz1[i] = Rz[i] + Vz0[i] * dt + 1./(2. * massa[i]) * fz[i] * dt * dt;
	}


	for (int i = 0; i < 10; i++){
		for (int j = 0; j < 10; j++){
			if (i != j){
				fx1[i] = fx1[i] + (G * massa[i] * massa[j] * (Rx1[j] - Rx1[i])/pow(pow(Rx1[i] - Rx1[j],2) + pow(Ry1[i] - Ry1[j],2) + pow(Rz1[i] - Rz1[j],2),3./2.));
				fy1[i] = fy1[i] + (G * massa[i] * massa[j] * (Ry1[j] - Ry1[i])/pow(pow(Rx1[i] - Rx1[j],2) + pow(Ry1[i] - Ry1[j],2) + pow(Rz1[i] - Rz1[j],2),3./2.));
				fz1[i] = fz1[i] + (G * massa[i] * massa[j] * (Rz1[j] - Rz1[i])/pow(pow(Rx1[i] - Rx1[j],2) + pow(Ry1[i] - Ry1[j],2) + pow(Rz1[i] - Rz1[j],2),3./2.));
				energiap[i] = energiap[i] - 0.5 * (G * massa[i] * massa[j]/sqrt(pow(Rx[i] - Rx[j],2) + pow(Ry[i] - Ry[j],2) + pow(Rz[i] - Rz[j],2))) ;
			}
		}
	}

	for (int i = 0; i < 10; i++){
		Vx01[i] = Vx0[i] + 1./(2 * massa[i]) * (fx[i] + fx1[i]) * dt;
		Vy01[i] = Vy0[i] + 1./(2 * massa[i]) * (fy[i] + fy1[i]) * dt;
		Vz01[i] = Vz0[i] + 1./(2 * massa[i]) * (fz[i] + fz1[i]) * dt;
		energiac[i] = 0.5 * massa[i] * (pow(Vx0[i],2) + pow(Vy0[i],2) + pow(Vz0[i],2));
		energ = energ + (energiap[i] + energiac[i]);
	}

	if (T % 100 ==0){
	oFileSo << T << "\t" << Rx[0] << "\t" << Ry[0] << "\t" << Rz[0] << "\n";
	oFileMe << T << "\t" << Rx[1] << "\t" << Ry[1] << "\t" << Rz[1] << "\n";
	oFileV << T << "\t" << Rx[2] << "\t" << Ry[2] << "\t" << Rz[2] << "\n";
	oFileT << T << "\t" << Rx[3] << "\t" << Ry[3] << "\t" << Rz[3] << "\n";
	oFileMa << T << "\t" << Rx[4] << "\t" << Ry[4] << "\t" << Rz[4] << "\n";
	oFileG << T << "\t" << Rx[5] << "\t" << Ry[5] << "\t" << Rz[5] << "\n";
	oFileSa << T << "\t" << Rx[6] << "\t" << Ry[6] << "\t" << Rz[6] << "\n";
	oFileU << T << "\t" << Rx[7] << "\t" << Ry[7] << "\t" << Rz[7] << "\n";
	oFileN << T << "\t" << Rx[8] << "\t" << Ry[8] << "\t" << Rz[8] << "\n";
	oFileH << T << "\t" << Rx[9] << "\t" << Ry[9] << "\t" << Rz[9] << "\n";
	oFileEnergia << T << "\t" << energ <<"\n";
	}


	T += 1;

	for (int i = 0; i < 10; i++){
		Rx[i] = Rx1[i]; Ry[i] = Ry1[i]; Rz[i] = Rz1[i];
		Vx0[i] = Vx01[i]; Vy0[i] = Vy01[i]; Vz0[i] = Vz01[i];
	}

	} while (T < 100000);

	return 0;
}