#include <stdio.h>
#include <math.h>

#include <string>
#include <string.h>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <sstream>  
#include <ctime>   
#include <functional> 
//#include "pch.h"
#define PUT_LINE() puts("\n---------------\n")

using namespace std;

extern "C"
{
	__declspec(dllexport)  void DisplayHelloFromDLL()
	{
		printf("Hello from DLL truong dep trai !\n");
		
	}
	__declspec(dllexport)  int Addd(int a, int b)
	 {
		return a + b;
	 }
	__declspec(dllexport)  float subb(int a, int b)
	 {
		 return a - b;
	 }
	const int NumGen = 100;
	const int PopSize = 50;
	const int N = 10; //Number of patients
	const int M = 5; //Number of machines
	const int MaxN = 5; //Maximun therapies that each patient is treated
	const double Pc = 0.8;
	const double Pm = 1 - Pc;
	const int Z = 10000;	//Big number
	const int penalty = 100;	//penalty for over serive time
	const int Lmax = 1000;		// maximum service time
	const double lamda = 0.9;

	/* Operations include therapies and conjunctive therapies */

	int no[] = { 3, 3, 2, 3, 3, 3, 3, 3,3,3 };//Number of operations of each patient; the first value is of patient 1.  matrix of number of operations. index is patient and value is number of operation (after conjuncting)
	int n[N] = { 4,4,3,4,4,4,4,4,4,4 };	//matrix of number of therapies of each patient, n[patient]=?, index +1 is patient name
	int PcPointer;
	int PmPointer;
	signed int x[2 * PopSize][N][MaxN + 1];		// operations sequence
	signed int xd[2 * PopSize][N][MaxN + 1];	//therapies sequence
	int Wtt[PopSize * 2];
	int Cmax[PopSize * 2];
	int k;
	int P_st[PopSize * 2][N];					//starting tme of patient
	int P_ct[PopSize * 2][N];					//completion time of patient
	int start_time[PopSize * 2][N][MaxN];		//starting time of therapy
	int complete_time[PopSize * 2][N][MaxN];	//completion time of therapy
	int machine[PopSize * 2][N][MaxN];	// contain machine name which therapy is assigned
	double e[PopSize * 2];
	double Best_value, Best_chro;
	double t0;
	// fitness values
	int tem, ran;

	void Initiation();
	void Crossover();
	void Mutation();
	void SAneighbor();
	void Selection();
	void Decoding();
	void Initial_temperature();

	__declspec(dllexport) int maiin() {// int argc, char **argv

		srand((unsigned)time(NULL));

		try {
			int gen = 0;
			Initiation();
			for (k = 0; k < PopSize; k++) {
				Decoding();
			}

			Initial_temperature();
			int tk = t0;
			double Best_previous = 1;
			int L = 1;
			while (L > 30) {
				cout << "  Generation  " << gen << endl;
				Crossover();
				Mutation();
				for (k = 0; k < PopSize * 2; k++) {
					Decoding();
				}
				Selection();

				for (k = 0; k < PopSize; k++) {
					int k_origin = PopSize + 1;
					e[k_origin] = e[k];
					for (int i = 0; i < N; i++) {
						for (int j = 0; j <= MaxN; j++) {
							x[k_origin][i][j] = x[k][i][j];

						}
					}
					for (int n = 0; n < N*MaxN; n++) {
						SAneighbor();	//a neighbor chro will be located on PopSize position
						int t = k;
						k = PopSize;	// position of neighbor
						Decoding();
						int de = e[k] - e[t];
						if (exp(-de / tk) > (double)rand() / (double)(RAND_MAX)) {
							e[t] = e[k];
							for (int i = 0; i < N; i++) {
								for (int j = 0; j <= MaxN; j++) {
									x[t][i][j] = x[k][i][j];
								}
							}
						}
						k = t;
					}
					//=============== update local best solution =============
					if (e[k_origin] < e[k]) {
						e[k] = e[k_origin];
						for (int i = 0; i < N; i++) {
							for (int j = 0; j <= MaxN; j++) x[k][i][j] = x[k_origin][i][j];
						}
					}
				}

				double Best_current = e[0];
				Best_chro = 0;
				for (int chro = 1; chro < PopSize; chro++) {
					if (e[chro] < Best_current) {
						Best_current = e[chro];
						Best_chro = chro;
					}
				}


				if (Best_previous == Best_current) L++;
				else L = 1;
				Best_previous = Best_current;

				tk = lamda * tk;
				gen++;
			}

			cout << "Best solution = " << Best_chro << endl;
			k = Best_chro;
			Decoding();
			for (int i = 0; i < N; i++) {
				for (int j = 0; j <= MaxN; j++) {
					cout << xd[k][i][j] << "   ";
				}
				cout << endl;
			}
			for (int i = 0; i < N; i++) {
				cout << "Position " << i + 1 << "   Patient " << xd[k][i][0] << " St " << P_st[k][i] << P_ct[k][i] << endl;
				for (int j = 1; j <= n[xd[k][i][0] - 1]; j++) {
					cout << "T " << xd[k][i][j] << " M " << machine[k][i][j] << "  St " << start_time[k][i][j] << "ct" << complete_time[k][i][j];
				}
				cout << "---------------" << endl;
			}
			cout << "Cmax = " << Cmax[k] << " Wtt" << Wtt[k];

			/*
				for (k=0; k<PopSize; k++){
					cout<<"chro k = "<<k<<endl;
					for (int i=0; i<N; i++){
						for (int j=0; j<=MaxN; j++){
							cout<<x[k][i][j]<<"   ";
						}
						cout<<endl;
					}
					cout<<"---------------"<<endl;
				}
			*/
		}
		catch (string & errmsg) {
			cerr << "Error: " << errmsg << endl;
		}
		catch (...) {
			cerr << "Error: unkown error caught!" << endl;
		}
		//system("PAUSE");
		return 12345;
	}
	/*========================================================================*/
	void Initiation() {
		for (int k = 0; k < PopSize; k++) {
			signed int y[N];
			for (int i = 0; i < N; i++) {
				y[i] = i + 1;			// a temporary vector			
			}
			random_shuffle(y, y + N);
			for (int i = 0; i < N; i++) {	//generate a random patients sequence
				x[k][i][0] = y[i];
			}
			for (int i = 0; i < N; i++) {
				signed int v[MaxN];
				for (int j = 1; j <= no[y[i] - 1]; j++) {
					v[j - 1] = j;				//temporary vector
				}
				random_shuffle(v, v + no[y[i] - 1]);
				for (int j = 1; j <= no[y[i] - 1]; j++) {
					x[k][i][j] = v[j - 1];
				}
			}
		}
	}
	/*========================================================================*/
	void Crossover() {
		PcPointer = PopSize;
		for (int c = 1; c <= PopSize * Pc / 2; c++) {
			int k1 = rand() % PopSize;
			int k2 = rand() % PopSize;
			int P1 = rand() % N;
			int P2 = rand() % N;
			int CutP1 = min(P1, P2);
			int CutP2 = max(P1, P2);
			int Copy;
			/*
			cout<<"k1 = "<<k1<<endl;
			cout<<"k2 = "<<k2<<endl;
			cout<<"CutP1 = "<<CutP1<<endl;
			cout<<"CutP2 = "<<CutP2<<endl;

			/*Child 1*/
			for (int i = CutP1; i < CutP2; i++) {
				for (int j = 0; j <= MaxN; j++) {
					x[PcPointer][i][j] = x[k1][i][j];
				}
			}
			int i = 0;
			for (int ik2 = 0; ik2 < N; ik2++) {
				Copy = 1;
				if (i >= CutP1 && i < CutP2) i = CutP2;
				for (int ik1 = CutP1; ik1 < CutP2; ik1++) {
					if (x[k2][ik2][0] == x[k1][ik1][0]) {
						Copy = 0;
					}
				}
				//cout<<"ik2  "<<ik2<<"Copy "<<Copy<<endl;
				if (Copy == 1) {
					for (int j = 0; j <= MaxN; j++) {
						x[PcPointer][i][j] = x[k2][ik2][j];
						//cout<<x[PcPointer][i][j]<<"  ";
					}
					i++;
				}
			}
			PcPointer++;
			/*Child 2*/
			for (int i = CutP1; i < CutP2; i++) {
				for (int j = 0; j <= MaxN; j++) {
					x[PcPointer][i][j] = x[k2][i][j];
				}
			}
			i = 0;
			for (int ik1 = 0; ik1 < N; ik1++) {
				Copy = 1;
				if (i >= CutP1 && i < CutP2) i = CutP2;
				for (int ik2 = CutP1; ik2 < CutP2; ik2++) {
					if (x[k1][ik1][0] == x[k2][ik2][0]) {
						Copy = 0;
					}
				}
				if (Copy == 1) {
					for (int j = 0; j <= MaxN; j++) {
						x[PcPointer][i][j] = x[k1][ik1][j];
					}
					i++;
				}
			}
			PcPointer++;
		}
	}
	/*=====================================================================*/
	void Mutation() {
		PmPointer = PcPointer;
		for (int m = PmPointer; m < PopSize * 2; m++) {
			int k1 = rand() % PopSize;
			int m1 = rand() % N;
			int m2 = rand() % N;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j <= MaxN; j++) {
					x[m][i][j] = x[k1][i][j];
				}
			}
			for (int j = 0; j <= MaxN; j++) {
				x[m][m1][j] = x[k1][m2][j];
				x[m][m2][j] = x[k1][m1][j];
			}

		}
	}
	/*=====================================================================*/
	void SAneighbor() {	//generate a neighbor solution from solution k
		int c = rand() % N;
		int m1 = rand() % no[x[k][c][0] - 1] + 1;
		int m2 = rand() % no[x[k][c][0] - 1] + 1;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j <= MaxN; j++) x[PopSize][i][j] = x[k][i][j];	// new chrome will be located on PopSize
		}
		x[PopSize][c][m1] = x[k][c][m2];
		x[PopSize][c][m2] = x[k][c][m1];
	}
	/*=====================================================================*/
	void Decoding() {	//decoding for chromosome k

		const int Wmax = 20;	// total maximum waiting time of each patient
		const int wmax = 5;	//maximum waiting time of precedence therapies
		const int v1 = 10;	// weighting factor for waiting time term
		const int v2 = 1;		// weighting factor for Cmax term


		int T[N][3][3] = { 1,1,0, 1,2,0, 2,3,4,  1,1,0, 2,2,3, 1,4,0,    2,1,2,  1,3,0,  0,0,0,  1,1,0, 2,2,3, 1,4,0,  1,1,0, 1,2,0, 2,3,4,
						  1,1,0, 2,2,3, 1,4,0,  2,1,2, 1,3,0, 1,4, 0,   1,1,0, 2,2,3, 1,4,0,   2,1,2, 1,3,0, 1,4,0,  1,1,0,  1,2,0, 2,3,4 };	// information matrix of operations and its therapies. Patient order, operation order, column 0 is number of therapies of operation, then name of therapies.


		int TP[N] = { 20,20,20, 17, 15, 17, 16, 21, 21, 29 };	//total processing time of each patient. column index is patient
		int RL[N] = { 15, 10, 5, 2, 6,  8,  6, 9, 3, 10 };	//total relaxing time of each patient. first column  is patient 1
		int P[N][4] = { 5,5,5,5,  5,5,5,5,  5,10,5,0,  3,5,3,6,  3,3,4,5, 3,4,5,5,  4,4,4,4,  5,5,5,6 ,  4,3,6,8, 8,8,5,8 };	// matrix of processing time. Each row index is patient name. row 0 is patient 1
		int rl[N][4] = { 5,5,0,5,  5,0,0,5,  0,0,5,0, 0,0,0,2,  3,3,0,0, 3,0,0,5, 0,2,2,2, 3,0,3,3,  0,0,0,3, 3,3,0,4 };  //matrix of relaxing time.  row i is patient i+1
		int PC[N][4] = { 0,0,0,1,   0,0,1,0,  0,1,0,0, 0,0,1,0,  0,0,0,1, 0,0,1,0,  0,1,0,0,  0,0,1,0,  0,1,0,0,  0,0,0,1 };	// matrix of precedence. first row is patient 1. column is therapy. value = 1 is has precedence therapy
		int Ma[N][4][3] = { 2,1,2,  1,4,0,  1,3,0, 1,5,0,    1,4,0, 2,1,2, 1,3,0, 1,5,0,   2,1,2, 1,3,0, 1,5,0, 0,0,0,  2,1,2, 1,3,0, 1,5,0, 1,4,0,   1,3,0,  1,2,0, 2,4,5, 1,1,0,
						   1,1,0,  1,3,0,  2,2,4,  1,5,0,    1,4,0, 2,1,2, 1,5,0, 1,3,0,  2,1,2,  1,3,0,  1,5,0,  1,4,0, 1,1,0,  1,2,0, 2,3,5, 1,4,0,   2,1,3, 1,2,0, 1,4,0, 1,5,0 };	// a 3D matrix of a set of machines that Oij can be assigned. patient, therapy, machine (column 0 is number of machine, then name of machines)
		Wtt[k] = 0;
		//substitute operations by its therapies
		for (int i = 0; i < N; i++) {
			int Pname = x[k][i][0] - 1;		// patient name
			int j = 1;
			int x1[N][MaxN + 1];
			for (int jo = 1; jo <= no[Pname]; jo++) {
				int Oname = x[k][i][jo] - 1;	// operation name
				for (int t = 1; t <= T[Pname][Oname][0]; t++) {	//all therapies of current operation
					xd[k][i][j] = T[Pname][Oname][t];	//substitute operation by its therapies
					j++;
				}
			}
			xd[k][i][0] = x[k][i][0];
		}
		//Find first timeslot of all machine m belongs to Mij
		int Mw[M + 1][100];
		for (int m = 1; m <= M; m++) {		// set initial time window for all machines
			Mw[m][0] = 1;
			Mw[m][1] = 0;
			Mw[m][2] = Z;
		}
		int dos;		//minimum starting time dos
		int ts[M];		//temporary variable for starting time
		int h;
		int Ps[MaxN];		// time window position which j be assigned
		for (int i = 0; i < N; i++) {

			int Pname = xd[k][i][0] - 1;
			int Wi = Z;
			if (i == 0) dos = 0;
			else dos = P_st[k][i - 1];
			while (Wi > Wmax) {
				int j = 1;
				int Ps_m[M + 1];
				while (j <= n[Pname]) {
					int Tname = xd[k][i][j] - 1;					//therapy name
					int Wm[M + 1][2];							//a timeslot on machine m for Oij

					for (int m = 1; m <= Ma[Pname][Tname][0]; m++) {
						//cout<<"m="<<m;
						int Mname = Ma[Pname][Tname][m];
						h = 1;
						int w = 1;
						//cout<<"ts = "<<ts[m];
						while (w <= Mw[Mname][0]) {
							if (Mw[Mname][2 * h] >= dos + P[Pname][Tname] && Mw[Mname][2 * h - 1] <= Mw[Mname][2 * h] - P[Pname][Tname]) {
								Wm[Mname][0] = Mw[Mname][2 * h - 1];
								Wm[Mname][1] = Mw[Mname][2 * h];
								w = Mw[Mname][0] + 1;				//stop finding timeslot
							}
							h++;								//time window h be assigned
							w++;
						}
						if (Wm[Mname][0] >= dos) ts[m] = Wm[Mname][0];
						else ts[m] = dos;
						Ps_m[Mname] = h - 1; 						//cout<<"ts = "<<ts[m];
					}
					int st = ts[1];
					int ma = Ma[Pname][Tname][1];				// machine be assigned
					for (int m = 1; m <= Ma[Pname][Tname][0]; m++) {	//find start time and machine for therapy jth
						if (ts[m] < st) {
							st = ts[m];
							ma = Ma[Pname][Tname][m];
						}
						else if (ts[m] == st && Wm[ma][1] > Wm[Ma[Pname][Tname][m]][1]) ma = Ma[Pname][Tname][m];
					}
					start_time[k][i][j] = st; 							//cout<<" start_time[k][i][j] = "<<start_time[k][i][j];
					complete_time[k][i][j] = st + P[Pname][Tname]; 			//cout<<"  complete_time[k][i][j]  = "<<complete_time[k][i][j];
					machine[k][i][j] = ma; 								//cout<<"  machine[k][i][j]  = "<<machine[k][i][j];
					dos = complete_time[k][i][j] + rl[Pname][Tname]; 	//cout<<"pname"<<Pname+1<<"tname"<<Tname+1; cout<<" dos = "<<dos;	//update dos
					Ps[j] = Ps_m[ma];											//cout<<" position "<<Ps[j]<<endl;		// time window position which j be assigned
					// Check precedence relation

					if (PC[Pname][xd[k][i][j] - 1] == 1) {

						int dw = start_time[k][i][j] - complete_time[k][i][j - 1];
						if (dw > wmax) {

							int dt = dw - wmax;
							j--;
							dos = start_time[k][i][j] + dt;
							j--;	// to compensate for the increase at out of this if statement							 
						}
					}

					j++;	//next therapy
				}	//finish all therapies of patient ith
				P_st[k][i] = start_time[k][i][1];														//cout<<"  P_st[k][i] "<<P_st[k][i];// starting time of patient ith
				P_ct[k][i] = complete_time[k][i][n[Pname]];											//cout<<"  P_ct[k][i] "<<P_ct[k][i];
				Wi = P_ct[k][i] - P_st[k][i] - TP[Pname] - RL[Pname] + rl[Pname][xd[k][i][n[Pname]] - 1]; 		//cout<<" wi "<<Wi<<endl;
				if (Wi > Wmax) dos = P_st[k][i] + Wi - Wmax;	//update dos		
			}
			// update Mw
			for (int j = 1; j <= n[Pname]; j++) {
				int ma = machine[k][i][j];
				int p = Ps[j];	// position of timeslot on machine ma
				if (Mw[ma][2 * p - 1] == start_time[k][i][j] && Mw[ma][2 * p] == complete_time[k][i][j]) { //absolutely fit time window
					Mw[ma][0]--;
					for (int w = p; w <= Mw[ma][0]; w++) {	// left shift 1 position of time window
						Mw[ma][2 * w - 1] = Mw[ma][2 * w + 1];
						Mw[ma][2 * w] = Mw[ma][2 * w + 2];
					}
				}
				else if (Mw[ma][2 * p - 1]<start_time[k][i][j] && Mw[ma][2 * p]>complete_time[k][i][j]) {//absolutely belongs to time window
					Mw[ma][0]++;
					for (int w = Mw[ma][0]; w >= p + 1; w--) {
						Mw[ma][2 * w] = Mw[ma][2 * w - 2];
						Mw[ma][2 * w - 1] = Mw[ma][2 * w - 3];
					}
					Mw[ma][2 * p] = start_time[k][i][j];	//a new upper bound of timw window p
					Mw[ma][2 * p + 1] = complete_time[k][i][j];	// a new lower bound of time window p+1				
				}
				else if (Mw[ma][2 * p - 1] == start_time[k][i][j] && Mw[ma][2 * p] > complete_time[k][i][j]) {
					Mw[ma][2 * p - 1] = complete_time[k][i][j];	// update new lower bound for time window p				
				}
				else Mw[ma][2 * p] = start_time[k][i][j];	// update new upper bound for time window p
				cout << "Mw:" << ma << "N:" << Mw[ma][0] << "p" << p;
				for (int m = 1; m <= Mw[ma][0]; m++) {
					cout << "TL:" << m << ": " << Mw[ma][2 * m - 1] << "-" << Mw[ma][2 * m];
				}
				cout << endl;
			}
			Wtt[k] += Wi;

			Cmax[k] = P_ct[k][0];
			for (int i = 1; i < N; i++) {
				if (Cmax[k] < P_ct[k][i]) Cmax[k] = P_ct[k][i];
			}
		}
		for (int chro = 0; chro < 2 * PopSize; chro++) {
			int p = 0;
			int n = 0;		//number of patient cannot complete before Lmax
			for (int i = 0; i < N; i++) {
				if (P_ct[chro][i] > Lmax) n++;
			}
			if (Cmax[chro] > Lmax) p = (Cmax[chro] - Lmax)*n*penalty;
			e[chro] = v1 * Wtt[chro] + v2 * Cmax[chro] + p;
		}
		/*// print
		for (int i=0;i<N;i++){
			cout<<"P "<<x[1][i][0]<<" P_st  "<<P_st[k][i]<<"  P_ct  "<<P_ct[k][i]<<endl;
		}
			*/
			//cout<<"Cmax = "<<Cmax[k]<<"  Wtt = "<<Wtt[k]<<endl;
		cout << endl;

	}
	void Selection() {
		int rank[PopSize * 2];	// a vector contains the name of chromosome in decreasing order
		int x_temp[PopSize][N][MaxN + 1];
		int e_temp[PopSize * 2];
		int temp[PopSize * 2];

		/*------------Sorting-----------*/
		for (int k = 0; k < PopSize * 2; k++) temp[k] = e[k];
		for (int i = 0; i < N; i++) rank[i] = i;
		for (int i = 0; i < PopSize * 2 - 1; i++) {
			for (int j = i + 1; j < PopSize * 2; j++) {
				

				if (temp[i] < temp[j]) {
					temp[i] = tem;
					temp[i] = temp[j];
					temp[j] = tem;
					rank[i] = ran;
					rank[i] = rank[j];
					rank[j] = ran;
				}

			}
		}
		/*-------------------------------*/
		for (int k = 0; k < PopSize*0.2; k++) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j <= MaxN; j++) {
					x_temp[k][i][j] = x[rank[k]][i][j];
				}
			}
		}
		// binary tournament without replacement

		int r[PopSize * 2];
		for (int k = ceil(PopSize*0.2); k < 2 * PopSize; k++) {
			r[k] = k;
		}
		int c = ceil(PopSize*0.2);
		random_shuffle(&r[c], &r[PopSize * 2]);		//random order of remaing positions
		int position = ceil(PopSize*0.2);
		for (int k = 0; k < PopSize - position; k++) {
			if (e[rank[r[position + 2 * k]]] < e[rank[r[position + 2 * k + 1]]]) {
				for (int i = 0; i < N; i++) {
					for (int j = 0; j <= MaxN; j++) {
						x_temp[position + k][i][j] = x[rank[r[position + 2 * k]]][i][j];
						e_temp[position + k] = e[rank[r[position + 2 * k]]];
					}
				}
			}
			else
				for (int i = 0; i < N; i++) {
					for (int j = 0; j <= MaxN; j++) {
						x_temp[position + k][i][j] = x[rank[r[position + 2 * k + 1]]][i][j];
						e_temp[position + k] = e[rank[r[position + 2 * k + 1]]];
					}
				}
		}
		for (int k = 0; k < PopSize; k++) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j <= MaxN; j++) {
					x[k][i][j] = x_temp[k][i][j];
					e[k] = e_temp[k];
				}
			}
		}
	}
	void Global_Best_Solution() {
		Best_value = e[0];
		Best_chro = 0;
		for (int k1 = 1; k1 < PopSize; k1++) {
			if (e[k1] < Best_value) {
				Best_value = e[k1];
				Best_chro = k1;
			}
		}
	}
	void Initial_temperature() {
		const double pr = 0.1;
		int e_max = e[0];
		int e_min = e[0];
		for (int chro = 1; chro < PopSize; chro++) {
			if (e[chro] > e_max) e_max = e[chro];
			if (e[chro] < e_min) e_min = e[chro];
		}
		t0 = -(e_max - e_min) / log(pr);
	}





















}