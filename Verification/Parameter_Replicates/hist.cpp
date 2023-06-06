#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <random>
#include <cstring>
using namespace std;

int main(int argc, char* argv[]){

	int dim = atoi(argv[1]);
	int num_dat;

   	char src[256];
   	sprintf(src,"../../src/%s/tmcmc.par",argv[3]);
   	FILE *f = fopen(src, "r");
    if (f == NULL) {
        return 1;
    }
	char line[256];
	int line_no = 0;
    while (fgets(line, 256, f)!= NULL) {
        line_no++;
		if ((line[0] == '#')||(strlen(line)==0)) {
            printf("ignoring line %d\n", line_no);
            continue;
        }
		else if (strstr(line, "PopSize")) {
            sscanf(line, "%*s %d", &num_dat);
        }
	}
	fclose(f);

   	int nbins=num_dat/4, i=0, nsamp=5000;
   	//double rs[num_dat],rr[num_dat],K[num_dat],sigma2[num_dat],logl[num_dat],ev[num_dat]; 
   	double par[num_dat][dim], tmp[num_dat];
   	double x_h[nbins+1][dim], p_hist[nbins+1][dim];
   	double rnd[nsamp][dim];
   	double dx, norm, u, sp;
   	int j,k,n;

   	srand(time(NULL)); // initialize random number generator
   	random_device rand_dev;
   	mt19937 generator(rand_dev());

   	char fname[1024];
   	char hname[1024];
   	char uname[1024];
   	sprintf(fname,"final_%s.txt",argv[2]);
   	sprintf(hname,"hist_%s.txt",argv[2]);
   	sprintf(uname,"rnd_%s.txt",argv[2]);

   	//Read sampling file
   	fstream ifile;
   	ifile.open(fname,ios::in);
   	ifile.precision(16);
   	cout.precision(16);
   	if(ifile.is_open())
   	{
   		while(i<num_dat){

      	for(j=0;j<dim;j++) ifile >> par[i][j];
	  	cout << par[i][0] << ", " << par[i][1] << ", " << par[i][2] << endl;
	  	i++;

		}
		ifile.close();
   	}

   	//Make histogram
   	for(j=0;j<dim;j++){
   
     	for(k=0;k<num_dat;k++) tmp[k] = par[k][j];
     	sort(tmp, tmp+num_dat);

  	 	dx = (tmp[num_dat-1]-tmp[0])/nbins;  

	 	for(k=0;k<nbins+1;k++){x_h[k][j] = tmp[0] + k*dx; p_hist[k][j]=0.;} 

	 	for(i=0; i<num_dat; i++){
	   	for(k=0;k<nbins+1;k++)
	     	if(tmp[i]<=x_h[k+1][j] && tmp[i]>=x_h[k][j]) p_hist[k][j] ++; 

	   	//if(tmp[i]<=x_h[nbins][j] && tmp[i]>=x_h[nbins-1][j]) p_hist[nbins-1][j] ++;
	 	}
	 	norm=0.;
	 	for(k=0;k<nbins+1;k++) norm += p_hist[k][j];
	 	for(k=0;k<nbins+1;k++) if(p_hist[k][j]>0) {p_hist[k][j] /= norm;}

	 	//Inverse transform sampling
     	for(n=0;n<nsamp;n++){
       	u=(double)rand()/(double)RAND_MAX;

       	for(k=0;k<nbins+1;k++){
         sp += p_hist[k][j];
         if((sp-p_hist[k][j])<u && u<=sp) rnd[n][j] = x_h[k][j]; //rfile << x_h[k] << endl;
       	}
	   	sp=0;
     	}

   	}


   	//Save Probabilities
   	ofstream ofile;
   	ofile.precision(16);
   	ofile.open(hname);
   	for(k=0;k<nbins;k++){
     for(j=0;j<dim;j++) ofile << x_h[k][j] << '\t' << p_hist[k][j] << '\t';
	 ofile << endl;
   	}

   	//Save replicates
   	ofstream rfile;
   	rfile.precision(16);
   	rfile.open(uname);
   	for(n=0;n<nsamp;n++){
     for(j=0;j<dim;j++) rfile << rnd[n][j] << '\t'; 
     rfile << endl;
   	}
   	rfile.close();

	return 0;
}//hist.cpp
