#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <unistd.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <math.h>


using namespace std;

struct SimData
{
	//double *dat;
	//double *tp;
	vector<double> dat;
	vector<double> tp;
};

struct ExpData
{
    double *sigma2;
    double *dat;
    int    *intrv;
};


void read_ExpDat( ExpData &exp );
void read_SimDat(SimData &sim, int rep, char* study, char* data, char* model);
void free_SimDat( SimData &sim );

#endif //declarations_h

