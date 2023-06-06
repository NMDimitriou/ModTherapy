#include "declarations.hpp"

//Read experimental data
void read_ExpDat(ExpData& expr)
{
        char hostname[HOST_NAME_MAX + 1];
        gethostname(hostname, HOST_NAME_MAX + 1);
        printf("hostname: %s\n", hostname);
        printf("\n");

        char filename[4096];
        const char* name = "inExp";
        sprintf(filename, "../../src/%s_%s.txt", name, hostname);

        printf("Importing Experimental data from shared memory\n");
        printf("from file %s\n", filename);

        int ind[3]; //element 0: key for sigma, element 2: key for data, element 1: key for time interval

        ifstream dat;
        dat.open(filename);
        if (dat.is_open())
        {

        	dat >> ind[0];
        	dat >> ind[1];
        	dat >> ind[2];

        	dat.close();
        }

    	expr.sigma2= (double *)shmat(ind[0],NULL,0);
    	expr.dat   = (double *)shmat(ind[1],NULL,0);
    	expr.intrv = (int    *)shmat(ind[2],NULL,0);
    	if (expr.sigma2== NULL) cerr << "Memory allocation for sigma has failed " << endl;
    	if (expr.dat   == NULL) cerr << "Memory allocation for experimental data has failed " << endl;
    	if (expr.intrv == NULL) cerr << "Memory allocation for time interval data has failed" << endl;
    	printf("Completed.\n");
		for(int i=0; i<10; i++) cout << expr.intrv[i] << '\t' << expr.dat[i] << endl;
    	printf(" \n");
}


void read_SimDat(SimData &sim, int rep, char* study, char* data, char* model)
{
	int i=0;
	double tmp0, tmp1, tmp2, tmp3;
	char filename[4096];
	sprintf(filename,"../output_rep_%s_%s_vs_%s/output_rep_%d.txt",study,data,model,rep);
//	cout << filename << endl;
	ifstream dat;
//	dat.precision(16);
//	cout.precision(16);
	dat.open(filename);
	if (dat.is_open())
	{
		while(!dat.eof())
		{
			dat >> tmp0;
			dat >> tmp1;
			dat >> tmp2;
			dat >> tmp3;
			sim.tp.push_back(tmp0);
			sim.dat.push_back(tmp3);
			//cout << "t=" << sim.tp[i] << '\t' << "dat=" << sim.dat[i] << endl;
			i++;
		}
		dat.close();
	}

}

void free_SimDat( SimData &sim )
{
	vector<double>().swap(sim.tp);
	vector<double>().swap(sim.dat);
}

