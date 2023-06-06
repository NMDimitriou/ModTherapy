#include "declarations.hpp"

ExpData expr;
SimData sim;


int main(int argc, char* argv[])
{
	
	int i,j;
	int reps = atoi(argv[1]);
	int total= atoi(argv[2]);
	int idx; 
	double tpdat[total];
	vector<double> temp; 

	//Load experimental data
	read_ExpDat(expr);

	for(i=0;i<total;i++) tpdat[i] = 0.0;
	
	//P-value for every time-point
	for(i=1;i<=reps;i++) 
	{
		//Load simulation data
		read_SimDat(sim,i,argv[3],argv[4],argv[5]);
		temp = sim.tp;
		//cout << "rep=" << i << endl;

		for(j=0; j<total; j++)
    	{
        	int value = expr.intrv[j];
        	for_each(sim.tp.begin(), sim.tp.end(), [value](auto& x) { x = fabs(x-value); });
        	idx    = min_element(sim.tp.begin(),sim.tp.end()) - sim.tp.begin();
        	sim.tp = temp;
			//cout << "t_exp=" << value << '\t' << "t_sim=" << sim.tp[idx] << '\t' <<"exp=" << expr.dat[j] << '\t' << "sim=" << sim.dat[idx] << endl;

			if(expr.dat[j] < sim.dat[idx]) tpdat[j] ++; 
    	}

		//Unload simulation data
		free_SimDat(sim);
	}

	//Write the output
	ofstream out;
	char outname[4096];
	sprintf(outname,"ppp_%s_%s_%s.txt",argv[3],argv[4],argv[5]);
    out.open(outname);
	for(j=0;j<total;j++)
	{
		out << "t_{" << j << "}" << '\t' << "&" << '\t' << tpdat[j]/reps << endl;
	}


}
