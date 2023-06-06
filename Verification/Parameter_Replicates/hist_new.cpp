#include "declarations.hpp"

int main(int argc, char* argv[]){

//	double *cdf;
//  	string filename;
  	int num_dat=1024;
	int reps=5000;
  	int dim= atoi(argv[1]);;
//  	double *pdf;
  	int seed;
//  	double *u;
//  	double *xy;
	double pdf_save[num_dat];
	double par[num_dat][dim], save_rnd[reps][dim];
	int i=0,j;

	char fname[1024];
	char uname[1024];
	sprintf(fname,"final_%s.txt",argv[2]);
	sprintf(uname,"rnd_%s.txt",argv[2]);

	//Read sample file
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

	for(int i=0; i<dim;i++)
	{
		double *cdf;
		double *pdf;
		double *u;
    	double *xy;

		for(j=0;j<num_dat;j++)
        {
            pdf_save[j] = par[j][i];
        }
		
		//  Construct a PDF from the data.
		//pdf = r8mat_copy_new( num_dat, 1, pdf_save );
		pdf = pdf_save;

		double total = r8mat_sum( num_dat, 1, pdf );

        cout << "\n";
        cout << "  PDF data sums to " << total << "\n";

        double scale = 1.0 / total;

        r8mat_scale( num_dat, 1, scale, pdf );

		// "Integrate" the data over rows and columns of the region to get the CDF.
		cdf = set_discrete_cdf ( num_dat, 1, pdf );

		//  Choose N CDF values at random.
  		seed = 123456789;

  		u = r8vec_uniform_01_new ( reps, &seed );

		//  Find the cell corresponding to each CDF value,
		//  and choose a random point in that cell.
  		xy = discrete_cdf_to_xy ( num_dat, 1, cdf, reps, u, &seed );
		cout << "am I here" << endl;
		for(j=0;j<reps;j++)
        {
            save_rnd[j][i] = xy[j];
        }
		//
    	//  Free memory.
    	//
    	delete [] cdf;
    	delete [] pdf;
    	delete [] u;
    	delete [] xy;
	}

	//Save random sample
	ofstream rfile;
   	rfile.precision(16);
   	rfile.open(uname);
	for(i=0;i<reps;i++){
		for(j=0;j<dim;j++) rfile << save_rnd[i][j] << '\t';
		rfile << endl;
	}
	rfile.close();

  	return 1;
}
