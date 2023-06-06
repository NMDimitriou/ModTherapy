#include "declarations.hpp"

//****************************************************************************80

double *discrete_cdf_to_xy ( int n1, int n2, double cdf[], int n, double u[],
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_CDF_TO_XY finds XY points corresponding to discrete CDF values.
//
//  Discussion:
//
//    This program is given a discrete CDF function and a set of N random
//    values U.  Each value of U corresponds to a particular (I,J) subregion
//    whose CDF value just exceeds the value of U.  Inside that subregion,
//    we pick a point at random - this is equivalent to assuming the PDF
//    is constant over the subregion.
//
//    This function is part of an example program, for which various
//    assumptions have been made.  In particular, the region is the unit
//    square, and the subregions are formed by a 20 by 20 grid of equal
//    subsquares.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the number of rows and columns of data.
//
//    Input, double CDF[N1*N2], the CDF values associated with each 
//    subcell.  A particular ordering has been given to the subcells so that the
//    CDF is a monotonoe function when the subcells are listed in that order.
//
//    Input, int N, the number of sample points.
//
//    Input, double U[N], N random values.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double DISCRETE_CDF_TO_XY[2*N], the sample points.
//
{
  	double high;
  	int i;
  	int j;
  	int k;
  	double low;
  	double *r;
  	double *xy;

  	xy = new double[2*n];

  	low = 0.0;
  	for ( j = 0; j < n2; j++ )
  	{
    	for ( i = 0; i < n1; i++ )
    	{
      		high = cdf[i+j*n1];
      		for ( k = 0; k < n; k++ )
      		{
        		if ( low <= u[k] && u[k] <= high )
        		{
          			r = r8vec_uniform_01_new ( 2, seed );
          			xy[0+k*2] = ( ( double ) ( i ) + r[0] ) / ( double ) n1;
          			xy[1+k*2] = ( ( double ) ( j ) + r[1] ) / ( double ) n2;
          			delete [] r;
        		}
      		}
      		low = high;
    	}
  	}
  	return xy;
}
//****************************************************************************80


double *r8mat_copy_new ( int m, int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
//
{
  	double *a2;
  	int i;
	int j;

  	a2 = new double[m*n];

  	for ( j = 0; j < n; j++ )
  	{
    	for ( i = 0; i < m; i++ )
    	{
      		a2[i+j*m] = a1[i+j*m];
    	}
  	}
  	return a2;
}
//****************************************************************************80


void r8mat_scale ( int m, int n, double s, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SCALE multiplies an R8MAT by a scalar.
//
//  Discussion:
//
//    An R8MAT is an array of R8 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double S, the scale factor.
//
//    Input/output, double A[M*N], the matrix to be scaled.
{
  	int i;
  	int j;

  	for ( j = 0; j < n; j++ )
  	{
    	for ( i = 0; i < m; i++ )
    	{
      		a[i+j*m] = a[i+j*m] * s;
    	}
  	}
  	return;
}
//****************************************************************************80

double r8mat_sum ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SUM returns the sum of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], the array.
//
//    Output, double R8MAT_SUM, the sum of the entries.
//
{
	int i;
  	int j;
  	double value;

  	value = 0.0;
  	for ( j = 0; j < n; j++ )
  	{
    	for ( i = 0; i < m; i++ )
    	{
      		value = value + a[i+j*m];
    	}
  	}
  	return value;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
{
  	int i;
  	int i4_huge = 2147483647;
  	int k;
  	double *r;

  	if ( *seed == 0 )
  	{
    	cerr << "\n";
    	cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    	cerr << "  Input value of SEED = 0.\n";
    	exit ( 1 );
  	}

  	r = new double[n];

  	for ( i = 0; i < n; i++ )
  	{
   		k = *seed / 127773;

    	*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    	if ( *seed < 0 )
    	{
      		*seed = *seed + i4_huge;
    	}

    	r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  	}

  	return r;
}
//****************************************************************************80

double *set_discrete_cdf ( int n1, int n2, double pdf[] )

//****************************************************************************80
//
//  Purpose:
//
//    SET_DISCRETE_CDF sets a CDF from a discrete PDF.
//
//  Discussion:
//
//    Here, we proceed from cell (1,1) to (2,1) to (20,1), (1,2), (2,2)...(20,20).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the number of rows and columns of data.
//
//    Input, double PDF[N1,N2], the discrete PDF for the cell (I,J),
//    normalized so that the sum over all cells is 1.
//
//    Output, double CDF[N1,N2], the discrete CDF for the cell (I,J).
//    The last entry of CDF should be 1.
//
{
  	double *cdf;
  	int i;
  	int j;
	double total;

  	cdf = new double[n1*n2];

  	total = 0.0;
  	for ( j = 0; j < n2; j++ )
  	{
    	for ( i = 0; i < n1; i++ )
    	{
      		total = total + pdf[i+j*n1];
      		cdf[i+j*n1] = total;
    	}	
  	}
  	return cdf;
}
//****************************************************************************80
