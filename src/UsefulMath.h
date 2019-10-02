#pragma once

#include <math.h>

#include "SafeOps.h"


class UsefulMath {

	UsefulMath() {};
	~UsefulMath() {};

public:

	static double dddotprod ( double * a, double * b ) {
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}


	static double didotprod ( double * a, int * b ) {
		return a[0]*(double)b[0] + a[1]*(double)b[1] + a[2]*(double)b[2];
	}


	static int iidotprod ( int * a, int * b ) {
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}

/*
	static double min ( double a, double b ) {
		if ( a > b ) 
			return b;
		else 
			return a;
	}
*/
	// invert a NxN matrix using routines from numerical recipes
	static void invert_matrix( int n, double **a, double **ai ) {

		int *indx;
		double *col,d;

		SafeOps::malloc( col,  n*sizeof(double), __LINE__, __FILE__ );
		SafeOps::malloc( indx, n*sizeof(int),    __LINE__, __FILE__ );
	
		LU_decomp( a, n, indx, &d );

		for( int i=0; i<n; i++ ){
			for( int j=0; j<n; j++ )
				col[j] = 0.;
			col[i] = 1.;
			LU_bksb( a, n, indx, col );
			for( int j=0; j<n; j++ )
				ai[j][i] = col[j];
		}

		free(col);
		free(indx);

	}

	// numerical recipes routines for inverting a general matrix 
	static void LU_decomp( double **a, int n, int *indx, double *d )
	{
		const double TINY = 1.0e-20;

		int imax=0;
		double big,dum,sum,temp,*vv;

		SafeOps::malloc( vv, n*sizeof(double), __LINE__, __FILE__ );
	
		*d = 1.0;
		for( int i=0; i<n; i++ ) {
			big=0.0;
			for( int j=0; j<n; j++ )
				if ((temp=fabs(a[i][j])) > big) 
					big = temp;
			// note big cannot be zero 
			if(big == 0.0) {
				Output::err( "LU_decomp: matrix to invert cannot be singular.\n" );
				throw attempted_singular_mtx_inv;
			}
			vv[i]=1.0/big;
		}

		for( int j=0; j<n; j++ ) {
			for( int i=0; i<j; i++ ) {
				sum=a[i][j];
				for( int k=0; k<i; k++ )
					sum -= a[i][k]*a[k][j];
				a[i][j]=sum;
			}
			big=0.0;

			for( int i=j; i<n; i++ ) {
				sum=a[i][j];
				for( int k=0; k<j; k++ )
					sum -= a[i][k]*a[k][j];
				a[i][j]=sum;
				if( (dum=vv[i]*fabs(sum)) >= big) {
					big=dum;
					imax=i;
				}
			}
			if( j != imax) {
				for( int k=0; k<n; k++ ) {
					dum = a[imax][k];
					a[imax][k] = a[j][k];
					a[j][k] = dum;
				}
				*d = -(*d);
				vv[imax] = vv[j];
			}
			indx[j]=imax;

			if (a[j][j] == 0.0)
				a[j][j]=TINY;
			if (j != n) {
				dum = 1.0/(a[j][j]);
			for( int i=j+1; i<n; i++)
				a[i][j] *= dum;
			}
		}
		free(vv);
	}


	static void LU_bksb(double **a,int n,int *indx,double b[])
	{
		int ii=-1, ip;
		double sum;

		for( int i=0;i<n;i++) {
			ip    = indx[i ];
			sum   = b   [ip];
			b[ip] = b   [i ];
			if( ii != -1 ) 
				for( int j=ii; j<=i-1; j++ )
					sum -= a[i][j]*b[j];
			else if(sum)
				ii=i;
			b[i]=sum;
		}

		for( int i=n-1; i>=0; i-- ) {
			sum=b[i];
			for( int j=i+1; j<n; j++ )
				sum -= a[i][j]*b[j];
			b[i]=sum/a[i][i];
		}
	}


	static double factorial( int n )
	{
		int i;
		double fac = 1.0;
		for (i=2;i<=n;i++)
			fac *= i;
		return fac;
	}

};