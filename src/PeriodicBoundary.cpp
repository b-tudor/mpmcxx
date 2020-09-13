#include "PeriodicBoundary.h"

#include <math.h>
#include <stdio.h>

#include "constants.h"
#include "UsefulMath.h"
#include "Output.h"


PeriodicBoundary::~PeriodicBoundary() {}
PeriodicBoundary:: PeriodicBoundary() {}


PeriodicBoundary & PeriodicBoundary::operator=( const PeriodicBoundary& rhs ) {

	cutoff = rhs.cutoff;
	volume = rhs.volume;

	for(int i=0; i<3; i++)
		for (int j = 0; j < 3; j++) {
			           basis[i][j] = rhs.basis[i][j];
			reciprocal_basis[i][j] = rhs.reciprocal_basis[i][j];
		}

	return *this;
}


// compute volume, cutoff and reciprocal
void PeriodicBoundary::update() {
	compute_volume();
	compute_cutoff();
	compute_reciprocal();
}



// calculates the min cutoff radius from the basis lattice (AKA shortest vector problem)
double PeriodicBoundary::compute_cutoff() {

	double    curr_mag;              // current magnitude
	double    short_mag = MAXVALUE;  // shortest magnitude
	const int MAX_VECT_COEF = 15;
	double    curr_vec[3];

	// Maximum coeff for each basis vector when searching for shortest vector
	if ( volume <= 0 ) 
		return MAXVALUE;

	for ( int i=-MAX_VECT_COEF; i<=MAX_VECT_COEF; i++ ) {
		for ( int j=-MAX_VECT_COEF; j<=MAX_VECT_COEF; j++ ) {
			for ( int k=-MAX_VECT_COEF; k<=MAX_VECT_COEF; k++ ) {
				if ( i == 0 && j == 0 && k == 0 ) 
					continue;
				for ( int p = 0; p < 3; p++ ) 
					curr_vec[p] = i*basis[0][p] + j*basis[1][p] + k*basis[2][p];
				curr_mag = sqrt( UsefulMath::dddotprod(curr_vec, curr_vec));
				if ( curr_mag < short_mag ) 
					short_mag = curr_mag;
			}
		}
	}
	cutoff = 0.5*short_mag;
	return cutoff;
}



// take the determinant of the basis matrix 
double PeriodicBoundary::compute_volume() {

	volume =  basis[0][0]*(basis[1][1]*basis[2][2] - basis[1][2]*basis[2][1]);
	volume += basis[0][1]*(basis[1][2]*basis[2][0] - basis[1][0]*basis[2][2]);
	volume += basis[0][2]*(basis[1][0]*basis[2][1] - basis[1][1]*basis[2][0]);

	return volume;
}



// get the reciprocal space basis 
void PeriodicBoundary::compute_reciprocal() {

	double inverse_volume;

	inverse_volume = 1.0/compute_volume();

	reciprocal_basis[0][0] = inverse_volume*(basis[1][1]*basis[2][2] - basis[1][2]*basis[2][1]);
	reciprocal_basis[0][1] = inverse_volume*(basis[0][2]*basis[2][1] - basis[0][1]*basis[2][2]);
	reciprocal_basis[0][2] = inverse_volume*(basis[0][1]*basis[1][2] - basis[0][2]*basis[1][1]);

	reciprocal_basis[1][0] = inverse_volume*(basis[1][2]*basis[2][0] - basis[1][0]*basis[2][2]);
	reciprocal_basis[1][1] = inverse_volume*(basis[0][0]*basis[2][2] - basis[0][2]*basis[2][0]);
	reciprocal_basis[1][2] = inverse_volume*(basis[0][2]*basis[1][0] - basis[0][0]*basis[1][2]);

	reciprocal_basis[2][0] = inverse_volume*(basis[1][0]*basis[2][1] - basis[1][1]*basis[2][0]);
	reciprocal_basis[2][1] = inverse_volume*(basis[0][1]*basis[2][0] - basis[0][0]*basis[2][1]);
	reciprocal_basis[2][2] = inverse_volume*(basis[0][0]*basis[1][1] - basis[0][1]*basis[1][0]);

}



// output data about the simulation box to stdout
void PeriodicBoundary::printboxdim() {

	char buffer[maxLine];

	sprintf(buffer,"PBC: unit cell volume = %.3f A^3 (cutoff = %.3f A)\n", volume, cutoff);
	Output::out1(buffer);

	sprintf(buffer,"PBC: Basis vectors { a = %8.3lf \tb = %8.3lf \tc = %8.3lf }\n", 
		sqrt(UsefulMath::dddotprod(basis[0], basis[0])),
		sqrt(UsefulMath::dddotprod(basis[1], basis[1])),
		sqrt(UsefulMath::dddotprod(basis[2], basis[2])));
	Output::out1(buffer);

	double alpha = 180.0/pi*acos( UsefulMath::dddotprod(basis[1],basis[2]) / sqrt( UsefulMath::dddotprod(basis[1], basis[1]) * UsefulMath::dddotprod(basis[2], basis[2]) ) );
	double beta  = 180.0/pi*acos( UsefulMath::dddotprod(basis[2],basis[0]) / sqrt( UsefulMath::dddotprod(basis[2], basis[2]) * UsefulMath::dddotprod(basis[0], basis[0]) ) );
	double gamma = 180.0/pi*acos( UsefulMath::dddotprod(basis[0],basis[1]) / sqrt( UsefulMath::dddotprod(basis[0], basis[0]) * UsefulMath::dddotprod(basis[1], basis[1]) ) );

	#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
		sprintf(buffer,"PBC: Basis angles  { alpha = %8.3lf \tbeta  = %8.3lf \tgamma = %8.3lf }\n", alpha, beta, gamma );
	#else
		sprintf(buffer,"PBC: Basis angles  { \u03B1 = %8.3lf \t\u03B2  = %8.3lf \t\u03B3 = %8.3lf }\n", alpha, beta, gamma );
	#endif
	Output::out1(buffer);
}
