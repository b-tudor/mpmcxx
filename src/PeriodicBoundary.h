#pragma once
class PeriodicBoundary
{
public:
	PeriodicBoundary();
	~PeriodicBoundary();

	double cutoff = 0;   // radial cutoff (A)
	double volume = 0;   // unit cell volume (A^3) 
	double            basis[3][3] = {0,0,0,0,0,0,0,0,0}; // unit cell lattice (A)
	double reciprocal_basis[3][3] = {0,0,0,0,0,0,0,0,0}; // reciprocal space lattice (1/A)

	PeriodicBoundary & operator=( const PeriodicBoundary &rhs );

	void   update();
	double compute_volume();      // takes the determinant of the basis matrix
	void   compute_reciprocal();  // computes the reciprocal space basis
	double compute_cutoff();      // calculates the min cutoff radius from the basis lattice (AKA shortest vector problem)
	void   printboxdim();         // outputs data about the simulation box to stdout

};

