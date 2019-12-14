#pragma once
#include "constants.h"

#include "Atom.h"
#include "Rando.h"
#include "Vector3D.h"

class PeriodicBoundary;

class Molecule
{
public:

	Molecule();
	Molecule( const  Molecule &orig );

	~Molecule();

	// void free_pairs();
	void rotate_rand(double scale); //, const PeriodicBoundary &pbc, std::mt19937 *mt_rand );
	void rotate( double x, double y, double z, double angle ); // params are axis-of-rotation and angle-of-rotation (in degrees)
	void orient(Vector3D orientation_vector, int orientation_site );
	void translate_rand_pbc( double scale, const PeriodicBoundary &pbc, std::mt19937 *mt_rand );
	void translate_rand_pbc( double scale, const PeriodicBoundary &pbc, double dice[6]);   // dice are 6 random numbers [0..1]
	inline void translate(const Vector3D v) { translate(v.x(), v.y(), v.z()); };
	void translate(const double x, const double y, const double z ); // translate a molecule's COM to (x,y,z)
	void move_to_( double x, double y, double z ); // move a molecule's COM to (x,y,z)
	void displace_gwp( double scale, std::mt19937 *mt_rand );
	void update_COM();
	void free_atoms(); // Free memory used by the molecule's linked list of atoms. 
	void wipe_pair_refs(); // zero out all the pair list pointers while keeping the lists themselves
	                       // (for use prior to deletion, when the pair list is to be preserved elsewhere)
private:
	void recursive_free_atoms(Atom * &atom);



public:
	int        id;
	char       moleculetype[maxID];
	double     mass;
	int        frozen, 
	           adiabatic,
	           spectre,
	           target;
	double     com[3],
	           wrapped_com[3]; //center of mass
	double     iCOM       [3]; // initial Center of Mass
	int        nuclear_spin;
	double     rot_partfunc_g, 
	           rot_partfunc_u,
	           rot_partfunc;
	Atom     * atoms;
	Molecule * next;


	#ifdef QM_ROTATION
		double     * quantum_rotational_energies;
		complex_t ** quantum_rotational_eigenvectors;
		int        * quantum_rotational_eigensymmetry;
		double       quantum_rotational_potential_grid[QUANTUM_ROTATION_GRID][QUANTUM_ROTATION_GRID];
	#endif // QM_ROTATION
	#ifdef XXX
		// XXX - vib work in progress
		double     * quantum_vibrational_energies;
		complex_t ** quantum_vibrational_eigenvectors;
		int        * quantum_vibrational_eigensymmetry;
	#endif // XXX
};

