#include <cstring>
#include <math.h>
#include <vector>

#include "Molecule.h"
#include "PeriodicBoundary.h"
#include "Pair.h"
#include "Quaternion.h"
#include "Rando.h"
#include "SafeOps.h"





Molecule::Molecule(){ }
Molecule::Molecule( const Molecule &other ) :
	id(             other.id),
	mass(           other.mass),
	frozen(         other.frozen),
	adiabatic(      other.adiabatic),
	spectre(        other.spectre),
	target(         other.target),
	nuclear_spin(   other.nuclear_spin),
	rot_partfunc_g( other.rot_partfunc_g),
	rot_partfunc_u( other.rot_partfunc_u),
	rot_partfunc(   other.rot_partfunc),
	atoms(          nullptr),
	next(           nullptr)
{
	strcpy(moleculetype, other.moleculetype);

	for( int i=0; i<3; i++ ) {
		com        [i]  =  other.com        [i];
		wrapped_com[i]  =  other.wrapped_com[i];
	}

	// Copy the atom list
	Atom *atom_ptr = other.atoms;

	if ( atom_ptr ) {
		// If the molecule has an atom list (it should), then we copy them
		atoms = new Atom(*atom_ptr); // create the first atom in this molecule...
		Atom* a = atoms;             // ... and create a pointer to traverse the list we are creating

		for(  atom_ptr=atom_ptr->next;  atom_ptr;  atom_ptr=atom_ptr->next  ) {  
			a->next = new Atom(*atom_ptr); // copy the next atom in the other list
			a = a->next;                   // move to the next position in this list
		}
		a->next = nullptr;
	}


	#ifdef QM_ROTATION
		// untouched--needs to be converted - bt
		int i, j;
		if(system->quantum_rotation) {
			dst->quantum_rotational_energies = calloc(system->quantum_rotation_level_max, sizeof(double));
			memnullcheck(dst->quantum_rotational_energies, system->quantum_rotation_level_max*sizeof(double),__LINE__-1, __FILE__);
			dst->quantum_rotational_eigenvectors = calloc(system->quantum_rotation_level_max, sizeof(complex_t *));
			memnullcheck(dst->quantum_rotational_eigenvectors,system->quantum_rotation_level_max*sizeof(complex_t *),__LINE__-1, __FILE__);
			for(i = 0; i < system->quantum_rotation_level_max; i++) {
				dst->quantum_rotational_eigenvectors[i] = calloc((system->quantum_rotation_l_max + 1)*(system->quantum_rotation_l_max + 1), sizeof(complex_t));
				memnullcheck(dst->quantum_rotational_eigenvectors[i],(system->quantum_rotation_l_max+1)*(system->quantum_rotation_l_max+1)*sizeof(complex_t),__LINE__-1, __FILE__);
			}
			dst->quantum_rotational_eigensymmetry = calloc(system->quantum_rotation_level_max, sizeof(int));
			memnullcheck(dst->quantum_rotational_eigensymmetry,system->quantum_rotation_level_max*sizeof(int),__LINE__-1, __FILE__);

			std::memcpy(dst->quantum_rotational_energies, src->quantum_rotational_energies, system->quantum_rotation_level_max*sizeof(double));
			for(i = 0; i < system->quantum_rotation_level_max; i++) {
				for(j = 0; j < (system->quantum_rotation_l_max + 1)*(system->quantum_rotation_l_max + 1); j++) {
					dst->quantum_rotational_eigenvectors[i][j].real = src->quantum_rotational_eigenvectors[i][j].real;
					dst->quantum_rotational_eigenvectors[i][j].imaginary = src->quantum_rotational_eigenvectors[i][j].imaginary;
				}
			}
			std::memcpy(dst->quantum_rotational_eigensymmetry, src->quantum_rotational_eigensymmetry, system->quantum_rotation_level_max*sizeof(int));
		}
	#endif // QM_ROTATION
}




Molecule::~Molecule() {
	
	#ifdef QM_ROTATION
		// untouched--needs to be converted - bt
		int i;
		if(system->quantum_rotation && !molecule->frozen) {

			SafeOps::free(molecule->quantum_rotational_energies);
			SafeOps::free(molecule->quantum_rotational_eigensymmetry);
			for(i = 0; i < system->quantum_rotation_level_max; i++)
				SafeOps::free(molecule->quantum_rotational_eigenvectors[i]);
			SafeOps::free(molecule->quantum_rotational_eigenvectors);

		}
	#endif // QM_ROTATION

	
	free_atoms();

}
void Molecule::free_atoms() {

	if (atoms)
		recursive_free_atoms(atoms);

	atoms = nullptr;
}
void Molecule::recursive_free_atoms(Atom* &atom) {

	if (atom->next)
		recursive_free_atoms(atom->next);

	atom->free_pairs();
	delete atom;
	atom = nullptr;
}
void Molecule::wipe_pair_refs() {
	for (Atom* pAtom = atoms; pAtom; pAtom = pAtom->next)
		pAtom->pairs = nullptr;
}



// perform a general random rotation now with quaternions - AH 
void Molecule::rotate_rand( double scale ) {    
	
	// To avoid bias, random numbers for rotate_rand(double, double[4]) must be pulled from a normal distribution.
	double x     = Rando::rand_normal();
	double y     = Rando::rand_normal();
	double z     = Rando::rand_normal();
	double angle = Rando::rand() * 360 * scale;
	
	rotate( x, y, z, angle );
}
void Molecule::rotate( double x, double y, double z, double angle ) {

	Atom   * atom_ptr        = nullptr;
	double * new_coord_array = nullptr;
	         
	Quaternion rand_rotation( x, y, z, angle, Quaternion::AXIS_ANGLE_DEGREE );
	Quaternion rand_rotation_conjugate = rand_rotation.conjugate(); // needed to transform coordinates
	
	// count the number of atoms in a molecule, and allocate new coords array 
	int n=0;
	for( atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next)
		++n;
	SafeOps::calloc( new_coord_array, (size_t)(n)*3, sizeof(double), __LINE__, __FILE__ );
	
	// translate the molecule to the origin 
	for(atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] -= com[0];
		atom_ptr->pos[1] -= com[1];
		atom_ptr->pos[2] -= com[2];
	}

	// quaternion multiply
	atom_ptr = atoms;
	for( int i=0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		int ii = i*3;
    

		//position_vector = position
		Quaternion position_vector( atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2], 0.0, Quaternion::XYZW );
		Quaternion answer = rand_rotation * (position_vector * rand_rotation_conjugate);


		// Original treatment: 
		Quaternion answer2 = position_vector * rand_rotation_conjugate;
		answer2 = rand_rotation * answer2;
		// Check against original treatment
		if ((fabs(answer.x() - answer2.x()) > 0.00001) && (fabs(answer.y() - answer2.y()) > 0.00001) && (fabs(answer.z() - answer2.z()) > 0.00001)) {
			Output::out("QUAT EXPRESSION WAS NOT THE SAME AS ORIGINAL!\n");
			throw internal_error;
		}
		


		//set the new coords
		new_coord_array[ii+0] = answer.x();
		new_coord_array[ii+1] = answer.y();
		new_coord_array[ii+2] = answer.z();
	}

	// set the new coordinates and then translate back from the origin
	atom_ptr = atoms;
	for( int i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		int ii = i*3;
		atom_ptr->pos[0] = new_coord_array[ii+0];
		atom_ptr->pos[1] = new_coord_array[ii+1];
		atom_ptr->pos[2] = new_coord_array[ii+2];

		atom_ptr->pos[0] += com[0];
		atom_ptr->pos[1] += com[1];
		atom_ptr->pos[2] += com[2];

	}

	// free our temporary array
	SafeOps::free(new_coord_array);

}




void Molecule::orient(Vector3D orientation, int orientation_site) {

	
	update_COM();
	Vector3D rCOM(com[0], com[1], com[2]);
	translate(-rCOM);

	// Get current orientation
	
	int site = 0;
	Atom *aPtr = nullptr;
	// Find the molecular site through which orientation will be determined
	for (aPtr = atoms; site != orientation_site; aPtr = aPtr->next)
		site++;
	Vector3D current_orientation(aPtr->pos[0], aPtr->pos[1], aPtr->pos[2]);


	// Rotate the orientation site through the center of mass (COM) so that it lies along the vector orientation

	// Compute the rotation angle
	current_orientation.normalize();
	double angle = acos((current_orientation * orientation) / (orientation.norm()));
	

	// Construct the rotation axis
	Vector3D rotationAxis = current_orientation.cross(orientation);

	// Construct the Quaterion that will perform a rotation of "angle" radians about "rotation axis"
	Quaternion rotation(rotationAxis.x(), rotationAxis.y(), rotationAxis.z(), angle, Quaternion::AXIS_ANGLE_RADIAN);

	// Apply the rotation to each site in this image of the molecule
	for (aPtr = atoms; aPtr; aPtr = aPtr->next) {
		Vector3D site_pos(aPtr->pos[0], aPtr->pos[1], aPtr->pos[2]);
		site_pos = rotation.rotate(site_pos);
		aPtr->pos[0] = site_pos.x();
		aPtr->pos[1] = site_pos.y();
		aPtr->pos[2] = site_pos.z();
	}

	// Move the molecule back to where it was
	translate(rCOM);

	return;
}



// update a molecules center of mass (COM)
void Molecule::update_COM() {

	Atom *atom_ptr = atoms;

	mass = 0; 

	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	
		
	for( atom_ptr=atoms; atom_ptr; atom_ptr=atom_ptr->next ) {

		mass += atom_ptr->mass;
		com[0] += atom_ptr->mass * atom_ptr->pos[0];
		com[1] += atom_ptr->mass * atom_ptr->pos[1];
		com[2] += atom_ptr->mass * atom_ptr->pos[2];
	}

	com[0] /= mass;
	com[1] /= mass;
	com[2] /= mass;
}



// perform a general random translation  (formerly "translate(...)" )
void Molecule::translate_rand_pbc( double scale,  const PeriodicBoundary &pbc, std::mt19937 *mt_rand ) {
	
	double dice_rolls[6];
	std::uniform_real_distribution<double> dist{ 0,1 };

	for (int i = 0; i < 6; i++)
		dice_rolls[i] = dist(*mt_rand);

	translate_rand_pbc( scale, pbc, dice_rolls );
}
void Molecule::translate_rand_pbc(double scale, const PeriodicBoundary &pbc, double dice[6] ) {

	Atom   * atom_ptr;
	double trans_x = 0;
	double trans_y = 0;
	double trans_z = 0;

	trans_x = scale * dice[0] * pbc.cutoff;
	trans_y = scale * dice[1] * pbc.cutoff;
	trans_z = scale * dice[2] * pbc.cutoff;

	if( dice[3] < 0.5 ) 
		trans_x *= -1.0;
	if( dice[4] < 0.5 ) 
		trans_y *= -1.0;
	if( dice[5] < 0.5 ) 
		trans_z *= -1.0;

	for (atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] += trans_x;
		atom_ptr->pos[1] += trans_y;
		atom_ptr->pos[2] += trans_z;
	}

	update_COM();
}



// move a molecule's COM to (x,y,z)
void Molecule::move_to_(double x, double y, double z) {
	translate(  x-com[0],  y-com[1],  z-com[2]  );
}



// translate a molecule's COM by (x,y,z)
void Molecule::translate(const double x, const double y, const double z) {
	
	com[0] += x;
	com[1] += y;
	com[2] += z;

	Atom   * atom_ptr = nullptr;
	for (atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] += x;
		atom_ptr->pos[1] += y;
		atom_ptr->pos[2] += z;
	}
}



// perform a perturbation to a gaussian width
void Molecule::displace_gwp( double scale, std::mt19937 *mt_rand ) {

	std::uniform_real_distribution<double> dist{0,1};

	Atom   * atom_ptr;
	double   perturb  = 0;

	for(atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next) {

		if(atom_ptr->gwp_spin) {
			perturb = scale * ( dist(* mt_rand) - 0.5);
			atom_ptr->gwp_alpha += perturb;
			// make sure the width remains positive 
			atom_ptr->gwp_alpha = fabs(atom_ptr->gwp_alpha);
		}
	}
}
