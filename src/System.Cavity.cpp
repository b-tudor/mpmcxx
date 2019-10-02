// @2007, Jonathan Belof
// Space Research Group
// Department of Chemistry
// University of South Florida

#include "Atom.h"
#include "Molecule.h"
#include "Pair.h"
#include "System.h"

#include "SafeOps.h"


// create a 3D histogram of atoms lying within a sphere centered at each grid point 
void System::cavity_update_grid() {

	Molecule *molecule_ptr;
	Atom     *atom_ptr;
	double    grid_component[3];
	double    grid_vector[3];
	double    r;

	
	// clear the grid 
	for(int i = 0; i < cavity_grid_size; i++ ) {

		for( int j = 0; j < cavity_grid_size; j++ ) {

			for( int k = 0; k < cavity_grid_size; k++ ) {

				cavity_grid[i][j][k].occupancy = 0;
				for( int p = 0; p < 3; p++ )
					cavity_grid[i][j][k].pos[p] = 0;
			}
		}
	}

	// loop over the grid, bin a sphere when needed 
	for( int i = 0;   i < cavity_grid_size;   i++) {
		for( int j = 0;   j < cavity_grid_size;   j++) {
			for( int k = 0;   k < cavity_grid_size;   k++) {

				// divide up each grid component
				grid_component[0] = ((double)(i + 1))  /  ((double)(cavity_grid_size + 1));
				grid_component[1] = ((double)(j + 1))  /  ((double)(cavity_grid_size + 1));
				grid_component[2] = ((double)(k + 1))  /  ((double)(cavity_grid_size + 1));

				// project the grid point onto our actual basis
				for(int p = 0;   p < 3;   p++ ) {
					grid_vector[p] = 0;
					for( int q = 0; q < 3; q++)
						grid_vector[p] += pbc.basis[q][p] * grid_component[q];
				}

				// put into real coordinates 
				for( int p = 0;   p < 3;   p++ )
					for( int q = 0;   q < 3;   q++ )
						grid_vector[p] -= 0.5 * pbc.basis[q][p];

				// if an atomic coordinate lies within a sphere centered on the grid point, then bin it 
				for( molecule_ptr = molecules;   molecule_ptr;   molecule_ptr = molecule_ptr->next ) {
					for( atom_ptr = molecule_ptr->atoms;   atom_ptr;   atom_ptr = atom_ptr->next ) {

						// get the displacement from the grid point
						r = 0;
						for(int p = 0; p < 3; p++)
							r += (grid_vector[p] - atom_ptr->wrapped_pos[p]) * (grid_vector[p] - atom_ptr->wrapped_pos[p]);
						r = sqrt(r);

						// inside the sphere?
						if( r < cavity_radius ) ++cavity_grid[i][j][k].occupancy;

					} // for atom 
				} // for molecule 

				// store the location of this grid point 
				for( int p = 0;  p < 3;  p++ )
					cavity_grid[i][j][k].pos[p] = grid_vector[p];

			} // for k 
		} // for j 
	} // for i 

	// update the cavity insertion probability estimate
	update_cavity_probability();
	// update the accessible insertion volume
	update_cavity_volume();

}



// probability of finding an empty cavity on the grid
void System::update_cavity_probability() {

	//int i, j, k;
	double probability  = 0;
	int    total_points = 0;

	// total number of potential cavities */
	total_points = cavity_grid_size * cavity_grid_size * cavity_grid_size;

	// find the number of open cavities
	cavities_open = 0;
	for( int i = 0;   i < cavity_grid_size;   i++)
		for( int j = 0;   j < cavity_grid_size;   j++)
			for( int k = 0;   k < cavity_grid_size;   k++)
				if( ! cavity_grid[i][j][k].occupancy )
					++cavities_open;

	// the overall probability ratio
	probability = ((double) cavities_open)  /  ((double) total_points);

	// update the observable
	nodestats->cavity_bias_probability = probability;
}




// total volume of accessible cavities via Monte Carlo integration 
void System::update_cavity_volume() {

	const double DARTSCALE = 0.1;
	int    hits = 0;
	int    num_darts = 0;
	double pos_vec[3] = { 0 };
	double grid_vec[3] = { 0 };
	double fraction_hits = 0;

	// good rule of thumb is 1 per 10 A^3 
	num_darts = (int)( pbc.volume * DARTSCALE );

	// throw random darts and count the number of hits
	for( int dart = 0; dart < num_darts; dart++ ) {

		// generate a random grid position 
		for( int p = 0; p < 3; p++ )
			grid_vec[p] = -0.5 + get_rand();

		// zero the coordinate vector
		for( int p = 0; p < 3; p++)
			pos_vec[p] = 0;

		// linear transform vector into real coordinates
		for( int p = 0; p < 3; p++ )
			for( int q = 0; q < 3; q++ )
				pos_vec[p] += pbc.basis[q][p] * grid_vec[q];

		// check if the random point lies within an empty cavity
		if(   is_point_empty(pos_vec[0], pos_vec[1], pos_vec[2])  )
			++hits;

	}

	// determine the percentage of free cavity space 
	fraction_hits = ((double)hits)/((double)num_darts);
	cavity_volume = fraction_hits * pbc.volume;	// normalize w.r.t. the cell volume

}



// check whether a point (x,y,z) lies within an empty cavity if so, return 1 
bool System::is_point_empty( double x, double y, double z) {

	double r = 0;

	for( int i = 0;   i < cavity_grid_size;   i++) {
		for( int j = 0;   j < cavity_grid_size;   j++) {
			for( int k = 0;   k < cavity_grid_size;   k++) {

				if( ! cavity_grid[i][j][k].occupancy )
				{
					r  = pow(  (x - cavity_grid[i][j][k].pos[0]),   2  );
					r += pow(  (y - cavity_grid[i][j][k].pos[1]),   2  );
					r += pow(  (z - cavity_grid[i][j][k].pos[2]),   2  );
					r  = sqrt(r);
					if( r < cavity_radius )
						return true;
				}

			} // end i 
		} // end j 
	} // end k

	return false;
}



// allocate the grid 
void System::setup_cavity_grid() {

	SafeOps::calloc( cavity_grid, cavity_grid_size, sizeof(cavity_t **), __LINE__, __FILE__ );
	
	for( int i = 0; i < cavity_grid_size; i++ ) {

		SafeOps::calloc( cavity_grid[i], cavity_grid_size, sizeof(cavity_t *), __LINE__, __FILE__ );
		
		for( int j = 0; j < cavity_grid_size; j++ ) 

			SafeOps::calloc( cavity_grid[i][j], cavity_grid_size, sizeof(cavity_t), __LINE__, __FILE__ );

	}
}



//check cavity_autoreject_absolute -- probably not the most efficient place to put this, but likely the safest
double System::cavity_absolute_check()
{
	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;
	
	for( molecule_ptr=molecules;  molecule_ptr;  molecule_ptr=molecule_ptr->next ) {
		for( atom_ptr=molecule_ptr->atoms;  atom_ptr;  atom_ptr=atom_ptr->next ) {
			for( pair_ptr=atom_ptr->pairs;  pair_ptr;  pair_ptr=pair_ptr->next ) {
				if( molecule_ptr == pair_ptr->molecule ) 
					continue; //skip if on the same molecule
				if( pair_ptr->rimg < cavity_autoreject_scale ) 
					return MAXVALUE;
			}
		}
	}
	return 0;
}