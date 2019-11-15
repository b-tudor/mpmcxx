#include "Atom.h"
#include "Pair.h"
#include "SafeOps.h"

#include <cstring>


Atom::Atom()
{
	
	id                   = 0;
	bond_id              = 0;
	atomtype[0]          = (char) 0;
	frozen               = 0;
	adiabatic            = 0;
	spectre              = 0;
	target               = 0;
	mass                 = 0.0;
	charge               = 0.0;
	polarizability       = 0.0;
	epsilon              = 0.0;
	sigma                = 0.0;
	omega                = 0.0;
	c6                   = 0;
	c8                   = 0;
	c10                  = 0;
	c9                   = 0;
	dipole_rrms          = 0.0;
	rank_metric          = 0.0;
	gwp_spin             = 0;
	gwp_alpha            = 0.0;
	site_neighbor_id     = 0; // dr fluctuations will be applied along the vector from this atom to the atom identified by this variable
	lrc_self             = 0.0;
	last_volume          = 0.0; // currently only used in disp_expansion.c
	es_self_point_energy = 0.0;

	for( int i=0; i<3; i++ ) {
		pos              [i] = 0.0;
		wrapped_pos      [i] = 0.0; //absolute and wrapped (into main unit cell) position
		ef_static        [i] = 0.0;
		ef_static_self   [i] = 0.0;
		ef_induced       [i] = 0.0;
		ef_induced_change[i] = 0.0;
		mu               [i] = 0.0;
		old_mu           [i] = 0.0;
		new_mu           [i] = 0.0;
	}
	

	// pair_t *pairs;

}





Atom::Atom(const Atom &other) {
	
	strcpy(atomtype, other.atomtype);

	id        = other.id;
	bond_id   = other.bond_id;
	frozen    = other.frozen;
	adiabatic = other.adiabatic;
	spectre   = other.spectre;
	target    = other.target;
	mass      = other.mass;
	charge    = other.charge;
	epsilon   = other.epsilon;
	sigma     = other.sigma;
	omega     = other.omega;
	c6        = other.c6;
	c8        = other.c8;
	c10       = other.c10;
	c9        = other.c9;
	polarizability           = other.polarizability;
	es_self_point_energy     = other.es_self_point_energy;
	dipole_rrms              = other.dipole_rrms;
	rank_metric              = other.rank_metric;
	gwp_alpha                = other.gwp_alpha;
	lrc_self                 = other.lrc_self;
	last_volume              = other.last_volume;
	gwp_spin                 = other.gwp_spin;
	site_neighbor_id         = other.site_neighbor_id;
	
	for (int i = 0; i < 3; i++) {
		pos[i]               = other.pos[i];
		wrapped_pos[i]       = other.wrapped_pos[i];
		ef_static[i]         = other.ef_static[i];
		ef_static_self[i]    = other.ef_static_self[i];
		ef_induced[i]        = other.ef_induced[i];
		ef_induced_change[i] = other.ef_induced_change[i];
		mu[i]                = other.mu[i];
		old_mu[i]            = other.old_mu[i];
		new_mu[i]            = other.new_mu[i];
	}
	
	
	next = nullptr;



	// Copy the pair list
	////////////////////////////////////////////////////////////////////
	{
		
		Pair *other_pair_ptr = other.pairs;   // position in the 'other' pair list

		if (other_pair_ptr) 
		// If there is a list... 
		{
			pairs = new Pair(*other_pair_ptr); // Create the first node for our new pairs list.
			Pair *this_pair_ptr = pairs;       // iterator for the list we are building

			// Traverse the 'other' pair list...
			for (other_pair_ptr = other_pair_ptr->next; other_pair_ptr; other_pair_ptr = other_pair_ptr->next) {
				// ...creating new nodes in our list for each pair node we discover in the 'other' list.
				this_pair_ptr->next = new Pair(*other_pair_ptr);
				this_pair_ptr = this_pair_ptr->next;
			}
			this_pair_ptr->next = nullptr; // Terminate our new list.

		} else {
			// If there is no pair list on the original, we leave the new pair list empty as well
			pairs = nullptr;
		}
	}
	
	
}




Atom::~Atom() {
	free_pairs();
}





void Atom::free_pairs() {

	if(pairs)
		recursive_free_pairs( pairs );

	pairs = nullptr;
	return;
}
void Atom::recursive_free_pairs(Pair *pPair) {
	
	if( pPair->next )
		recursive_free_pairs( pPair->next );

	SafeOps::free(pPair);
	return;
}