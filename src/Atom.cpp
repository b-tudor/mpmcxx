#include "Atom.h"
#include "Pair.h"
#include "SafeOps.h"

#include <cstring>


Atom::Atom() {}




Atom::Atom(const Atom &other) : 
	id(                   other.id), 
	frozen(               other.frozen),
	adiabatic(            other.adiabatic),
	spectre(              other.spectre),
	target(               other.target),
	mass(                 other.mass),
	charge(               other.charge),
	polarizability(       other.polarizability),
	epsilon(              other.epsilon),
	sigma(                other.sigma),
	omega(                other.omega),
	c6(                   other.c6),
	c8(                   other.c8),
	c10(                  other.c10),
	c9(                   other.c9),
	es_self_point_energy( other.es_self_point_energy),
	dipole_rrms(          other.dipole_rrms),
	rank_metric(          other.rank_metric),
	gwp_alpha(            other.gwp_alpha),
	lrc_self(             other.lrc_self),
	last_volume(          other.last_volume),
	gwp_spin(             other.gwp_spin),
	pairs(                nullptr),
	next(                 nullptr)
{
	
	strcpy(atomtype, other.atomtype);

	//\\///////////////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++) {                          
		pos[i] = other.pos[i];
		wrapped_pos[i] = other.wrapped_pos[i];
		ef_static[i] = other.ef_static[i];
		ef_static_self[i] = other.ef_static_self[i];
		ef_induced[i] = other.ef_induced[i];
		ef_induced_change[i] = other.ef_induced_change[i];
		mu[i] = other.mu[i];
		old_mu[i] = other.old_mu[i];
		new_mu[i] = other.new_mu[i];
	}//\\//////////////////////////////////////////////////////////////
	

	// Copy the pair list, if one exists
	////////////////////////////////////////////////////////////////////
			
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




Atom::~Atom() {
	free_pairs();
}





void Atom::free_pairs() {
	if(pairs)
		recursive_free_pairs( pairs );
}
void Atom::recursive_free_pairs(Pair * &pr) {
	
	if( pr->next )
		recursive_free_pairs( pr->next );

	delete pr;
	pr = nullptr;
}