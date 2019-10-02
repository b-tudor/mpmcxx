// @2007, Jonathan Belof
// Space Research Group
// Department of Chemistry
// University of South Florida


#include "Atom.h"
#include "Molecule.h"
#include "Pair.h"
#include "SafeOps.h"
#include "System.h"



// Allocate a series of linked lists cataloging pairs of sites, one entry per pair
//   E.g., pair linked list for A, B, C & D:
//     pairs[0]: B-> C-> D
//     pairs[1]: C-> D
//     pairs[2]: D

void System::allocate_pair_lists() {

	Pair *pair_ptr, *prev_pair_ptr;

	natoms = countNatoms();

	// build atom and molecule arrays for easy-access references
	rebuild_arrays();
	int n = natoms;

	// setup the pairs, top-right triangle (analogous to an upper/lower triangle matrix)
	for(int i = 0; i < (n - 1); i++) {

		SafeOps::calloc( atom_array[i]->pairs, 1, sizeof(Pair), __LINE__, __FILE__ );
		pair_ptr = atom_array[i]->pairs;
		prev_pair_ptr = pair_ptr;

		for(int j = (i + 1); j < n; j++) {
			SafeOps::calloc(pair_ptr->next, 1, sizeof(Pair), __LINE__, __FILE__ );
			prev_pair_ptr = pair_ptr;
			pair_ptr = pair_ptr->next;
		}

		prev_pair_ptr->next = nullptr;
		free(pair_ptr);
	}
}




// add new pairs for when a new molecule is created */
void System::update_pairs_insert() {

	int         n = 0;
	Molecule  * molecule_ptr;
	Atom      * atom_ptr;
	Pair      * pair_ptr;


	// count the number of atoms per molecule
	for( atom_ptr = checkpoint->molecule_altered->atoms; atom_ptr; atom_ptr = atom_ptr->next, n++);

	// add n number of pairs to altered and all molecules ahead of it in the list
 	for(molecule_ptr = molecules; molecule_ptr != checkpoint->tail; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			// go to the end of the pair list
			if(atom_ptr->pairs) {

				// go to the end of the pair list
				for(pair_ptr = atom_ptr->pairs; pair_ptr->next; pair_ptr = pair_ptr->next);

				// tag on the extra pairs
				for(int i = 0; i < n; i++) {
					SafeOps::calloc( pair_ptr->next, 1, sizeof(Pair), __LINE__, __FILE__);
					pair_ptr = pair_ptr->next;
				}

			} else {

				// needs a new list
				SafeOps::calloc( atom_ptr->pairs, 1, sizeof(Pair), __LINE__, __FILE__ );
				pair_ptr = atom_ptr->pairs;
				for( int i = 0; i < (n - 1); i++) {
					SafeOps::calloc( pair_ptr->next, 1, sizeof(Pair), __LINE__, __FILE__ );
					pair_ptr = pair_ptr->next;
				}

			}

		} // for atom 
	} // for molecule
}




// remove pairs when a molecule is deleted
void System::update_pairs_remove() {

	// When removing a molecule from the system, the number of atom-atom pairs in the system will
	// be reduced and so the quantity of Pair nodes required to document these pairings is likewise reduced. 
	// This function removes nodes from the atomic Pair lists of those atoms that will require fewer Pairs post-
	// molecule-removal. The references stored in each pair list are no longer valid as references to the removed
	// atoms are kept and references to extant atoms are deleted. Only the quantity of Pair nodes is
	// adjusted, and these nodes will have to be re-populated with valid data at a later time. 

	int         n            = 0,       // the number of atoms (or sites) in the removed molecule
		        nPairs       = 0;       // A counter for Pair nodes in a given Pair list
	Molecule  * molecule_ptr = nullptr;
	Atom      * atom_ptr     = nullptr;
	Pair      * pair_ptr     = nullptr,
	         ** pair_array   = nullptr;


	// count the number of atoms in the "backup" molecule, store result in n
	for( atom_ptr = checkpoint->molecule_backup->atoms, n=0; atom_ptr; atom_ptr = atom_ptr->next)
		++n;

	// remove n Pair nodes, from the Pair list for all molecules appearing before the Pair nodes
	// that reference the atoms in the molecule to be removed. 

	// Create an array for easy reference to the Pair nodes for the Pair list in question
	SafeOps::calloc( pair_array, 1, sizeof(Pair *), __LINE__, __FILE__ );

	for( molecule_ptr = molecules; molecule_ptr != checkpoint->tail; molecule_ptr = molecule_ptr->next) {
	// We only examine the molecules that appear (or appeared) before the molecule that is to be (or was) removed
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for( pair_ptr = atom_ptr->pairs, nPairs = 0; pair_ptr; pair_ptr = pair_ptr->next, nPairs++) {
				// update array length to reflect most current node count
				SafeOps::realloc( pair_array, sizeof(Pair *)*(nPairs + 1), __LINE__, __FILE__ ); 
				pair_array[nPairs] = pair_ptr; 
			}

			// delete n nodes from the end of the Pair list
			for(int i = (nPairs - n); i < nPairs; i++)
				free(pair_array[i]);

			// Terminate the newly truncated Pair list
			if((nPairs - n) > 0)
				pair_array[(nPairs - n - 1)]->next = nullptr;
			else
				atom_ptr->pairs = nullptr;

		}
	}

	// free our temporary array
	free(pair_array);

}