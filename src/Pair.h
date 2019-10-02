#pragma once

#include "Atom.h"
#include "Molecule.h"

// Each atom in a molecule has a pair list, which holds data relevant to each pair of atoms in the 
// system (so that it doesn't need to be recalculated for each simulation step--or perhaps a flag 
// indicating that the energetics for this pair *never* needs to be computed). If one of the pairs
// moves, the pair is flagged for recalculation, indicating that (some of) the stored data is no
// longer valid and needs to be recomputed.

// The first atom in the first molecule of the system will have a pair list with a pair entry for 
// every other atom in the system (across all other molecules). The 2nd atom in the system (i.e. 
// the 2nd atom in molecule 1, or [if "molecule" 1 only had 1 atom] atom 1 in molecule 2) will 
// have a pair list that contains a pair entry for every atom in the system, except atom 1. The 
// pairing of atom 1 with atom 2 is covered in the pair list of atom 1.

class Pair
{
public:
	Pair() {
		frozen               = 0;
		rd_excluded          = 0;
		es_excluded          = 0;
		attractive_only      = 0;
		recalculate_energy   = 0;
		lrc                  = 0;
		last_volume          = 0;
		epsilon, sigma       = 0;
		r                    = 0;
		rimg                 = 0;
		rd_energy            = 0;
		es_real_energy       = 0;
		es_self_intra_energy = 0;
		sigrep               = 0;
		c6                   = 0;
		c8                   = 0;
		c10                  = 0;

		for (int i = 0; i < 3; i++) {
			dimg  [i] = 0;
			d_prev[i] = 0;
		}

		atom     = nullptr;
		molecule = nullptr;
		next     = nullptr;

	};
	Pair(const Pair &other) {
		frozen               = other.frozen;
		rd_excluded          = other.rd_excluded;
		es_excluded          = other.es_excluded;
		attractive_only      = other.attractive_only;
		recalculate_energy   = other.recalculate_energy;
		lrc                  = other.lrc;
		last_volume          = other.last_volume;
		epsilon, sigma       = other.epsilon, sigma;
		r                    = other.r;
		rimg                 = other.rimg;
		rd_energy            = other.rd_energy;
		es_real_energy       = other.es_real_energy;
		es_self_intra_energy = other.es_self_intra_energy;
		sigrep               = other.sigrep;
		c6                   = other.c6;
		c8                   = other.c8;
		c10                  = other.c10;

		for (int i = 0; i < 3; i++) {
			dimg  [i] = other.dimg  [i];
			d_prev[i] = other.d_prev[i];
		}

		atom     = other.atom;
		molecule = other.molecule;

		next     = nullptr;
	};
	~Pair(){};

	
	int      frozen,               //are they both MOF atoms, for instance
	         rd_excluded, 
	         es_excluded,
	         attractive_only,
	         recalculate_energy;
	double   lrc,                  //LJ long-range correction
	         last_volume,          //what was the volume when we last calculated LRC? needed for NPT
	         epsilon, sigma,       //LJ
	         r, rimg, dimg[3],     //separation and separation with nearest image
	         d_prev[3],            //last known position
	         rd_energy, 
	         es_real_energy,
	         es_self_intra_energy,
	         sigrep,
	         c6, c8, c10;
	Atom     * atom;               // the other atom in the pairing
	Molecule * molecule;           // the molecule to which the other atom belongs
	Pair     * next;
	
};

