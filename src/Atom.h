#pragma once
#include "constants.h"

class Pair;





class Atom
{
public:
	Atom();
	Atom(const Atom &other);
	~Atom();
	void free_pairs();
private:
	void recursive_free_pairs( Pair * &pr);

	
public:
	char   atomtype[maxLine]    = { 0 };
	int    id                   = 0;
	int    frozen               = 0, 
	       adiabatic            = 0,
	       spectre              = 0,
	       target               = 0;
	double mass                 = 0,
	       charge               = 0,
	       polarizability       = 0,
	       epsilon              = 0,
	       sigma                = 0,
	       omega                = 0,
	       c6                   = 0,
	       c8                   = 0,
	       c10                  = 0,
	       c9                   = 0,
	       es_self_point_energy = 0,
	       pos[3]               = { 0,0,0 }, // absolute position
	       wrapped_pos[3]       = { 0,0,0 }, // wrapped (into main unit cell) position
	       ef_static[3]         = { 0,0,0 },
	       ef_static_self[3]    = { 0,0,0 },
	       ef_induced[3]        = { 0,0,0 },
	       ef_induced_change[3] = { 0,0,0 },
	       mu[3]                = { 0,0,0 },
	       old_mu[3]            = { 0,0,0 },
	       new_mu[3]            = { 0,0,0 },
	       dipole_rrms          = 0,
	       rank_metric          = 0,
	       gwp_alpha            = 0,
	       lrc_self             = 0,
	       last_volume          = 0;
	int    gwp_spin             = 0;
	Pair   *pairs               = nullptr;
	Atom   *next                = nullptr;
};

