#include <cstring>
#include <iostream>

extern int rank, size;
extern bool mpi;
#ifdef _MPI
#include <mpi.h>
#endif

#include "constants.h"
#include "Pair.h"
#include "PeriodicBoundary.h"
#include "SafeOps.h"
#include "System.h"




 
bool System::mc() {
// Implements the Markov chain

	double  initial_energy  = 0, 
	        final_energy    = 0,
	        current_energy  = 0;
	
	if( sorbateCount > 1 ) SafeOps::calloc( sorbateGlobal, sorbateCount, sizeof(sorbateAverages_t), __LINE__, __FILE__ );
	if( cavity_bias      ) cavity_update_grid(); // update the grid for the first time 
	observables->volume = pbc.volume; // set volume observable
	initial_energy = mc_initial_energy();
	mpiData mpi = setup_mpi();
	
	// save the initial state 
	do_checkpoint();

	
	// main MC loop 
	for( step = 1; step <= numsteps; step++ ) {

		// restore the last accepted energy 
		initial_energy = observables->energy;

		// perturb the system 
		make_move();

		// calculate the energy change 
		final_energy = energy();

		#ifdef QM_ROTATION
			// solve for the rotational energy levels 
			if(system->quantum_rotation && (system->checkpoint->movetype == MOVETYPE_SPINFLIP))
				quantum_system_rotational_energies(system);
		#endif // QM_ROTATION 

		// treat a bad contact as a reject 
		if( !std::isfinite(final_energy) ) {
			observables->energy = MAXVALUE;
			nodestats->boltzmann_factor = 0;
		} 
		else 
			boltzmann_factor( initial_energy, final_energy );

		// Metropolis function 
		if(  (get_rand() < nodestats->boltzmann_factor)   &&   ! iterator_failed  ) {

			/////////// ACCEPT

			current_energy = final_energy;

			// checkpoint 
			do_checkpoint();
			register_accept();

			// SA 
			if( simulated_annealing )
			{
				if( simulated_annealing_linear )
				{
					temperature = temperature + (simulated_annealing_target - temperature)/(numsteps - step);
					if(  (numsteps - step)   ==   0   )
						temperature = simulated_annealing_target;
				}
				else
					temperature = simulated_annealing_target + (temperature - simulated_annealing_target)*simulated_annealing_schedule;
			}

		} else { 
			// REJECT ///////////////////////////////////////////////////////////////
			//reset the polar iterative failure flag and restore from last checkpoint
			current_energy = initial_energy; //used in parallel tempering
			iterator_failed = 0;
			restore();
			register_reject();
		}

		// perform parallel_tempering
		if ( (parallel_tempering) && ((step % ptemp_freq) == 0) )
			temper_system( current_energy );

		// track the acceptance_rate 
		track_ar( nodestats );

		// each node calculates its stats 
		update_nodestats( nodestats, avg_nodestats );


		// do this every correlation time, and at the very end 
		if(  !(step % corrtime)  ||  (step == numsteps)  ) 
			do_corrtime_bookkeeping(mpi);

	} // main loop 
	


	// output restart files for each node
	if(  write_molecules_wrapper(pqr_output) < 0  ) {
		Output::err( "MC: could not write final state to disk\n" );
		throw unknown_file_error;
	}

	if(!rank) {
		free( mpi.rcv_strct   );
		free( mpi.temperature );
	}

	if( sorbateCount > 1 ) 
		free( mpi.sinfo );
	
	free( mpi.snd_strct     );
	free( mpi.observables   );
	free( mpi.avg_nodestats );

	return ok;
}




void System::output_file_data() {

	// write the averages to stdout 
	if( fp_histogram )
		write_histogram( fp_histogram, grids->avg_histogram->grid );
	
	if( write_performance( step ) < 0) {
		Output::err( "MC: could not write performance data to stdout\n" );
		throw unknown_file_error;
	}
	if( write_averages() < 0) {
		Output::err( "MC: could not write statistics to stdout\n" );
		throw unknown_file_error;
	}
}




double System::mc_initial_energy() {

	// get the initial energy of the system 
	double init_energy = energy();

	#ifdef QM_ROTATION
		// solve for the rotational energy levels 
		if(system->quantum_rotation) quantum_system_rotational_energies(system);
	#endif // QM_ROTATION 

	// be a bit forgiving of the initial state 
	if( !std::isfinite(init_energy) ) 
		init_energy = observables->energy = MAXVALUE;

	return init_energy;
}




System::mpiData System::setup_mpi() {

	mpiData mpi;
	setup_mpi_dataStructs( mpi );

	// write initial observables to stdout and logs
	if( ! rank ) {

		open_files(); // open output files 

		calc_system_mass();
		// average in the initial values once  (we don't want to double-count the initial state when using MPI)
		update_root_averages( observables );
		// average in the initial sorbate values
		if ( sorbateCount > 1 ) {
			update_sorbate_info(); //local update
			update_root_sorb_averages( sorbateInfo ); //global update
		}
		// write initial observables exactly once
		if( fp_energy ) 
			write_observables( fp_energy, observables, temperature );
		if( fp_energy_csv ) 
			write_observables_csv( fp_energy_csv, observables, temperature );
		Output::out( "MC: initial values:\n" );
		write_averages();
	}

	return mpi;
}



void System::setup_mpi_dataStructs(int qty) { setup_mpi_dataStructs(mpi_data, qty); }
void System::setup_mpi_dataStructs( mpiData &md ) {	setup_mpi_dataStructs(md, size);}
void System::setup_mpi_dataStructs( mpiData& md, int qty ) {

	// allocate the statistics structures 
	SafeOps::calloc(md.observables, 1, sizeof(observables_t), __LINE__, __FILE__);
	SafeOps::calloc(md.avg_nodestats, 1, sizeof(avg_nodestats_t), __LINE__, __FILE__);

	// if multiple-sorbates, allocate sorbate statistics structs
	if (sorbateCount > 1) {
		SafeOps::calloc(md.sinfo, sorbateCount, sizeof(sorbateInfo_t), __LINE__, __FILE__);
		SafeOps::calloc(sorbateGlobal, sorbateCount, sizeof(sorbateAverages_t), __LINE__, __FILE__);
	}

	// compute MPI message size	(for reporting data to head node)
	md.msgsize = sizeof(observables_t) + sizeof(avg_nodestats_t);
	if (calc_hist)
		md.msgsize += n_histogram_bins * sizeof(int);
	if (sorbateCount > 1)
		md.msgsize += sorbateCount * sizeof(sorbateInfo_t);

	#ifdef _MPI
		if (mpi) {
			// MPI_Datatype msgtype;
			MPI_Type_contiguous(md.msgsize, MPI_BYTE, &msgtype);
			MPI_Type_commit(&msgtype);
		}
	#endif

	// allocate MPI structures 
	SafeOps::calloc(md.snd_strct, md.msgsize, 1, __LINE__, __FILE__);

	if ((!mpi) || (!rank)) {
		// These structs will be allocated on all systems[] in a single-threaded run
		SafeOps::calloc(   md.rcv_strct, qty,     md.msgsize, __LINE__, __FILE__ );
		SafeOps::calloc( md.temperature, qty, sizeof(double), __LINE__, __FILE__ ); //temperature list for parallel tempering
	}
}



 void System::do_checkpoint() {
// this function
//   (a) backs up the state to be restore upon rejection
//   (b) determines what move will be made next time make_move() is called
//   (c) picks a molecule to which the move in (a) will be applied
//   (d) creates pointers to said molecule's head/tail list locations (easing access later)
//   
//   IT CANNOT pick moves for coordinated ensembles (e.g. ENSEMBLE_NVT_GIBBS || ENSEMBLE_PATH_INTEGRAL_NVT)
//   so these ensembles must have their own dedicated checkpoint()-type functions.

	int          i_exchange                     = 0,
	             i_adiabatic                    = 0,
	             num_molecules_exchange         = 0,
	             num_molecules_adiabatic        = 0,
	             alt                            = 0,
	             altered                        = 0;
	double       num_molecules_adiabatic_double = 0;
	Molecule  ** ptr_array_exchange             = nullptr,
	          ** ptr_array_adiabatic            = nullptr,
	           * molecule_ptr                   = nullptr,
	           * prev_molecule_ptr              = nullptr;


	// save the current observables 
	/////////////////////////////////////////////////////////////////////////////////////////////
	std::memcpy( checkpoint->observables, observables, sizeof(observables_t) );

	// Count exchangeable and adiabatic molecules, then allocate an array whose size is 
	// determined by said count. Populate the array with pointers to the molecules that
	// are eligible to be displaced, removed or altered. 
	/////////////////////////////////////////////////////////////////////////////////////////////
	{
		num_molecules_exchange  = 0;
		num_molecules_adiabatic = 0;

		// count molecule types
		for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
			if(!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target)) 
				++num_molecules_exchange;
			if(molecule_ptr->adiabatic) 
				++num_molecules_adiabatic;
		}

		// allocate array mem
		if( num_molecules_exchange > 0 )
			SafeOps::calloc( ptr_array_exchange,  num_molecules_exchange,  sizeof(Molecule *), __LINE__, __FILE__ );
		if( num_molecules_adiabatic > 0 )
			SafeOps::calloc( ptr_array_adiabatic, num_molecules_adiabatic, sizeof(Molecule *), __LINE__, __FILE__ );

		// populate arrays with molecule pointers
		i_exchange = 0;
		i_adiabatic = 0;
		for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
			if(!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target)) {
				ptr_array_exchange[i_exchange] = molecule_ptr;
				++i_exchange;
			}

			if( molecule_ptr->adiabatic ) {
				ptr_array_adiabatic[i_adiabatic] = molecule_ptr;
				++i_adiabatic;
			}
		}
	}
	

	// Determine what kind of move to do 
	//////////////////////////////////////////////////////////////////////////////////////////////////
	switch( ensemble ) {
	
		case ENSEMBLE_UVT :	
	
			if( get_rand() < insert_probability ) {	// do a particle insertion/deletion move
	
				if(get_rand() < 0.5) {  
					// INSERT:
					checkpoint->movetype = MOVETYPE_INSERT;
				} else {
					// REMOVE:
					checkpoint->movetype = MOVETYPE_REMOVE;
				}
	
			} else if( quantum_rotation ) {
	
				if( get_rand() < spinflip_probability )
					// SPINFLIP
					checkpoint->movetype = MOVETYPE_SPINFLIP;
				else {
					if(num_molecules_adiabatic && (get_rand() < 0.5))
						// DISPLACE
						checkpoint->movetype = MOVETYPE_ADIABATIC;	// for the adiabatic mole fraction
					else
						checkpoint->movetype = MOVETYPE_DISPLACE;
				}
	
			} else { 
				// DISPLACE
				if(num_molecules_adiabatic && (get_rand() < 0.5))
					// DISPLACE
					checkpoint->movetype = MOVETYPE_ADIABATIC;	// for the adiabatic mole fraction
				else
					checkpoint->movetype = MOVETYPE_DISPLACE;
			}
			break;


		case ENSEMBLE_NVT :
		case ENSEMBLE_NVE :

			if( quantum_rotation && (get_rand() < spinflip_probability))
				checkpoint->movetype = MOVETYPE_SPINFLIP;
			else
				checkpoint->movetype = MOVETYPE_DISPLACE;
			break;


		case ENSEMBLE_NPT :	
	
			if( volume_probability == 0.0 ) { //if volume probability isn't set, then do volume moves with prob = 1/nmolecules
				if ( get_rand() < 1.0 / observables->N ) 
					checkpoint->movetype = MOVETYPE_VOLUME;
				else 
					checkpoint->movetype = MOVETYPE_DISPLACE;
			}
			else { //if volume probability IS set
				if ( get_rand() < volume_probability )
					checkpoint->movetype = MOVETYPE_VOLUME;
				else
					checkpoint->movetype = MOVETYPE_DISPLACE;
			}
			break;


		case ENSEMBLE_NVT_GIBBS:
			// see System.MonteCarlo.cpp:
			// System::backup_observables(systems);
			// System::pick_Gibbs_move(systems);
		case ENSEMBLE_PATH_INTEGRAL_NVT:
			// see ...
			
		default:
			Output::err("CHECKPOINT: invalid ensemble\n");
			throw invalid_ensemble;
		
	}
	
	
	// Randomly pick a (moveable) atom to remove, alter or duplicate (in the case of an insert)
	///////////////////////////////////////////////////////////////////////////////////////////////
	{

		if( checkpoint->movetype == MOVETYPE_ADIABATIC ) {
			// if we have any adiabatic molecules, then we have to do some special stuff
			--num_molecules_adiabatic;
			num_molecules_adiabatic_double = (double)num_molecules_adiabatic;
			altered = num_molecules_adiabatic - (int)rint(  num_molecules_adiabatic_double * get_rand()  );
			checkpoint->molecule_altered = ptr_array_adiabatic[altered];

		} else {
			
			if(   num_insertion_molecules   &&   checkpoint->movetype == MOVETYPE_INSERT   )
			{
				// If the move is an insertion and a molecule insertion list
				// has been set up, then select a molecule from the insertion list
				alt = (int)floor(   get_rand() * (double)num_insertion_molecules   );
				checkpoint->molecule_altered = insertion_molecules_array[ alt ];
				//needed to calculate boltzmann factor
				sorbateInsert = alt;

			}	// multi-sorbate + insert
			else
			{
				// Otherwise, pick one of the system molecules to alter
				--num_molecules_exchange;
				altered = (int)floor(  get_rand() * observables->N  );
				checkpoint->molecule_altered = ptr_array_exchange[altered];

				// if removing a molecule in a multi sorbate system, we also need to record the type
				if(  num_insertion_molecules   &&   checkpoint->movetype == MOVETYPE_REMOVE  ) {
					alt = 0;
					for( int j=0; j<sorbateCount; j++ ) {
						if ( SafeOps::iequals( sorbateInfo[j].id, ptr_array_exchange[altered]->moleculetype ) ) {
							sorbateInsert = alt;
							break;
						}
						alt++;
					}
				}	// multi-sorbate remove
			} // MOVETYPE_DISPLACE / REMOVE / SPINFLIP / VOLUME
		} //end non-adiabatic


		// MPMC must have at least one molecule in the system. If the move is to remove
		// the last molecule, then change this move to a SPINFLIP or DISPLACE
		///////////////////////////////////////////////////////////////////////////////////////////
		if(  !num_molecules_exchange  &&  checkpoint->movetype == MOVETYPE_REMOVE  ) {
			if( quantum_rotation && (get_rand() < spinflip_probability) )
				checkpoint->movetype = MOVETYPE_SPINFLIP;
			else
				checkpoint->movetype = MOVETYPE_DISPLACE;
		}


		// Determine the head and tail of the selected molecule                              
		///////////////////////////////////////////////////////////////////////////////////////////
		if( num_insertion_molecules  &&  checkpoint->movetype == MOVETYPE_INSERT ) {

			// When using an insertion list, we will always insert at the end of the list.
		
			// Advance to the end of the list
			molecule_ptr = molecules;
			while( molecule_ptr->next )
				molecule_ptr = molecule_ptr->next; 

			// Set head to the last molecule in the list, tail to the list terminator.
			checkpoint->head = molecule_ptr;
			checkpoint->tail = nullptr;

		} else {
			// If this is not a insertion using an insertion list, perform the original MPMC treatment:
			prev_molecule_ptr = nullptr;
			for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
				if(molecule_ptr == checkpoint->molecule_altered ) {
					checkpoint->head = prev_molecule_ptr;
					checkpoint->tail = molecule_ptr->next;
					break;
				}
				prev_molecule_ptr = molecule_ptr;
			}
		}
		

	} // End "pick a molecule to alter"


	// Free memory and make a backup of the molecule that will be altered
	///////////////////////////////////////////////////////////////////////////////////////////////
	{
		// free our temporary eligibility lists
		if( ptr_array_exchange )
			free(ptr_array_exchange);
		if(ptr_array_adiabatic)
			free(ptr_array_adiabatic);

		// if we have a molecule already backed up (from a previous accept), go ahead and free it
		if( checkpoint->molecule_backup )
			delete checkpoint->molecule_backup;
		// backup the state that will be altered
		checkpoint->molecule_backup = new Molecule( *checkpoint->molecule_altered );
	}
}




void System::backup_observables( std::vector<System*> &sys ) {
// Backup the observables for systems (intended for use in Gibbs ensemble)
	int nSystems = (int) sys.size();
	for (int i = 0; i < nSystems; i++)
		std::memcpy( sys[i]->checkpoint->observables, sys[i]->observables, sizeof(System::observables_t));
}




int System::pick_Gibbs_move( std::vector<System*> &sys ) {
// this function determines what move will be made next time make_move() is called and selects
// the molecule to which said move will be applied (and creates pointers to its list location).

	int          i_exchange                        = 0,
	             i_adiabatic                       = 0,
	             num_molecules_exchange [2]        = { 0 },
	             num_molecules_adiabatic[2]        = { 0 },
	             alt                               = 0,
	             altered                           = 0;
	double       num_molecules_adiabatic_double[2] = { 0 };
	Molecule  ** ptr_array_exchange[2]             = { nullptr },
	          ** ptr_array_adiabatic[2]            = { nullptr },
	           * molecule_ptr                      = { nullptr },
	           * prev_molecule_ptr                 = { nullptr };


	//|  Count exchangeable and adiabatic molecules, then allocate an array whose size is 
	//|  determined by the count. Populate the array with pointers to the molecules that
	//|  are eligible to be displaced, removed or altered. 
	//|_______________________________________________________________________________________________________

	// allocate array mem
	for (int i = 0; i < 2; i++) {

		num_molecules_exchange [i] = 0;
		num_molecules_adiabatic[i] = 0;

		// count molecule types
		for (molecule_ptr = sys[i]->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
			if (!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target))
				++num_molecules_exchange[i];
			if (molecule_ptr->adiabatic)
				++num_molecules_adiabatic[i];
		}

		if (num_molecules_exchange[i] > 0)
			SafeOps::calloc(ptr_array_exchange[i], num_molecules_exchange[i], sizeof(Molecule *), __LINE__, __FILE__);
		if (num_molecules_adiabatic[i] > 0)
			SafeOps::calloc(ptr_array_adiabatic[i], num_molecules_adiabatic[i], sizeof(Molecule *), __LINE__, __FILE__);

		// populate arrays with molecule pointers
		i_exchange  = 0;
		i_adiabatic = 0;
		for (molecule_ptr = sys[i]->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
			if (!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target)) {
				ptr_array_exchange[i][i_exchange] = molecule_ptr;
				++i_exchange;
			}
			if (molecule_ptr->adiabatic) {
				ptr_array_adiabatic[i][i_adiabatic] = molecule_ptr;
				++i_adiabatic;
			}
		}
	}  // for (create arrays)
	
	
	//|  Determine what kind of move to do 
	//|_______________________________________________________________________________________________________
	{
		double spinflip_prob = sys[0]->spinflip_probability;
		double volume_prob   = sys[0]->volume_probability + spinflip_prob;
		double transfer_prob = sys[0]->transfer_probability + volume_prob;
		double dice_roll     = sys[0]->get_rand();
		
		if( sys[0]->quantum_rotation && (dice_roll < spinflip_prob)) {
			sys[0]->checkpoint->movetype = MOVETYPE_SPINFLIP;
			sys[1]->checkpoint->movetype = MOVETYPE_SPINFLIP;
		} else {
			if ( dice_roll < volume_prob) {
				// volume change move
				sys[0]->checkpoint->movetype = MOVETYPE_VOLUME;
				sys[1]->checkpoint->movetype = MOVETYPE_VOLUME;
			} else {
				if( dice_roll < transfer_prob) {
					if (sys[0]->get_rand() < 0.5) {
						// transfer particle from system 0 to system 1
						sys[0]->checkpoint->movetype = MOVETYPE_REMOVE;
						sys[1]->checkpoint->movetype = MOVETYPE_INSERT;
					} else {
						// transfer particle from system 1 to system 0
						sys[0]->checkpoint->movetype = MOVETYPE_INSERT;
						sys[1]->checkpoint->movetype = MOVETYPE_REMOVE;
					}
				} else {
					// Default is to do a displacement move
					sys[0]->checkpoint->movetype = MOVETYPE_DISPLACE;
					sys[1]->checkpoint->movetype = MOVETYPE_DISPLACE;
				}
			}
		}
	}
	
		
	//|  Randomly pick a (moveable) atom to remove, alter or duplicate (in the case of an insert)
	//|_______________________________________________________________________________________________________

	for(int i=0; i<2; i++ ){

		if (sys[i]->checkpoint->movetype == MOVETYPE_ADIABATIC) {
			// if we have any adiabatic molecules, then we have to do some special stuff
			--num_molecules_adiabatic[i];
			num_molecules_adiabatic_double[i] = (double)num_molecules_adiabatic[i];
			altered = num_molecules_adiabatic[i] - (int)rint(num_molecules_adiabatic_double[i] * sys[i]->get_rand());
			sys[i]->checkpoint->molecule_altered = ptr_array_adiabatic[i][altered];

		}
		else {
			
			if (sys[i]->num_insertion_molecules   &&   sys[i]->checkpoint->movetype == MOVETYPE_INSERT)
			{
				// multi-sorbate + insert
				// If the move is an insertion and a molecule insertion list
				// has been set up, then select a molecule from the insertion list
				alt = (int)floor(sys[i]->get_rand() * (double)sys[i]->num_insertion_molecules);
				sys[i]->checkpoint->molecule_altered = sys[i]->insertion_molecules_array[alt];
				//needed to calculate boltzmann factor
				sys[i]->sorbateInsert = alt;

			}	
			else
			{
				// Otherwise, pick one of the system molecules to alter
				--num_molecules_exchange[i];
				altered = (int) floor( sys[i]->get_rand() * sys[i]->observables->N );
				sys[i]->checkpoint->molecule_altered = ptr_array_exchange[i][altered];

				// if removing a molecule in a multi sorbate system, we also need to record the type
				if (sys[i]->num_insertion_molecules   &&   sys[i]->checkpoint->movetype == MOVETYPE_REMOVE) {
					alt = 0;
					for (int j = 0; j<sys[i]->sorbateCount; j++) {
						if (SafeOps::iequals(sys[i]->sorbateInfo[j].id, ptr_array_exchange[i][altered]->moleculetype)) {
							sys[i]->sorbateInsert = alt;
							break;
						}
						alt++;
					}
				}	// multi-sorbate remove
			} // MOVETYPE_DISPLACE / REMOVE / SPINFLIP / VOLUME
		} //end non-adiabatic



		// MPMC must have at least one molecule in the system. If the move is to remove
		// the last molecule, then change this move to a SPINFLIP or DISPLACE
		///////////////////////////////////////////////////////////////////////////////////////
		if ((!num_molecules_exchange[0] && sys[0]->checkpoint->movetype == MOVETYPE_REMOVE) || (!num_molecules_exchange[1] && sys[1]->checkpoint->movetype == MOVETYPE_REMOVE)) {
			for (int l = 0; l<2; l++)
				if (sys[l]->quantum_rotation && (sys[l]->get_rand() < sys[l]->spinflip_probability))
					sys[l]->checkpoint->movetype = MOVETYPE_SPINFLIP;
				else
					sys[l]->checkpoint->movetype = MOVETYPE_DISPLACE;
		}



		// Determine the head and tail of the selected molecule                              
		//////////////////////////////////////////////////////////////////////////////////////
		if( sys[i]->num_insertion_molecules  &&  sys[i]->checkpoint->movetype == MOVETYPE_INSERT) {

			// When using an insertion list, we will always insert at the end of the list.

			// Advance to the end of the list
			molecule_ptr = sys[i]->molecules;
			while (molecule_ptr->next)
				molecule_ptr = molecule_ptr->next;

			// Set head to the last molecule in the list, tail to the list terminator.
			sys[i]->checkpoint->head = molecule_ptr;
			sys[i]->checkpoint->tail = nullptr;

		}
		else {
			// If this is not a insertion using an insertion list, perform the original MPMC treatment:
			prev_molecule_ptr = nullptr;
			for (molecule_ptr = sys[i]->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
				if (molecule_ptr == sys[i]->checkpoint->molecule_altered) {
					sys[i]->checkpoint->head = prev_molecule_ptr;
					sys[i]->checkpoint->tail = molecule_ptr->next;
				}
				prev_molecule_ptr = molecule_ptr;
			}
		}
	} // End "pick a molecule to alter"


	// Clean up
	/////////////////////////////////////////////////////////////////////////////////////////////
	{
		for (int i = 0; i < 2; i++) {
			// free our temporary eligibility lists
			if (ptr_array_exchange[i])
				free(ptr_array_exchange[i]);
			if (ptr_array_adiabatic[i])
				free(ptr_array_adiabatic[i]);

			// if we have a molecule already backed up (from a previous accept), go ahead and free it
			if( sys[i]->checkpoint->molecule_backup)
				delete sys[i]->checkpoint->molecule_backup;
			// backup the state that will be altered
			sys[i]->checkpoint->molecule_backup = new Molecule(*sys[i]->checkpoint->molecule_altered);
		}
	}

	return sys[0]->checkpoint->movetype;
}




void System::make_move() {
// apply the move that was selected in the checkpoint

	int         cavities_array_counter =  0,
	            random_index           =  0;
	double      com [3]                = {0},
	            rand[3]                = {0};
	cavity_t  * cavities_array         = nullptr;
	Molecule  * molecule_ptr           = nullptr;
	Atom      * atom_ptr               = nullptr;
	Pair      * pair_ptr               = nullptr;

	// update the cavity grid prior to making a move 
	if( cavity_bias ) {
		cavity_update_grid();
		checkpoint->biased_move = 0;
	}
	
	switch( checkpoint->movetype ) {


		case MOVETYPE_INSERT :  // insert a molecule at a random pos and orientation 
			// umbrella sampling 
			if( cavity_bias && cavities_open ) {
				// doing a biased move - this flag lets mc.c know about it 
				checkpoint->biased_move = 1;
				// make an array of possible insertion points
				SafeOps::calloc( cavities_array, cavities_open, sizeof(cavity_t), __LINE__, __FILE__ );
				for( int i = 0; i < cavity_grid_size; i++) {
					for( int j = 0; j < cavity_grid_size; j++) {
						for(int k = 0; k < cavity_grid_size; k++) {
							if( ! cavity_grid[i][j][k].occupancy ) {
								for( int p = 0; p < 3; p++ )
									cavities_array[cavities_array_counter].pos[p] = cavity_grid[i][j][k].pos[p];
								++cavities_array_counter;
							}
						} // end k 
					} // end j 
				} // end i 
				// insert randomly at one of the free cavity points 
				random_index = (cavities_open - 1) - (int)rint(((double)(cavities_open - 1))*get_rand());
				for( int p = 0; p < 3; p++ )
					com[p] = cavities_array[random_index].pos[p];
				// free the insertion array 
				free(cavities_array);
			} // end umbrella

			else {
				// insert the molecule to a random location within the unit cell
				for( int p = 0; p < 3; p++ )
					rand[p] = 0.5 - get_rand();
				for( int p = 0; p < 3; p++ ) {
					com[p] = 0;
					for( int q = 0; q < 3; q++ )
						com[p] += pbc.basis[q][p]*rand[q];
				}
			}

			// process the inserted molecule 
			for( atom_ptr = checkpoint->molecule_backup->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
				// move the molecule back to the origin and then assign it to com 
				for( int p = 0; p < 3; p++ )
					atom_ptr->pos[p] += com[p] - checkpoint->molecule_backup->com[p];
			}

			// update the molecular com 
			for( int p = 0; p < 3; p++ )
				checkpoint->molecule_backup->com[p] = com[p];
			// give it a random orientation 
			checkpoint->molecule_backup->rotate_rand(1.0); // , pbc, &mt_rand );

			// insert into the list 
			if( num_insertion_molecules ) {
				// IS THIS RIGHT? Looks like an insertion at the point where the old molecule
				// was, followed by a truncation of the list at that point. 
				// ====== original comment follows ======
				// If inserting a molecule from an insertion list, we will always insert at the end
				checkpoint->head->next = checkpoint->molecule_backup;
				checkpoint->molecule_backup->next = nullptr;
			} else {
				if( ! checkpoint->head ) {      // if we're at the start of the list:
					molecules = checkpoint->molecule_backup;		
				} else {
					checkpoint->head->next = checkpoint->molecule_backup;
				}
				checkpoint->molecule_backup->next = checkpoint->molecule_altered;
			}

			// set new altered and tail to reflect the insertion 
			checkpoint->molecule_altered = checkpoint->molecule_backup;
			checkpoint->tail             = checkpoint->molecule_altered->next;
			checkpoint->molecule_backup  = nullptr;

			if( num_insertion_molecules ) { //multi sorbate
				// Free all pair memory in the list
				for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
					for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
						pair_ptr = atom_ptr->pairs;
						while( pair_ptr ) {
							Pair *temp = pair_ptr;
							pair_ptr = pair_ptr->next;
							free( temp );
						}
					}
				}
				// Generate new pairs lists for all atoms in system
				allocate_pair_lists();
			} // only one sorbate
			else		
				update_pairs_insert();

			//reset atom and molecule id's
			enumerate_particles();

			break;


		case MOVETYPE_REMOVE : // remove a randomly chosen molecule 
	
			if( cavity_bias ) {
				if(get_rand() < pow((1.0 - avg_observables->cavity_bias_probability), ((double)cavity_grid_size*cavity_grid_size*cavity_grid_size)))
					checkpoint->biased_move = 0;
				else
					checkpoint->biased_move = 1;
			}
	
			// remove 'altered' from the list 
			if( ! checkpoint->head ) {	// handle the case where we're removing from the start of the list 
				checkpoint->molecule_altered = molecules;
				molecules = molecules->next;
			} else {
				checkpoint->head->next = checkpoint->tail;
			}
			//free_molecule( system, system->checkpoint->molecule_altered );
			delete checkpoint->molecule_altered;
			update_pairs_remove();

			//reset atom and molecule id's
			enumerate_particles();

			break;

		case MOVETYPE_DISPLACE :	
	
			// change coords of 'altered' 
			if( rd_anharmonic )
				displace_1D(checkpoint->molecule_altered, move_factor);
			else if( spectre )
				spectre_displace( checkpoint->molecule_altered, move_factor, spectre_max_charge, spectre_max_target );
			else if( gwp ) {
				if( checkpoint->molecule_altered->atoms->gwp_spin ) {
					displace( checkpoint->molecule_altered, pbc, gwp_probability, rot_factor);
					checkpoint->molecule_altered->displace_gwp( gwp_probability, &mt_rand );
				} else
					displace( checkpoint->molecule_altered, pbc, move_factor, rot_factor );
			} else
				displace( checkpoint->molecule_altered, pbc, move_factor, rot_factor);
		break;

		case MOVETYPE_ADIABATIC :
			// change coords of 'altered' 
			displace( checkpoint->molecule_altered, pbc, adiabatic_probability, 1.0 );
	
		break;
		case MOVETYPE_SPINFLIP :
	
			// XXX - should have separate function do spinfip move 
			if( checkpoint->molecule_altered->nuclear_spin == NUCLEAR_SPIN_PARA )
				checkpoint->molecule_altered->nuclear_spin = NUCLEAR_SPIN_ORTHO;
			else
				checkpoint->molecule_altered->nuclear_spin = NUCLEAR_SPIN_PARA;
	
		break;
		case MOVETYPE_VOLUME :
			volume_change(); // I don't want to contribute to the god damned mess -- kmclaugh
		break;

		default:
			Output::err("MC_MOVES: invalid mc move\n");
			throw invalid_monte_carlo_move;
	}
}




void System::make_move_Gibbs(std::vector<System*> &sys) {
	
	// update the cavity grids prior to making a move 
	for( int i=0; i<2; i++ )
		if( sys[i]->cavity_bias ) {
			sys[i]->cavity_update_grid();
			sys[i]->checkpoint->biased_move = 0;
		}



	switch( sys[0]->checkpoint->movetype ) {

		case MOVETYPE_DISPLACE:
		
			for( int i=0; i<2; i++ )
				// change coords of 'altered' 
				if( sys[i]->rd_anharmonic )
					sys[i]->displace_1D( sys[i]->checkpoint->molecule_altered, sys[i]->move_factor );
				else if( sys[i]->spectre )
					sys[i]->spectre_displace( sys[i]->checkpoint->molecule_altered, sys[i]->move_factor, sys[i]->spectre_max_charge, sys[i]->spectre_max_target );
				else if( sys[i]->gwp ) {
					if( sys[i]->checkpoint->molecule_altered->atoms->gwp_spin ) {
						sys[i]->displace( sys[i]->checkpoint->molecule_altered, sys[i]->pbc, sys[i]->gwp_probability, sys[i]->rot_factor);
						sys[i]->checkpoint->molecule_altered->displace_gwp( sys[i]->gwp_probability, &sys[i]->mt_rand );
					}
					else
						sys[i]->displace( sys[i]->checkpoint->molecule_altered, sys[i]->pbc, sys[i]->move_factor, sys[i]->rot_factor );
				}
				else
					sys[i]->displace( sys[i]->checkpoint->molecule_altered, sys[i]->pbc, sys[i]->move_factor, sys[i]->rot_factor );
			break;
			

		case MOVETYPE_SPINFLIP:
			for (int i = 0; i < 2; i++) 
				if( sys[i]->checkpoint->molecule_altered->nuclear_spin == NUCLEAR_SPIN_PARA )
					sys[i]->checkpoint->molecule_altered->nuclear_spin = NUCLEAR_SPIN_ORTHO;
				else
					sys[i]->checkpoint->molecule_altered->nuclear_spin = NUCLEAR_SPIN_PARA;
			break;


		case MOVETYPE_VOLUME:
			volume_change_Gibbs( sys );
			break;

		// Gibbs Transfer Move
		case MOVETYPE_INSERT:
		case MOVETYPE_REMOVE:
			

			for (int s = 0; s < 2; s++)
			{
				if (sys[s]->checkpoint->movetype == MOVETYPE_INSERT)
				{

					cavity_t  * cavities_array = nullptr;
					Molecule  * molecule_ptr   = nullptr;
					Atom      * atom_ptr       = nullptr;
					Pair      * pair_ptr       = nullptr;

					int         random_index           = 0;
					double      com[3]                 = { 0 },
					            rand[3]                = { 0 };
					




					// SYSTEM s
					////////////////////////////////////////////////////////////////////////////////////////////////////////
					// insert a molecule in at a random pos and orientation 
					// umbrella sampling 
					if (sys[s]->cavity_bias && sys[s]->cavities_open) {
						// doing a biased move - this flag lets mc.c know about it 
						sys[s]->checkpoint->biased_move = 1;
						// make an array of possible insertion points
						SafeOps::calloc(cavities_array, sys[s]->cavities_open, sizeof(cavity_t), __LINE__, __FILE__);
						int cavities_array_counter = 0; // keeps track of position in cavities_array
						for (int i = 0; i < sys[s]->cavity_grid_size; i++) {
							for (int j = 0; j < sys[s]->cavity_grid_size; j++) {
								for (int k = 0; k < sys[s]->cavity_grid_size; k++) {
									if (!sys[s]->cavity_grid[i][j][k].occupancy) {
										for (int p = 0; p < 3; p++)
											cavities_array[cavities_array_counter].pos[p] = sys[s]->cavity_grid[i][j][k].pos[p];
										++cavities_array_counter;
									}
								} // end k 
							} // end j 
						} // end i 

						// insert randomly at one of the free cavity points 
						random_index = (sys[s]->cavities_open - 1) - (int)rint(((double)(sys[s]->cavities_open - 1)) * sys[s]->get_rand());
						for (int p = 0; p < 3; p++)
							com[p] = cavities_array[random_index].pos[p];
						// free the insertion array 
						free(cavities_array);
					} // end umbrella

					else {
						// insert the molecule to a random location within the unit cell
						for (int p = 0; p < 3; p++)
							rand[p] = 0.5 - sys[s]->get_rand();
						for (int p = 0; p < 3; p++) {
							com[p] = 0;
							for (int q = 0; q < 3; q++)
								com[p] += sys[s]->pbc.basis[q][p] * rand[q];
						}
					}

					// process the inserted molecule 
					for (atom_ptr = sys[s]->checkpoint->molecule_backup->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
						// move the molecule back to the origin and then assign it to com 
						for (int p = 0; p < 3; p++)
							atom_ptr->pos[p] += com[p] - sys[s]->checkpoint->molecule_backup->com[p];
					}

					// update the molecular com 
					for (int p = 0; p < 3; p++)
						sys[s]->checkpoint->molecule_backup->com[p] = com[p];
					// give it a random orientation 
					sys[s]->checkpoint->molecule_backup->rotate_rand(1.0); //  , sys[s]->pbc, &sys[s]->mt_rand);

					// insert into the list 
					if (sys[s]->num_insertion_molecules) {

						// If inserting a molecule from an insertion list, we will always insert at the end
						sys[s]->checkpoint->head->next = sys[s]->checkpoint->molecule_backup;
						sys[s]->checkpoint->molecule_backup->next = nullptr;

					}
					else {

						if (!sys[s]->checkpoint->head) {      // if we're at the start of the list:
							sys[s]->molecules = sys[s]->checkpoint->molecule_backup;
						}
						else {
							sys[s]->checkpoint->head->next = sys[s]->checkpoint->molecule_backup;
						}
						sys[s]->checkpoint->molecule_backup->next = sys[s]->checkpoint->molecule_altered;
					}

					// set new altered and tail to reflect the insertion 
					sys[s]->checkpoint->molecule_altered = sys[s]->checkpoint->molecule_backup;
					sys[s]->checkpoint->tail = sys[s]->checkpoint->molecule_altered->next;
					sys[s]->checkpoint->molecule_backup = nullptr;

					if (sys[s]->num_insertion_molecules) { //multi sorbate
						// Free all pair memory in the list
						for (molecule_ptr = sys[s]->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
							for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
								pair_ptr = atom_ptr->pairs;
								while (pair_ptr) {
									Pair *temp = pair_ptr;
									pair_ptr = pair_ptr->next;
									free(temp);
								}
							}
						}
						// Generate new pairs lists for all atoms in system
						sys[s]->allocate_pair_lists();
					} // only one sorbate
					else
						sys[s]->update_pairs_insert();

					//reset atom and molecule id's
					sys[s]->enumerate_particles();




					// SYSTEM !s  (The OTHER system)
					////////////////////////////////////////////////////////////////////////////////////////////////////////
					// remove a randomly chosen molecule 

					if (sys[!s]->cavity_bias) {
						if (sys[!s]->get_rand() < pow((1.0 - sys[!s]->avg_observables->cavity_bias_probability), ((double)sys[!s]->cavity_grid_size * sys[!s]->cavity_grid_size* sys[!s]->cavity_grid_size)))
							sys[!s]->checkpoint->biased_move = 0;
						else
							sys[!s]->checkpoint->biased_move = 1;
					}

					// remove 'altered' from the list 
					if (!sys[!s]->checkpoint->head) {	// handle the case where we're removing from the start of the list 
						sys[!s]->checkpoint->molecule_altered = sys[!s]->molecules;
						sys[!s]->molecules = sys[!s]->molecules->next;
					}
					else {
						sys[!s]->checkpoint->head->next = sys[!s]->checkpoint->tail;
					}
					//free_molecule( system, system->checkpoint->molecule_altered );
					delete sys[!s]->checkpoint->molecule_altered;
					sys[!s]->update_pairs_remove();

					//reset atom and molecule id's
					sys[!s]->enumerate_particles();
				}
			}
			break;



		default:
			Output::err("MC_MOVES: invalid mc move\n");
			throw invalid_monte_carlo_move;
	}
}




void System::enumerate_particles () {
// re-enumerate atoms and molecules -> set atom and molecule id's which get messed up in UVT runs
	Molecule  * mptr;
	Atom      * aptr;
	int         aid = 1,
	            mid = 1;
	
	for( mptr = molecules; mptr; mptr=mptr->next ) {
		mptr->id = mid++;
		for( aptr = mptr->atoms; aptr; aptr=aptr->next )
			aptr->id = aid++;
	}
}




void System::displace_1D( Molecule *molecule, double scale ) {
// perform a 1D translation without periodic boundaries

	Atom   * atom_ptr = nullptr;
	double   trans    = 0;

	trans = scale * get_rand();
	if(get_rand() < 0.5) trans *= -1.0;
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next)
		atom_ptr->pos[0] += trans;

	molecule->com[0] += trans;

}	




void System::spectre_displace( Molecule *molecule, double trans_scale, double max_charge, double max_target ) {
// perform a random translation/rotation of a molecule for SPECTRE algorithm
	Atom  * atom_ptr     = nullptr;
	int     p            =  0;
	double  delta_charge =  0,
	        trans[3]     = {0};
	

	// randomly translate 
	for(p = 0; p < 3; p++) {
		trans[p] = trans_scale * get_rand() * max_target; 
		if(get_rand() < 0.5) 
			trans[p] *= -1.0;
	}
	
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

		// position reassignment
		for(p = 0; p < 3; p++)
			atom_ptr->pos[p] += trans[p];

		// charge reassignment
		do {
			delta_charge = get_rand(); if(get_rand() < 0.5) delta_charge *= -1.0;
		} while (fabs(atom_ptr->charge + delta_charge) > max_charge);
		atom_ptr->charge += delta_charge;


	}

	// restrict the SPECTRE charges to the domain
	spectre_wrapall();

	// renormalize all SPECTRE charges
	spectre_charge_renormalize();

}




void System::spectre_charge_renormalize() {
// ensure neutrality for the SPECTRE system 

	Molecule * molecule_ptr    = nullptr;
	Atom     * atom_ptr        = nullptr;
	int        num_spectre     = 0;
	double     residual_charge = 0,
	           frac_charge     = 0;

	// calculate the net charge
	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			if(atom_ptr->spectre) {
				++num_spectre;
				residual_charge += atom_ptr->charge;
			}
		}
	}

	// spread that charge out among the other SPECTRE sites
	frac_charge = -1.0*residual_charge/((double)num_spectre);

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next )
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next )
			if(atom_ptr->spectre)
				atom_ptr->charge += frac_charge;

	return;
}




void System::displace( Molecule *molecule, const PeriodicBoundary &PBC, double trans_scale, double rot_scale ) {
// perform a random translation/rotation of a molecule		
	molecule->translate_rand_pbc( trans_scale, PBC, &mt_rand );
	molecule->rotate_rand( rot_scale ); 
}




void System::volume_change() {
// scale the volume : used in NPT ensemble	
	Molecule * m;
	Atom     * a;
	double     new_volume         =  0,
	           log_new_volume     =  0,
	           basis_scale_factor =  0;
	double     new_com[3]         = {0},  
	           old_com[3]         = {0}, 
	           delta_pos[3]       = {0};
	//int i,j;

	// figure out what the new volume will be
	if( ensemble == ENSEMBLE_REPLAY )
		//if ensemble replay, then we're just trying to calculate the pressure via dV change
		new_volume = pbc.volume + calc_pressure_dv;
	else {
		log_new_volume=log( pbc.volume )   +   (get_rand()-0.5) * volume_change_factor;
		new_volume=exp( log_new_volume );
	}

	//scale basis
	basis_scale_factor = pow(new_volume/pbc.volume,1.0/3.0);
	for ( int i=0; i<3; i++ )
		for ( int j=0; j<3; j++ )
			pbc.basis[i][j] *= basis_scale_factor;

	//recalculate PBC stuff (volume/reciprocal basis/cutoff)
	update_pbc();
	observables->volume = pbc.volume;

	//scale molecule positions
	for ( m = molecules; m; m = m->next ) {
		for ( int i=0; i<3; i++ ) {
			old_com[i]   = m-> com[i]; //molecule ptr's com will be udated during pairs routine
			new_com[i]   = m-> com[i] * basis_scale_factor;
			delta_pos[i] = new_com[i] - old_com[i];
		}
		for ( a = m->atoms; a; a=a->next ) { //calculate new atomic positions based on new com
			for ( int i=0; i<3; i++ ) {
				a->pos[i]         += delta_pos[i];
				a->wrapped_pos[i] += delta_pos[i];
			}
		}
	}		

	return;
}




void System::volume_change_Gibbs(std::vector<System*> &sys) {

	Molecule * m;
	Atom     * a;
	double     new_volume[2] = { 0 },
	           log_new_volume = 0;
	

	// figure out what the new volume will be
	if (sys[0]->ensemble == ENSEMBLE_REPLAY) {
		//if ensemble replay, then we're just trying to calculate the pressure via dV change
		new_volume[0] = sys[0]->pbc.volume + sys[0]->calc_pressure_dv;
		new_volume[1] = sys[1]->pbc.volume + sys[1]->calc_pressure_dv;
	}
	else do {
		log_new_volume = log(sys[0]->pbc.volume) + (sys[0]->get_rand() - 0.5) * sys[0]->volume_change_factor;
		new_volume[0] = exp(log_new_volume);
		new_volume[1] = sys[1]->pbc.volume + sys[0]->pbc.volume - new_volume[0];
	} while( new_volume[1] <= 0.0 );


	for (int s = 0; s < 2; s++) {

		//scale basis
		double basis_scale_factor = pow(new_volume[s] / sys[s]->pbc.volume, 1.0 / 3.0);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				sys[s]->pbc.basis[i][j] *= basis_scale_factor;

		//recalculate PBC stuff (volume/reciprocal basis/cutoff)
		sys[s]->update_pbc();
		sys[s]->observables->volume = sys[s]->pbc.volume;


		double new_com  [3] = { 0 },
		       old_com  [3] = { 0 },
		       delta_pos[3] = { 0 };
		
		//scale molecule positions
		for (m = sys[s]->molecules; m; m = m->next) {
			for (int i = 0; i < 3; i++) {
				old_com[i] = m->com[i]; //molecule ptr's com will be udated during pairs routine
				new_com[i] = m->com[i] * basis_scale_factor;
				delta_pos[i] = new_com[i] - old_com[i];
			}
			for (a = m->atoms; a; a = a->next) { //calculate new atomic positions based on new com
				for (int i = 0; i < 3; i++) {
					a->pos[i] += delta_pos[i];
					a->wrapped_pos[i] += delta_pos[i];
				}
			}
		}
	}
}




void System::boltzmann_factor( double initial_energy, double final_energy) {
// the prime quantity of interest 
	double delta_energy   = final_energy - initial_energy,
	       u              = 0,
	       g              = 0, 
	       partfunc_ratio = 0,
	       v_new          = 0,
	       v_old          = 0,
	       fugacity       = 0;

	
	switch ( ensemble ) {

		case ENSEMBLE_UVT :
			//obtain the correct fugacity value
			if( h2_fugacity || co2_fugacity || ch4_fugacity || n2_fugacity )
				fugacity = fugacities[0];
			else if ( user_fugacities )
				fugacity = fugacities[sorbateInsert];
			else
				fugacity = pressure;

			// if biased_move not set, no cavity available so do normal evaluation below
			if( cavity_bias && checkpoint->biased_move ) {
				// modified metropolis function
				switch ( checkpoint->movetype ) {
					case MOVETYPE_INSERT :
						nodestats->boltzmann_factor = 
							(cavity_volume * avg_nodestats->cavity_bias_probability * fugacity * ATM2REDUCED / 
							(temperature * observables->N )) * exp(-delta_energy/temperature) * (double) sorbateCount;
						break;
					case MOVETYPE_REMOVE :
						nodestats->boltzmann_factor = 
							(temperature * (observables->N + 1.0)) / 
							(cavity_volume * avg_nodestats->cavity_bias_probability * fugacity*ATM2REDUCED) *
							exp(-delta_energy/temperature) / (double) sorbateCount;
						break;
					case MOVETYPE_DISPLACE :
						nodestats->boltzmann_factor = exp(-delta_energy/temperature);
						break;
					default :
						Output::err("MC: invalid mc move (not implemented for biased moves?)\n");
						throw invalid_monte_carlo_move;
				} //mc move switch
			} //end biased move
			else {
				switch ( checkpoint->movetype ) {
					case MOVETYPE_INSERT :
						nodestats->boltzmann_factor =
							pbc.volume * fugacity * ATM2REDUCED/(temperature * (double)(observables->N)) *
							exp(-delta_energy/temperature) *
							(double)(sorbateCount); //for add/delete bias
					break;
					case MOVETYPE_REMOVE :
						nodestats->boltzmann_factor =
							temperature * ((double)(observables->N) + 1.0)/(pbc.volume*fugacity*ATM2REDUCED) *
							exp(-delta_energy/temperature) /
						  (double)(sorbateCount); //for add/delete bias
					break;
					case MOVETYPE_DISPLACE :
						nodestats->boltzmann_factor = exp(-delta_energy/temperature);
					break;
					case MOVETYPE_SPINFLIP :
						g = checkpoint->molecule_altered->rot_partfunc_g;
						u = checkpoint->molecule_altered->rot_partfunc_u;
						if( checkpoint->molecule_altered->nuclear_spin == NUCLEAR_SPIN_PARA ) 
							partfunc_ratio = g/(g+u);
						else
							partfunc_ratio = u/(g+u);
						// set the boltz factor, including ratio of partfuncs for different symmetry rotational levels
						nodestats->boltzmann_factor = partfunc_ratio;
						break;
					default:
						Output::err("MC: invalid mc move (not implemented for binary mixtures?)\n");
						throw invalid_monte_carlo_move;
				}//MC move switch
			} //end biased or not?
		break; //end UVT

		case ENSEMBLE_NVT :
			switch ( checkpoint->movetype ) {
				case MOVETYPE_SPINFLIP : 
					g = checkpoint->molecule_altered->rot_partfunc_g;
					u = checkpoint->molecule_altered->rot_partfunc_u;
					if ( checkpoint->molecule_altered->nuclear_spin == NUCLEAR_SPIN_PARA ) 
						partfunc_ratio = g/(g+u);
					else
						partfunc_ratio = u/(g+u);
					// set the boltz factor, including ratio of partfuncs for different symmetry rotational levels
					nodestats->boltzmann_factor = partfunc_ratio;
				break;
				default: // DISPLACE
					nodestats->boltzmann_factor = exp(-delta_energy/temperature);
			}
		break; //end NVT

		case ENSEMBLE_NPT : 
			//either displace or change volume
			switch ( checkpoint->movetype ) {
				case MOVETYPE_VOLUME : 
			
					v_old = checkpoint->observables->volume;
					v_new = observables->volume;
					nodestats->boltzmann_factor = 
						exp(-( delta_energy 
						    + pressure * ATM2REDUCED * (v_new-v_old) 
						    - (observables->N + 1) * temperature * log(v_new/v_old)
						)/temperature);
				break;
				default:	// displace
					nodestats->boltzmann_factor = exp(-delta_energy/temperature);
			}
		break;

		case ENSEMBLE_NVE :
			nodestats->boltzmann_factor  = pow((total_energy - final_energy  ), 3.0*N/2.0 );
			nodestats->boltzmann_factor /= pow((total_energy - initial_energy), 3.0*N/2.0 );
		break;

		default:
			Output::err("MC: invalid ensemble. aborting.\n");
			throw invalid_ensemble;
		}

	return;
}




void System::register_accept() {
// keep track of which specific moves were accepted
	++nodestats->accept;
	switch( checkpoint->movetype ) {

		case MOVETYPE_INSERT:
			++nodestats->accept_insert;
			break;
		case MOVETYPE_REMOVE:
			++nodestats->accept_remove;
			break;
		case MOVETYPE_DISPLACE:
			++nodestats->accept_displace;
			break;
		case MOVETYPE_ADIABATIC:
			++nodestats->accept_adiabatic;
			break;
		case MOVETYPE_SPINFLIP:
			++nodestats->accept_spinflip;
			break;
		case MOVETYPE_VOLUME:
			++nodestats->accept_volume;
			break;
		case MOVETYPE_PERTURB_BEADS: 
			++nodestats->accept_beadPerturb;
	}
}




void System::restore() {
// this function 
//    (a) undoes what make_move() did and...
//    (b) determines the next move sequence by calling do_checkpoint()
	
	// restore the remaining observables 
	std::memcpy( observables, checkpoint->observables, sizeof(observables_t) );

	// restore state by undoing the steps of make_move()
	switch ( checkpoint->movetype ) {

		case MOVETYPE_INSERT : 
	
			// take altered out of the list
			if( ! checkpoint->head ) { 
				// handle the case where we inserted at the head of the list
				molecules = molecules->next;
			} else {
				checkpoint->head->next = checkpoint->tail;
			}
			unupdate_pairs_insert();
			delete checkpoint->molecule_altered;
			// Previously, delete was: free_molecule(system, system->checkpoint->molecule_altered);

			//reset atom and molecule id's
			enumerate_particles();

		break;
		case MOVETYPE_REMOVE :
	
			// put backup back into the list
			if( ! checkpoint->head ) {
				molecules = checkpoint->molecule_backup;
			} else {
				checkpoint->head->next = checkpoint->molecule_backup;
			}
			checkpoint->molecule_backup->next = checkpoint->tail;
			unupdate_pairs_remove();
			checkpoint->molecule_backup = nullptr;

			//reset atom and molecule id's
			enumerate_particles();

		break;
		case MOVETYPE_VOLUME :
	
			revert_volume_change();

		break;
		case MOVETYPE_DISPLACE :
		case MOVETYPE_PERTURB_BEADS:
		default :

			// link the backup  molecule back into the working list
			// checkpoint->head will be nullptr if the first item in the system list was the backed up molecule
			if( ! checkpoint->head )
				molecules = checkpoint->molecule_backup;
			else
				checkpoint->head->next = checkpoint->molecule_backup;
			checkpoint->molecule_backup->next = checkpoint->tail;
			// Delete the molecule configuration that was rejected. The pair list pointer is zeroed in order to keep
			// the pair list intact as it was transferred to the accepted molecule and is needed. 
			checkpoint->molecule_altered->wipe_pair_refs();
			delete checkpoint->molecule_altered;
			checkpoint->molecule_altered = nullptr;
			// wipe "backup" reference, since this molecule is now linked into the system list
			checkpoint->molecule_backup  = nullptr;
	}	
	
	// renormalize charges
	if( spectre ) spectre_charge_renormalize();

	// eject early for the following ensembles:
	if((ENSEMBLE_PATH_INTEGRAL_NVT==ensemble)  ||  (ENSEMBLE_NVT_GIBBS==ensemble)) return;

	// establish the previous checkpoint again
	do_checkpoint();
}




void System::unupdate_pairs_insert() {
// if an insert move is rejected, remove the pairs that were previously added
	int        n            = 0,
	           m            = 0;
	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	Pair     * pair_ptr     = nullptr, 
	        ** pair_array   = nullptr; 


	// count the number of atoms per molecule
	for(atom_ptr = checkpoint->molecule_altered->atoms; atom_ptr; atom_ptr = atom_ptr->next)
		++n;

	// remove n number of pairs for all molecules ahead of the removal point
	SafeOps::calloc( pair_array, 1, sizeof(Pair *), __LINE__, __FILE__ );

	for( molecule_ptr = molecules;   molecule_ptr != checkpoint->tail;   molecule_ptr = molecule_ptr->next ) {
		for( atom_ptr = molecule_ptr->atoms;   atom_ptr;   atom_ptr = atom_ptr->next ) {

			// build the pair pointer array
			for( m=0, pair_ptr = atom_ptr->pairs;  pair_ptr;  pair_ptr = pair_ptr->next, m++ ) {

				SafeOps::realloc( pair_array, sizeof(Pair *)*(m + 1), __LINE__, __FILE__ );
				pair_array[m] = pair_ptr;
			}

			for( int i = (m - n); i < m; i++)
				free(pair_array[i]);

			// handle the end of the list
			if((m - n) > 0)
				pair_array[(m - n - 1)]->next = nullptr;
			else
				atom_ptr->pairs = nullptr;

		}
	}

	// free our temporary array
	free(pair_array);
}




void System::unupdate_pairs_remove() {
// if a remove is rejected and a molecule that was removed from the system has been added back
// to the molecule list, then add back the number of pair nodes that were previously removed.

	int         n            = 0;       // The number of atoms/sites in the molecule that is to be added.
	Molecule  * molecule_ptr = nullptr;
	Atom      * atom_ptr     = nullptr;
	Pair      * pair_ptr     = nullptr;


	// count the number of atoms per molecule
	for( atom_ptr = checkpoint->molecule_backup->atoms; atom_ptr; atom_ptr = atom_ptr->next)
		++n;

	// add n number of pairs to altered and all molecules ahead of it in the list
	for( molecule_ptr = molecules; molecule_ptr != checkpoint->molecule_backup; molecule_ptr = molecule_ptr->next ) {
		for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {

			// go to the end of the pair list
			if( atom_ptr->pairs ) {

				// go to the end of the pair list 
				for( pair_ptr = atom_ptr->pairs; pair_ptr->next; pair_ptr = pair_ptr->next);

				// tag on the extra pairs
				for( int i = 0; i < n; i++) {
					SafeOps::calloc( pair_ptr->next, 1, sizeof(Pair), __LINE__, __FILE__ );
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




void System::revert_volume_change() {
// revert changes to volume
	Molecule  * m                  = nullptr; 
	Atom      * a                  = nullptr;
	double      new_volume         = checkpoint->observables->volume,
	            basis_scale_factor = pow(new_volume/pbc.volume,1.0/3.0), // pbc.volume still has rejected value until update_pbc() is called
	            old_com[3]         = {0},
	            new_com[3]         = {0},
	            delta_pos[3]       = {0};
	
	// int i,j;
	
	for ( int i=0; i<3; i++ )
		for ( int j=0; j<3; j++ )
			pbc.basis[i][j] *= basis_scale_factor;

	//recalc PBC stuff
	update_pbc();
	observables->volume = pbc.volume;

	//scale molecule positions
	for ( m = molecules; m; m = m->next ) {
		for ( int i=0; i<3; i++ ) {
			old_com  [i] = m-> com[i]; //molecule ptr's com will ned to be updated since this is a revert (no energy call following)
			new_com  [i] = m-> com[i] * basis_scale_factor;
			m-> com  [i] = new_com[i];
			delta_pos[i] = new_com[i] - old_com[i];
		}
		for ( a = m->atoms; a; a=a->next ) { //calculate new atomic positions based on new com
			for( int i=0; i<3; i++ ) {
				a->pos        [i] += delta_pos[i];
				a->wrapped_pos[i] += delta_pos[i];
			}
		}
	}

	return;
}




void System::register_reject() {
// keep track of which specific moves were rejected 

	++nodestats->reject;
	switch( checkpoint->movetype ) {

		case MOVETYPE_INSERT:
			++nodestats->reject_insert;
			break;
		case MOVETYPE_REMOVE:
			++nodestats->reject_remove;
			break;
		case MOVETYPE_DISPLACE:
			++nodestats->reject_displace;
			break;
		case MOVETYPE_ADIABATIC:
			++nodestats->reject_adiabatic;
			break;
		case MOVETYPE_SPINFLIP:
			++nodestats->reject_spinflip;
			break;
		case MOVETYPE_VOLUME:
			++nodestats->reject_volume;
			break;
		case MOVETYPE_PERTURB_BEADS:
			++nodestats->reject_beadPerturb;

	}

}




void System::temper_system( double current_energy ) {
//  parallel tempering stuff
//  parallel tempering is treated differently than the other types of MC moves, because it is done
//  IN ADDITION to any other moves: i.e. some random move is performed and the energy is calculated
//  before and after. once that is all figured out, we then perform a temper step. The boltzmann factor
//  that is calculated will NOT be averaged into the BF quantity. It can be, but care must be taken so
//  there is no double-counting, and it really isn't so important to do so. 

	/*
		// system->ptemp->index[j] is a mapping from core -> bath_index (bath_index 0 is lowest temperature, etc.)
		// bath2core maps from bath_index -> core 
		// only half of the cores will be responsible for carrying out the calculations
		// partner_list maps from core -> swap partner, if core is a master. core -> -1, if core is a slave.
		int  * bath2core        = nullptr,
		     * partner_list     = nullptr,
		     * is_master        = nullptr,
		       master           = 0,
		       slave            = 0,
		     * update_index     = nullptr,
		       new_index_val    = 0;
		double slave_energy     = 0,
		       boltzmann_factor = 0,
		     * update_templist  = nullptr,
		       new_templist_val = 0;
		int    lucky            = 0,
		       accept_move      = 0;
		MPI_Status status;

		//make the bath index
		SafeOps::calloc( bath2core, size, sizeof(int), __LINE__, __FILE__ );
		for ( int i=0; i<size; i++ )
			for ( int j=0; j<size; j++ )
				if ( ptemp->index[j] == i ) 
					bath2core[i] = j;

		// choose the lucky bath. it's not really lucky.. this is just the first bath that we consider for swapping
		if ( !rank ) 
			lucky = (int) floor(get_rand() * size);
		MPI_Bcast( &lucky, 1, MPI_INT, 0, MPI_COMM_WORLD );

		// we will use this array to designate whether a particular core is a master or slave
		SafeOps::calloc( is_master, size, sizeof(int), __LINE__, __FILE__ );
		
		//build the partner list
		SafeOps::calloc( partner_list, size, sizeof(int), __LINE__, __FILE__ );
		master = lucky;
		slave = (lucky+1) % size; //this assigns the slave bath index to the variable
		do {
			partner_list[bath2core[slave]] = bath2core[master]; //partner_list maps slave to master core
			partner_list[bath2core[master]]= bath2core[slave]; //partner_list maps master to slave
			is_master[bath2core[master]] = 1; //master flag is on
			is_master[bath2core[slave]] = 0; //master flag is off
			//now generate the next master
			master = (master+2) % size;
			slave = (slave+2) % size;
		} while ( (master != lucky) && (slave != lucky) );
		if ( slave == lucky )  {
			partner_list[bath2core[master]] = -1; //if odd # of cores, we have one member unassigned
			is_master[bath2core[master]] = -1; //neither master nor slave
		}
		//all cores (except 1 if size is odd) are mapped

		//communicate the energy to the master nodes
		if ( ! is_master[rank] ) {
			MPI_Send(&current_energy, 1, MPI_DOUBLE, partner_list[rank], 0, MPI_COMM_WORLD); //slave sends it's energy
			MPI_Recv(&accept_move, 1, MPI_INT, partner_list[rank], 1, MPI_COMM_WORLD, &status); //receive result (accept/reject)
		}
		else if ( is_master[rank] == 1 ) {
			//master receives energy from slave
			MPI_Recv(&slave_energy, 1, MPI_DOUBLE, partner_list[rank], 0, MPI_COMM_WORLD, &status); 
			//calculate boltzmann factor exp(dE*dB)
			boltzmann_factor = exp((current_energy-slave_energy)*
				(1.0/ptemp->templist[rank] - 1.0/ptemp->templist[partner_list[rank]]));
			if ( get_rand() < boltzmann_factor ) 
				accept_move = 1;
			else
				accept_move = 0;
			//communicate the move result to slave
			MPI_Send(&accept_move, 1, MPI_INT, partner_list[rank], 1, MPI_COMM_WORLD);
		} else
			accept_move = 0; //no partner
	
		if ( accept_move ) {
			//reassign local temperature
			temperature = ptemp->templist[partner_list[rank]];
			//update our temperature and index values to prepare for transmission to root
			new_templist_val = ptemp->templist[partner_list[rank]];
			new_index_val    = ptemp->index   [partner_list[rank]];
			//reassign fugacities
			if ( ! user_fugacities && fugacities )	{
					MPI_Send(fugacities, 1, MPI_DOUBLE, partner_list[rank], 2, MPI_COMM_WORLD); //send fugacity
					MPI_Recv(fugacities, 1, MPI_DOUBLE, partner_list[rank], 2, MPI_COMM_WORLD, &status); //receive fugacity
			}				
		} else { //reject
			new_templist_val = ptemp->templist[rank];
			new_index_val    = ptemp->index[rank];
		}

		//now we need to update the templist and index on each core
		SafeOps::calloc( update_templist, size, sizeof(double), __LINE__, __FILE__ );
		SafeOps::calloc( update_index,    size, sizeof(int),    __LINE__, __FILE__ );
		
		// create updated arrays on head node
		MPI_Gather(&new_templist_val, 1, MPI_DOUBLE, update_templist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&new_index_val, 1, MPI_INT, update_index, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// build the updated array
		if ( ! rank ) {
			for ( int i=0; i<size; i++ ) {
				ptemp->templist[i] = update_templist[i];
				ptemp->index   [i] = update_index   [i];
			}
		}
		// transmit to each core
		MPI_Bcast(ptemp->templist, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(ptemp->index, size, MPI_INT, 0, MPI_COMM_WORLD);

		if ( accept_move ) 
			nodestats->accept_ptemp++;
		else if ( is_master[rank] != -1 )
			nodestats->reject_ptemp++;

	
		free(update_templist);
		free(update_index);
		free(bath2core);
		free(is_master);
		free(partner_list);
	
	*/
	return;
}




void System::do_corrtime_bookkeeping(mpiData &mpi) {
	
	// copy observables and avgs to the mpi send buffer 
	// histogram array is at the end of the message 
	if( calc_hist ) {
		zero_grid( grids->histogram->grid );
		population_histogram();
	}

	// update frozen and total system mass
	calc_system_mass();

	// update sorbate info on each node
	if ( sorbateCount > 1 )
		update_sorbate_info();

	//write trajectory files for each node -> one at a time to avoid disk congestion
	#ifdef _MPI
		for ( int j=0; j<size; j++ ) {
			MPI_Barrier(MPI_COMM_WORLD);
			if( j == rank )
				write_states();
		}
	#else
		write_states();
	#endif

	// restart files for each node -> one at a time to avoid disk congestion
	if( write_molecules_wrapper(pqr_restart) < 0 ) {
		Output::err("MC: could not write restart state to disk\n");
		throw unknown_file_error;
	}

	// dipole/field data for each node -> one at a time to avoid disk congestion
	#ifdef _MPI
		if ( polarization ) {
			for ( int j=0; j<size; j++ ) {
				MPI_Barrier(MPI_COMM_WORLD);
				if ( j == rank ) {
					write_dipole();
					write_field();
				}
			}
		}
	#else
		if ( polarization ) {
			write_dipole();
			write_field();
		}
	#endif

	// zero the send buffer 
	std::memset( mpi.snd_strct, 0, mpi.msgsize );
	std::memcpy( mpi.snd_strct, observables, sizeof(observables_t) );
	std::memcpy( mpi.snd_strct + sizeof(observables_t), avg_nodestats, sizeof(avg_nodestats_t) );
	if( calc_hist )
		mpi_copy_histogram_to_sendbuffer(
			mpi.snd_strct + sizeof(observables_t) + sizeof(avg_nodestats_t), 
			grids->histogram->grid
		);
	if ( sorbateCount > 1 )
		std::memcpy(
			mpi.snd_strct   +   sizeof(observables_t)   +   sizeof(avg_nodestats_t)   +   calc_hist * n_histogram_bins * sizeof(int), //compensate for the size of hist data, if neccessary
			sorbateInfo,
			sorbateCount * sizeof( sorbateInfo_t )
		);

	if( ! rank )
		std::memset(mpi.rcv_strct, 0, size * mpi.msgsize);

	#ifdef _MPI
		MPI_Gather( mpi.snd_strct, 1, msgtype, mpi.rcv_strct, 1, msgtype, 0, MPI_COMM_WORLD );
		MPI_Gather( &temperature, 1, MPI_DOUBLE, mpi.temperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		// need to gather data for sorbate stats also
	#else
		std::memcpy( mpi.rcv_strct, mpi.snd_strct, mpi.msgsize );
		mpi.temperature[0] = temperature;
	#endif // MPI 

	// head node collects all observables and averages 
	if( !rank ) {
		// clear avg_nodestats to avoid double-counting 
		clear_avg_nodestats();
		//loop for each core -> shift data into variable_mpi, then average into avg_observables
		for( int j = 0; j < size; j++ ) { 
			// copy from the mpi buffer 
			std::memcpy( mpi.observables,   mpi.rcv_strct + j*mpi.msgsize,                         sizeof(observables_t  ));
			std::memcpy( mpi.avg_nodestats, mpi.rcv_strct + j*mpi.msgsize + sizeof(observables_t), sizeof(avg_nodestats_t));
			if( calc_hist )
				mpi_copy_rcv_histogram_to_data( mpi.rcv_strct + j*mpi.msgsize + sizeof(observables_t) + sizeof(avg_nodestats_t), grids->histogram->grid );
			if( sorbateCount > 1 )
				std::memcpy( mpi.sinfo, 
					mpi.rcv_strct   +   j*mpi.msgsize   +   sizeof(observables_t)   +   sizeof(avg_nodestats_t)   +   calc_hist * n_histogram_bins * sizeof(int), //compensate for the size of hist data, if neccessary
					sorbateCount * sizeof(sorbateInfo_t)
				);

			// write observables 
			if( fp_energy ) 
				write_observables( fp_energy, mpi.observables, mpi.temperature[j] );
			if( fp_energy_csv ) 
				write_observables_csv( fp_energy_csv, mpi.observables, mpi.temperature[j] );
			// collect the averages 
			// if parallel tempering, we will collect obserables from the coldest bath. this can't be done for
			// nodestats though, since nodestats are averaged over each corrtime, rather than based on a single 
			// taken at the corrtime 
			update_root_nodestats( mpi.avg_nodestats, avg_observables );
			if ( ! parallel_tempering ) {
				update_root_averages( mpi.observables );
				if( calc_hist )
					update_root_histogram();
				if( sorbateCount > 1 ) 
					update_root_sorb_averages( mpi.sinfo );
			}
			else if ( ptemp->index[j] == 0 ) {
				update_root_averages( mpi.observables );
				if( calc_hist )
					update_root_histogram();
				if( sorbateCount > 1 )
					update_root_sorb_averages( mpi.sinfo );
			}
		}

		if( ensemble != ENSEMBLE_NVT_GIBBS )
			output_file_data();

	} // !rank 
}