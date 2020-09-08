#include <algorithm>
#include <cstring>
#include <iostream>
#include <map>
#include <time.h>
#include <vector>

#include "Output.h"
#include "Quaternion.h"
#include "SafeOps.h"
#include "SimulationControl.h"
#include "Vector3D.h"

#ifdef _MPI
#include <mpi.h>
#endif

#ifdef _OMP
#include <omp.h>
#endif





extern int rank, size;
extern bool mpi;

 
//  In PI ensembles, the sys variable (the main and typically only system variable in traditional runs)
//  is used to track the aggregate observables that are derived from the simulation results generated
//  in the array of systems (referenced in the vector variable systems[]). 
bool SimulationControl::PI_nvt_mc() {
	
	// For ease of reference--to aid in self-documentation
	// (the single-system 'sys' variable is used to track aggregate quantities of the entire PI system
	System                *system                    =  systems[rank];                   // System assigned to this MPI thread
	System::observables_t &PI_sys_observables        = *sys.observables;                 // observables for the aggregate PI system
	System::observables_t &PI_checkpoint_observables = *sys.checkpoint->observables;     // ^^ ditto-ish
	double                &boltzmann_factor          =  sys.nodestats->boltzmann_factor; // Boltzmann factor for current sys config
	unsigned       int    &step                      =  sys.step;                        // current MC step
	const unsigned int     nSteps                    =  sys.numsteps;                    // total MC steps to perform
	const unsigned int     sample_time               =  sys.corrtime;                    // steps between samples (i.e. correlation time)
	int                    move                      =  0;                               // current MC move 
	

	for( System* SYS : systems ) {
		SYS->observables->temperature = sys.temperature;   // distribute initial temperature to all systems
		if (SYS->cavity_bias) SYS->cavity_update_grid();   // update the grid for the first time
		SYS->observables->volume = SYS->pbc.volume;        // set volume observable
	}
	
	// If not using a config file from a previous PI run, we perform an initial bead perturbation
	// so as not to drop the molecules into the system far out of equilibrium wrt thermal wavelength
	if( ! sys.parallel_restarts )
		PI_perturb_bead_COMs_ENTIRE_SYSTEM();
	
	// Perform initial energy calculations
	PI_calculate_energy();
	
	// write initial observables to stdout and logs
	if(  ! rank  ) {

		sys.open_files(); // open output files 
		PI_calc_system_mass();
		
		// average in the initial observables values once  (we don't want to double-count the initial state when using MPI)
		average_current_observables_into_PI_avgObservables();

		// write initial observables exactly once
		if (sys.fp_energy)
			sys.write_observables(sys.fp_energy, sys.observables, sys.temperature);
		if (sys.fp_energy_csv)
			sys.write_observables_csv(sys.fp_energy_csv, sys.observables, sys.temperature);

		// output observables to stdout
		Output::out("MC: initial values:\n");
		sys.display_averages();
	}
	
	// Pick the next MC move and target molecule
	move = PI_pick_NVT_move();

	// Backup the original observables set
	backup_observables_ALL_SYSTEMS();

	// Compute the initial values that will contribute to the Boltzmann factor
	BFC.potential.current = sys.observables->potential(); //  \/ be a bit forgiving of the initial state   \/ 
	if (!std::isfinite(BFC.potential.current)) { sys.observables->energy = BFC.potential.current = MAXVALUE; }

	// The mass/mu_len2 BF contributors will only ever contribute for steps whose moves are such that PI Beads are moved relative
	// to each other. Further, the beads in one atom/molecule are completely independent wrt the beads in any other atom/molecule.
	// So... we could track the current values if we computed a sum of the PI energy for every species in the system. But since
	// only the molecule whose beads have changed will contribute, and then so only during moves that perturb bead positions
	// relative to each other, it is not necessary to track the 'current' value--we will deal with these on an as-needed
	// per-molecule;per-move basis. I.e. it is a waste of time to compute this every time, especially considering on translations
	// and rotations this will contribute nothing.
	BFC.chain_mass_len2.current = 0; // We will never use this
	BFC. orient_mu_len2.current = 0; // We will never use this
	
	
	


	 //\    Main MC Loop 
	//  \_________________________________________________________________________________________________________________

	for( step=1; step <= nSteps; step++ ) {

		// restore the last accepted energy, chain-length and orientation-chain metrics
		BFC.potential.init       = BFC.potential.current; 
		//                         We will only bother with len2 BF Contributors if beads are being perturbed
		BFC.chain_mass_len2.init = (move == MOVETYPE_PERTURB_BEADS) ? PI_chain_mass_length2()       : 0;
		BFC. orient_mu_len2.init = (move == MOVETYPE_PERTURB_BEADS) ? PI_orientational_mu_length2() : 0;
		
		// perturb the system 
		PI_make_move( move );

		// calculate new energy change, new PI chain length & update obervables
		BFC.potential.trial = PI_calculate_potential();
		BFC.chain_mass_len2.trial = (move==MOVETYPE_PERTURB_BEADS)?       PI_chain_mass_length2() : 0;
		BFC. orient_mu_len2.trial = (move==MOVETYPE_PERTURB_BEADS)? PI_orientational_mu_length2() : 0;
		
		#ifdef QM_ROTATION
			// solve for the rotational energy levels 
			if (system->quantum_rotation && (system->checkpoint->movetype == MOVETYPE_SPINFLIP))
				quantum_system_rotational_energies(system);
		#endif // QM_ROTATION 

		// treat a bad contact as a reject 
		if( ! std::isfinite(BFC.potential.trial) ) {
			BFC.potential.trial = sys.observables->energy = MAXVALUE;
			boltzmann_factor = 0;  
		} else 
			boltzmann_factor = PI_NVT_boltzmann_factor(BFC);
		
			
		
		 //\   Metropolis function 
		//  \_____________________________________________________________________________________________________________

		if(   (Rando::rand() < boltzmann_factor)   &&   (system->iterator_failed == 0)   ) {
			  
			 //\//  ACCEPT  //////////////////////////////////////////////////////////////////////////////////////////////

			sys.register_accept(move); // register an accepted move (for computing move acceptance rates)

			BFC.potential.current = BFC.potential.trial;
			
			PI_calculate_energy();
			backup_observables_ALL_SYSTEMS();
			
			// Simulated Annealing
			if (sys.simulated_annealing) {
				if (sys.simulated_annealing_linear)	{
					system->temperature = system->temperature + (system->simulated_annealing_target - system->temperature) / ((double)(nSteps) - system->step);
					if (nSteps == system->step)
						systems[rank]->temperature = systems[rank]->simulated_annealing_target;
				}
				else
					systems[rank]->temperature = systems[rank]->simulated_annealing_target + (systems[rank]->temperature - systems[rank]->simulated_annealing_target)*systems[rank]->simulated_annealing_schedule;
			}
		} else {

			//\//  REJECT: restore from last checkpoint  ///////////////////////////////////////////////////////////////
			restore_PI_systems();  // restore pre-move system configs and the associated observables
			PI_sys_observables = PI_checkpoint_observables;  // restore the aggregate PI observables 
			sys.register_reject(move);  // register a rejected move (for computing move acceptance rates)
		}

		
		sys.compile_MC_algorithm_stats(); // track the acceptance_rates/BF and compute associated stats
		move = PI_pick_NVT_move();        // Pick the move (and target) for the next MC step
	


		//\//  Do this every correlation time, and at the very end   /////////////////////////////////////////////////////
		if(  !(step % sample_time)   ||   (step == nSteps)) {
			do_PI_corrtime_bookkeeping();
		}
	} // main loop 

	// Simulation has finished...
	if (!rank) {
		int sid = 0;
		for( System *SYS : systems ) {
			if (SYS->write_molecules_wrapper(SYS->pqr_output) < 0) {
				char linebuff[2*maxLine];
				sprintf(linebuff, "System: %d\nFilename: %s\n", sid, SYS->pqr_output);
				Output::err(linebuff);
				Output::err("MC: could not write final state to disk\n");
				throw unknown_file_error;
			}
			sid++;
		}
	}
	return ok;
}




void SimulationControl::restore_PI_systems() {
	for( System *SYS : systems ) {
		SYS->iterator_failed = 0; // reset the polar iterative failure flag
		SYS->restore();           // restore old system configuration and associated observables
	}
}




void SimulationControl::average_current_observables_into_PI_avgObservables() {
	
	// Create a reference to the "main system" geometry/simulation-box in the system variable that tracks the
	// Path Integral observables (some of the observables use this information to compute said observables)
	sys.molecules = systems[rank]->molecules;
	
	// Update periodic boundary conditions. Only needed once in in PI_NVT, but would require updating in, say PI_NPT
	sys.pbc = systems[rank]->pbc;

	// Update values in the PI observable tracker that are not updated during the course of the simulation...
	sys.observables->N           = systems[rank]->observables->N;
	sys.observables->volume      = systems[rank]->observables->volume;
	sys.observables->temperature = systems[rank]->observables->temperature;
	sys.observables->spin_ratio  = systems[rank]->observables->spin_ratio;
	sys.observables->NU          = systems[rank]->observables->NU;
	
	// ...and average them (and the other observables) into the mix:
	sys.update_root_averages(sys.observables);

	// Dereference the system geometry from the PI tracker, so the memory doesn't get accidentally corrupted or freed.
	sys.molecules = nullptr;
}




void SimulationControl::do_PI_corrtime_bookkeeping() {

	System *system = systems[rank];

	// copy observables and avgs to the mpi send buffer; histogram array is at the end of the message
	for( System *SYS : systems ) {
		if (SYS->calc_hist) {
			SYS->zero_grid(SYS->grids->histogram->grid);
			SYS->population_histogram();
		}

		// update frozen and total system mass
		SYS->calc_system_mass();

		// update sorbate info on each node
		if (SYS->sorbateCount > 1)
			SYS->update_sorbate_info();
	}
	sys.observables->total_mass  = system->observables->total_mass;
	sys.observables->frozen_mass = system->observables->frozen_mass;

	
	if(  ! rank   &&   PI_xyz_corrtime_frames_requested()) 
		write_PI_frame();

	// write observables
	if (sys.fp_energy)
		sys.append_observables_to_output_file();
	if (sys.fp_energy_csv)
		sys.append_observables_to_csv_file();

	sys.clear_avg_nodestats();
	sys.update_root_nodestats();
	average_current_observables_into_PI_avgObservables();
	
	if (sys.calc_hist)
		sys.update_root_histogram();

	if (sys.sorbateCount > 1)
		sys.update_root_sorb_averages(system->mpi_data.sinfo);

	sys.output_file_data();

	

	// Write the state file, restart files, and dipole/field files
	for (int sysID = 0; sysID < nSys; sysID++) {

		#ifdef _MPI
			// Write output files, one at a time to avoid disk congestion
			if (mpi) { MPI_Barrier(MPI_COMM_WORLD); }
		#endif // _MPI

		// *THIS* thread will write the states for every system on non-MPI runs
		// OTOH, MPI threads will only write the state for *its* thread
		if (  ! mpi   ||   sysID==rank  ) {

			// Determine the system of interst for this iteration
			// systems[sysID] should normally work for both cases, but just in case...
			System* SYS = (mpi ? system : systems[sysID]); 
			
			// write the trajectory 
			SYS->write_states();

			// write the restart geometry
			if (SYS->write_molecules_wrapper(SYS->pqr_restart) < 0) {
				Output::err("MC: could not write restart state to disk\n");
				throw unknown_file_error;
			}

			// write the dipole/field data for each node
			if (sys.polarization) {
				SYS->write_dipole();
				SYS->write_field();
			}
		} 
	}
	
	
	//  Original treatment follows:
	/*

	// Write the state file, restart files, and dipole/field files
	if (mpi) {

		#ifdef _MPI
			for (int j = 0; j < size; j++) {
			
				MPI_Barrier(MPI_COMM_WORLD);
				// Write output files, one at a time to avoid disk congestion
				if (j == rank) {

					// write trajectory files for each node
					system->write_states();

					// write restart files for each node
					if (system->write_molecules_wrapper(system->pqr_restart) < 0) {
						Output::err("MC: could not write restart state to disk\n");
						throw unknown_file_error;
					}

					// dipole/field data for each node
					if (sys.polarization) {
						systems[j]->write_dipole();
						systems[j]->write_field();
					}
				}
			}
		#endif
	}
	else {

		system->write_states();

		std::for_each(systems.begin(), systems.end(), [](System *SYS) {

			// write restart files for each system
			if (SYS->write_molecules_wrapper(SYS->pqr_restart) < 0) {
				Output::err("MC: could not write restart state to disk\n");
				throw unknown_file_error;
			}

			// write dipole/field data for each system
			if (SYS->polarization) {
				SYS->write_dipole();
				SYS->write_field();
			}
		});
	}


	// Copy the data we wish to share into the send buffers for the respective systems
	if (mpi) {
		// Copy the observables and the average nodestats for this system to its send buffer
		std::memset(system->mpi_data.snd_strct, 0, system->mpi_data.msgsize);
		std::memcpy(system->mpi_data.snd_strct, system->observables, sizeof(System::observables_t));
		std::memcpy(system->mpi_data.snd_strct + sizeof(System::observables_t), system->avg_nodestats, sizeof(System::avg_nodestats_t));

		// Ditto for the histogram and sorbate data (if we are using these systems) 
		if (sys.calc_hist)
			system->mpi_copy_histogram_to_sendbuffer(
				system->mpi_data.snd_strct + sizeof(System::observables_t) + sizeof(System::avg_nodestats_t),
				system->grids->histogram->grid
			);
		if (sys.sorbateCount > 1)
			std::memcpy(
				system->mpi_data.snd_strct + sizeof(System::observables_t) + sizeof(System::avg_nodestats_t) + (size_t) system->calc_hist * system->n_histogram_bins * sizeof(int), //compensate for the size of hist data, if neccessary
				system->sorbateInfo,
				system->sorbateCount * sizeof(System::sorbateInfo_t)
			);
	} 
	else std::for_each(systems.begin(), systems.end(), [](System *SYS) {
		
		// Copy the observables and average nodestats for this system into its send buffer
		std::memset(SYS->mpi_data.snd_strct, 0, SYS->mpi_data.msgsize);
		std::memcpy(SYS->mpi_data.snd_strct, SYS->observables, sizeof(System::observables_t));
		std::memcpy(SYS->mpi_data.snd_strct + sizeof(System::observables_t), SYS->avg_nodestats, sizeof(System::avg_nodestats_t));

		// Ditto for the histogram data and sorbate data (if we are using these systems)
		if (SYS->calc_hist)
			SYS->mpi_copy_histogram_to_sendbuffer(
				SYS->mpi_data.snd_strct + sizeof(System::observables_t) + sizeof(System::avg_nodestats_t),
				SYS->grids->histogram->grid
			);
		if (SYS->sorbateCount > 1)
			std::memcpy(
				SYS->mpi_data.snd_strct + sizeof(System::observables_t) + sizeof(System::avg_nodestats_t) + (size_t) SYS->calc_hist * SYS->n_histogram_bins * sizeof(int), //compensate for the size of hist data, if neccessary
				SYS->sorbateInfo,
				SYS->sorbateCount * sizeof(System::sorbateInfo_t)
			);
	});
	
	// Clear the receive buffer
	if( ! rank )
		std::memset(system->mpi_data.rcv_strct, 0, (size_t) (mpi?size:nSys) * system->mpi_data.msgsize);

	// Copy observables, average nodestats, histogram data, etc into the receive buffer
	if (mpi) {
		#ifdef _MPI
			// copy all data into the receive struct of the head node...
			MPI_Gather(   system->mpi_data.snd_strct, 1, system->msgtype, system->mpi_data.rcv_strct,   1, system->msgtype, 0, MPI_COMM_WORLD);
			MPI_Gather( &(system->temperature),       1,      MPI_DOUBLE, system->mpi_data.temperature, 1,      MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif
	} else {
		for (size_t i = 0; i < nSys; i++) {
			// ...or if single-threaded, copy data from the each of the send structs into the appropriate part of the receive struct manually
			std::memcpy( system->mpi_data.rcv_strct + i * system->mpi_data.msgsize, systems[i]->mpi_data.snd_strct, system->mpi_data.msgsize);
			if(sys.parallel_tempering)
				system->mpi_data.temperature[i] = systems[i]->temperature;
		}
	}

	// head node collects all observables and averages
	if( ! rank ) {

		// clear avg_nodestats to avoid double-counting
		system->clear_avg_nodestats();

		// loop for each core -> shift data into variable_mpi, then average into avg_observables
		for( size_t j = 0; j <  (mpi?size:nSys); j++ ) {
			// copy from the mpi buffer
			std::memcpy( systems[0]->mpi_data.observables,   systems[0]->mpi_data.rcv_strct  +  j * systems[0]->mpi_data.msgsize,  sizeof(System::observables_t )); // copy observable data into observables struct, for one system at a time
			std::memcpy( systems[0]->mpi_data.avg_nodestats, systems[0]->mpi_data.rcv_strct  +  j * systems[0]->mpi_data.msgsize + sizeof(System::observables_t), sizeof(System::avg_nodestats_t )); // ditto for avg_nodestates

			if( sys.calc_hist )
				system->mpi_copy_rcv_histogram_to_data( systems[0]->mpi_data.rcv_strct  +  j * systems[0]->mpi_data.msgsize  +  sizeof(System::observables_t) + sizeof(System::avg_nodestats_t), systems[0]->grids->histogram->grid ); // ditto for histogram grid
			if( sys.sorbateCount > 1 )
				std::memcpy( systems[0]->mpi_data.sinfo, // ditto for sorbateInfo
					systems[0]->mpi_data.rcv_strct   +   j * systems[0]->mpi_data.msgsize   +   sizeof(System::observables_t)   +   sizeof(System::avg_nodestats_t)   + sizeof(int) * sys.calc_hist * sys.n_histogram_bins, //compensate for the size of hist data, if neccessary
					sys.sorbateCount * sizeof(System::sorbateInfo_t)
				);
			
			// write observables
			if( system->fp_energy )
				system->write_observables( system->fp_energy, system->mpi_data.observables, system->mpi_data.temperature[j] );
			if( system->fp_energy_csv )
				system->write_observables_csv( system->fp_energy_csv, system->mpi_data.observables, system->mpi_data.temperature[j] );

			// collect the averages
			// if parallel tempering, we will collect obserables from the coldest bath. this can't be done for nodestats though,
			// since nodestats are averaged over each corrtime, rather than based on a single sample taken at the corrtime
			system->update_root_nodestats( system->mpi_data.avg_nodestats, system->avg_observables );
			
			if ( ! sys.parallel_tempering ) {
				/////////////////////////////////////////////////////////////////////////////////////////////
				// This is what needs to be updated for PI functionality:   /////////////////////////////////
				system->update_root_averages(system->mpi_data.observables); // system->avg_observables );
				
				if( sys.calc_hist )
					system->update_root_histogram();
				if( sys.sorbateCount > 1 )
					system->update_root_sorb_averages( system->mpi_data.sinfo );
				/////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////////////////////
			}
			else if ( system->ptemp->index[j] == 0 ) {
				system->update_root_averages( system->mpi_data.observables ); // , system->avg_observables );
				if( sys.calc_hist )
					system->update_root_histogram();
				if( sys.sorbateCount > 1 )
					system->update_root_sorb_averages( system->mpi_data.sinfo );
			}
		}

		system->output_file_data();

	} // !rank
	*/
}




double SimulationControl::PI_NVT_boltzmann_factor( PI_NVT_BFContributors BF ) {

	double delta_energy  = BF.potential.change();
	double delta_chain   = BF.chain_mass_len2.change();
	double delta_orient  = BF.orient_mu_len2.change();

	const int     P = nSys; // Trotter number: number of beads (or "systems") in PI representation
	const double  T = sys.temperature;
	double        boltzmann_factor = 0;


	switch (systems[rank]->checkpoint->movetype) {

	case MOVETYPE_PERTURB_BEADS:
	{
		// Conversion factor that will cast the PI chain length [(dAB)^2 + (dBC)^2 + (dCD)^2 ...] to energy (in Kelvin)
		// Part of this factor is (1/thermal_wavelength), with the mass factored out, because the COM contribution 
		// depends on the molecular mass, while the orientational factor depends on the reduced mass. Said masses are
		// baked into delta_chain and delta_orient. 
		const double PIchain_2_K = (P * pi * pi * kB * T) / (2.0 * h * h); //  (thermal wavelength)^(-2)  Without mass/reduced_mass
		                                                                   //  (bc it is baked into the delta variables)
		double potential_contrib      = delta_energy / T;           //  potential energy contribution
		double PI_COM_contrib         = delta_chain * PIchain_2_K;  //  PI COM energy contribution
		double PI_orientation_contrib = 0;                          //  PI orientational energy contribution
		
		std::map<std::string, int>::iterator it;
		it = sorbate_data_index.find(systems[rank]->checkpoint->molecule_altered->moleculetype);
		if (it != sorbate_data_index.end()) {
			// calculation for PIs with orientational degree of freedom
			double reduced_mass = sorbate_data[it->second].reduced_mass;
			PI_orientation_contrib = delta_orient * PIchain_2_K;
		}

		boltzmann_factor = exp( -potential_contrib - PI_COM_contrib - PI_orientation_contrib );
	}
	break;

	case MOVETYPE_SPINFLIP:
	{
		double partfunc_ratio = 0;
		double g = systems[rank]->checkpoint->molecule_altered->rot_partfunc_g;
		double u = systems[rank]->checkpoint->molecule_altered->rot_partfunc_u;
		if (systems[rank]->checkpoint->molecule_altered->nuclear_spin == NUCLEAR_SPIN_PARA)
			partfunc_ratio = g / (g + u);
		else
			partfunc_ratio = u / (g + u);

		// set the boltz factor, including ratio of partfuncs for different symmetry rotational levels
		boltzmann_factor = partfunc_ratio;
	}
	break;

	default: // DISPLACE
		boltzmann_factor = exp(-delta_energy/T);

	}
	
	return boltzmann_factor;
}




bool SimulationControl::check_PI_options() {
	
	//  Check to see if the bead count is a power of 2
	//  (i.e. count the number of ones in the binary representation of 'size' or nSys and make sure there is only 1)

	unsigned int bits = 8 * (mpi ? sizeof(size) : sizeof(nSys));
	unsigned int bitcount = 0;
	unsigned long long int bitmask = 1; 
	char linebuf[maxLine];

	for( unsigned int i = 0; i < bits; i++ ) {
		if( nSys & bitmask )
			bitcount++;
		bitmask = bitmask << 1;
	}

	if(  (nSys<4) || (bitcount != 1)  ) {
		if (mpi) {
			Output::out("SIMULATION CONTROL: Path Integrals require at least 4 MPI processes to run. One process per PI bead and a total of 2^N 'beads' (N >= 2).\n");
			sprintf(linebuf, "SIMULATION CONTROL: MPI-Reported world-size: %d.\n", size);
			Output::err(linebuf);
		} else {
			Output::err("Path intergrals require at least 4 systems to run--Trotter number (P) must be set to a power of 2 greater than 4. The\n");
			Output::err("Set the trotter number(P) with the - P X command line option. E.g.: mpmc++ -P 8 my_input_file\n");
		}
		Output::err("The Trotter number (P) sets the number of individual simulations that will be conducted that, in aggregate, will represent\n");
		Output::err("the single quantum system. Due to constraints in the orientation algorithm, the Trotter number must be a power of 2.\n");
		throw invalid_MPI_size_for_PI;
	}

	// Make sure the desired length of trial chains for bead perturbations is set appropriately

	if( ! PI_trial_chain_length ) {
		Output::err("SIMULATION CONTROL: PI_trial_chain_length must be set when using Path Integral ensembles.\n");
		throw invalid_setting;
	}
	if( (PI_trial_chain_length < 0)  ||  (PI_trial_chain_length >= nSys) ){
		Output::err( "SIM_CONTROL: PI_trial_chain_length must be in [1..P-1], where P is the Trotter number,\n"   );
		Output::err( "             i.e. the number of 'beads' (1 bead per MPI thread). For a non-MPI runs, the\n" );
		Output::err( "             bead count is set with the -P option on the command line.\n" );
		if(mpi) 
			sprintf(linebuf, "SIM_CONTROL: MPI-Reported world-size (i.e. PI bead count): %d.\n", size);
		else 
			sprintf(linebuf, "SIM_CONTROL: user requested bead count: %d\n", nSys );
		Output::err(linebuf);
		sprintf(linebuf, "SIM_CONTROL: requested length of trial chain: %d.\n", PI_trial_chain_length);
		Output::err(linebuf);
		throw invalid_setting;
	}
	
	return ok;
}




void SimulationControl::initialize_PI_NVT_Systems() {

	char linebuf[maxLine] = {0};
	Output::out1("\n");

	// Create a system object for each bead on the Path Integral loop, and populate those systems with the system geometry. 

	for( int i=0; i<nSys; i++ ) {
		
		systems.push_back( new System(sys) );
		
		if( sys.parallel_restarts)
			strcpy(systems[i]->pqr_input, pqr_input_filenames[i].c_str());
		strcpy(systems[i]->pqr_output,    pqr_final_filenames[i].c_str());
		strcpy(systems[i]->pqr_restart, pqr_restart_filenames[i].c_str());
		sprintf(linebuf, "SIM_CONTROL: SYSTEM[ %d ] Instantiated.\nSIM_CONTROL->SYSTEM[ %d ]: Constructing simulation box.\n", i, i );
		Output::out1(linebuf);
		
		systems[i]->setup_simulation_box();
		sprintf(linebuf, "SIM_CONTROL->SYSTEM[ %d ]: simulation box configured.\n\n", i);
		Output::out1(linebuf);
	
		systems[i]->allocateStatisticsMem();

		if (systems[i]->calc_hist)
			systems[i]->setup_histogram();

		if (systems[i]->cavity_bias)
			systems[i]->setup_cavity_grid();

		// ensure that all SPECTRE charges lie within the restricted domain 
		if (systems[i]->spectre)
			systems[i]->spectre_wrapall();

		if (mpi && (i != rank)) continue;
		// The initialization that follows supports energy calculations, which will only be done on 
		// rank's system (for multi-threaded runs) and so will only be instantiated on rank's data
		// structure. For single threaded runs these features will need to be enabled on all systems.
		///////////////////////////////////////////////////////////////////////////////////////////////////
		

		// Set up the pairs list etc. This must happen for every system on a single thread, but only for the system (i.e. 
		// systems[rank]) that this controller will be manipulating & for which it will be performing energy computations. 
		systems[i]->allocate_pair_lists();
		systems[i]->pairs();          // get all of the pairwise interactions, exclusions, etc.
		systems[i]->flag_all_pairs(); // set all pairs to initially have their energies calculated 
		
		
		// if polarization active, allocate the necessary matrices 
		if(  systems[i]->polarization  &&  ( ! systems[i]->cuda )  &&  ( ! systems[i]->polar_zodid )  )
			systems[i]->thole_resize_matrices();		
	}

	// Allocate memory for MPI-energy-transfers / multi-system-energy-bookkeeping / MC algo stats
	//\\///////////////////////////////////////////////////////////////////////////////////////////////////////
	SafeOps::calloc(          rd_energies, nSys, sizeof(double), __LINE__, __FILE__);
	SafeOps::calloc(   coulombic_energies, nSys, sizeof(double), __LINE__, __FILE__);
	SafeOps::calloc(polarization_energies, nSys, sizeof(double), __LINE__, __FILE__);
	SafeOps::calloc(         vdw_energies, nSys, sizeof(double), __LINE__, __FILE__);
	
	SafeOps::calloc( sys.observables,             1, sizeof(System::observables_t),     __LINE__, __FILE__);
	SafeOps::calloc( sys.avg_observables,         1, sizeof(System::avg_observables_t), __LINE__, __FILE__);
	SafeOps::calloc( sys.nodestats,               1, sizeof(System::nodestats_t),       __LINE__, __FILE__);
	SafeOps::calloc( sys.avg_nodestats,           1, sizeof(System::avg_nodestats_t),   __LINE__, __FILE__);
	SafeOps::calloc( sys.checkpoint,              1, sizeof(System::checkpoint_t),      __LINE__, __FILE__);
	SafeOps::calloc( sys.checkpoint->observables, 1, sizeof(System::observables_t),     __LINE__, __FILE__);
	//\\///////////////////////////////////////////////////////////////////////////////////////////////////////


	#ifdef _MPI
		if (mpi) MPI_Barrier(MPI_COMM_WORLD);
	#endif	

	Output::out1("SIM_CONTROL: finished allocating pair lists\n");
	Output::out1("SIM_CONTROL: finished calculating pairwise interactions\n");
	if( ! sys.use_sg || sys.rd_only ) {
		sprintf(linebuf, "SIM_CONTROL: Ewald gaussian width = %f A\n", sys.ewald_alpha);
		Output::out1(linebuf);
		sprintf(linebuf, "SIM_CONTROL: Ewald kmax = %d\n", sys.ewald_kmax);
		Output::out1(linebuf);
	}
}




void SimulationControl::write_PI_frame() {
// output an XYZ format frame that records positions of sites in all systems

	int nSites = nSys * systems[0]->countNatoms();
	
	// we will append to the existing file...
	char write_mode[2] = { 'a', '\0' };
	
	// ... unless this is the first frame, in which case we will start a new file
	const int very_first_frame_number = 1;
	static int frame_number = very_first_frame_number;
	if( frame_number == very_first_frame_number )
		write_mode[0] = 'w';

	// generate the filename and open the file
	FILE *outFile = fopen( PI_frames_filename, write_mode );

	fprintf(outFile, "%d\nFrame: %d\n", nSites, frame_number);

	++frame_number;

	// Record the position of every atom in every molecule in every PI system
	Molecule * m;
	Atom * a;
	for (int s = 0; s < nSys; s++) 
		for (m = systems[s]->molecules; m; m = m->next) 
			for (a = m->atoms; a; a = a->next) 
				fprintf(outFile, "%s     %0.4lf     %0.4lf     %0.4lf\n", a->atomtype, a->pos[0], a->pos[1], a->pos[2]);
	

	fclose(outFile);
	
}




double SimulationControl::PI_calculate_energy() {
	
	// Tuckerman (2010) Statistical Mechanics: Theory and Molecular Simulation.
	double kinetic   = PI_calculate_kinetic  ();    // terms 1 & 2, energy estimator.  (12.5.12)
	double potential = PI_calculate_potential();    // term 3, energy estimator.       (12.5.12)

	#ifdef QM_ROTATION
		// solve for the rotational energy levels 
		// if (systems[rank]->quantum_rotation) 
		//     quantum_system_rotational_energies(systems[rank]);
	#endif
	sys.observables->energy = kinetic + potential;

	return sys.observables->energy;
}



double SimulationControl::PI_calculate_potential() {

	// generic shorthand for whatever set of observables is relevant at the moment:
	System::observables_t* obs; 

	if (mpi) {

		systems[rank]->energy();
		obs = systems[rank]->observables;
		
		#ifdef _MPI
			MPI_Allgather(           &obs->rd_energy, 1, MPI_DOUBLE,           rd_energies, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Allgather(    &obs->coulombic_energy, 1, MPI_DOUBLE,    coulombic_energies, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Allgather( &obs->polarization_energy, 1, MPI_DOUBLE, polarization_energies, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Allgather(          &obs->vdw_energy, 1, MPI_DOUBLE,          vdw_energies, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		#endif
		
	} else {
		// For single-MPI-thread systems, energy computations happen on every system
		
		#pragma omp parallel for num_threads(nSys)
		for (int s = 0; s < nSys; s++) {
			systems[s]->energy();
			rd_energies[s]           = systems[s]->observables->rd_energy;
			coulombic_energies[s]    = systems[s]->observables->coulombic_energy;
			polarization_energies[s] = systems[s]->observables->polarization_energy;
			vdw_energies[s]          = systems[s]->observables->vdw_energy;
		}
	}

	// Observables for the aggregate Path Integral (PI) system will be accumulated in sys.observables
	
	obs = sys.observables;

	obs->rd_energy           = 0;
	obs->coulombic_energy    = 0;
	obs->polarization_energy = 0;
	obs->vdw_energy          = 0;

	for( int s=0; s < nSys; s++) {
		obs->rd_energy           += rd_energies[s];
		obs->coulombic_energy    += coulombic_energies[s];
		obs->polarization_energy += polarization_energies[s];
		obs->vdw_energy          += vdw_energies[s];
	}
	
	obs->rd_energy           /= nSys;
	obs->coulombic_energy    /= nSys;
	obs->polarization_energy /= nSys;
	obs->vdw_energy          /= nSys;
	
	return   sys.observables->rd_energy  + sys.observables->coulombic_energy
	       + sys.observables->vdw_energy + sys.observables->polarization_energy;
}




double SimulationControl::PI_calculate_kinetic() {
	const double d      = 3.0;                          // dimensionality of the system
	const double N      = systems[0]->countN();         // Number of sorbate (or moveable) molecules in the system
	const double P      = nSys;
	const double T      = sys.temperature;
	const double beta   = 1.0 / (kB * T);
	const double omega2 = P / (beta * beta * hBar2);    // omega squared

	double chain_mass_len2 =       PI_chain_mass_length2_ENTIRE_SYSTEM();
	double  orient_mu_len2 = PI_orientational_mu_length2_ENTIRE_SYSTEM();
	
	// Tuckerman (2010) Statistical Mechanics: Theory and Molecular Simulation.   Eq (12.5.12)
	double energy_estimator_term_1 = 0.5 * d * N * kB * T  *  P;     // equipartition energy [12.5.12]
	double energy_estimator_term_2 = 0.5 * omega2 * chain_mass_len2; // [12.5.12]

	sys.observables->kinetic_energy = (1.0/kB) * (energy_estimator_term_1 - energy_estimator_term_2); // converted to Kelvin
	
	return sys.observables->kinetic_energy;
} 




void SimulationControl::PI_calc_system_mass() {
	systems[rank]->calc_system_mass();
	sys.observables->frozen_mass = systems[0]->observables->frozen_mass;
	sys.observables-> total_mass = systems[0]->observables-> total_mass;
}




void SimulationControl::assert_all_perturb_targets_exist() {
	// Ensure that a molecule has been targeted for perturbation.
	for_each(systems.begin(), systems.end(), [](System* s) {
		if (s->checkpoint->molecule_altered == nullptr) {
			Output::err("System corrupted. No target molecule specified for PI chain measurement request.");
			throw internal_error;
		}
		});
}




// In the following functions, "chain mass length2" is a slight misnomer. If A, B and C are Molecule objects which, together, 
// comprise a PI/quantum representation of a single molecule, AND if dAB is the distance between COMs of molecules A & B, AND
// if M is the mass of said molecule, then the central quantity in these functions is M * (dAB*dAB + dBC*dBC + dCA*dCA).
// PI_chain_mass_length2_ENTIRE_SYSTEM() computes the sum of these values for every molecule in the system. 
double SimulationControl::PI_chain_mass_length2_ENTIRE_SYSTEM() {

	double sum = 0;
	Molecule* molPtr = nullptr;
	std::vector<Molecule*> molecules_ptr(nSys);  // this array of ptrs tracks the current molecule in each system 

	// Populate our vector with pointers to the first molecule in each of the images comprising the system. 
	for (int s = 0; s < nSys; s++) 
		molecules_ptr[s] = systems[s]->molecules;
	
	// Step through each perturb-able molecule in lockstep.
	//   That is, with each loop iteration, the entire vector of pointers should reference a different instance of
	//   what is conceptually the same physical molecule in each system image. E.g., if molecules_ptr[0] points to
	//   molecule 3 in the first system (0), then molecules_ptr[5] should point to sorbate 3 in the sixth system (5).
	//   In the next iteration, these same pointers should all point to sorbate 4 in their respective systems.

	while (molPtr = molecules_ptr[0]) {

		// Check to see if the current molecule should be perturbed
		// Corresponding molecules *should* be the same across all systems, so only the version of the molecule in
		// system 0 (molecules_ptr[0] AKA molptr) is checked for validity
		if (!(molPtr->frozen || molPtr->adiabatic || molPtr->target)) {

			// If the molecule in image/sys 0 looks good, verify that each member is pointing to something (anything!)
			for (int s = 0; s < nSys; s++) {
				if (nullptr == molecules_ptr[s]) {
					Output::err("ERROR: System images are not consistent, head system has more molecules than at least one other system image.\n");
					throw internal_error;
				}
			}
			// Compute the chain length for this molecule and add it to our running total for the entire system
			sum += PI_chain_mass_length2(molecules_ptr);
		}

		// advance our collection of molecule pointers to the next molecule in their respective system image
		for (int s = 0; s < nSys; s++)
			molecules_ptr[s] = molecules_ptr[s]->next;
	}

	return sum;
}
double SimulationControl::PI_chain_mass_length2() {
// This function only computes the chain length of the molecule targeted for perturbation, that is, the molecule
// pointed to by system->checkpoint->molecule_altered. For Boltzman Factor computations, this is the only measurement
// required, as the lengths of all non-changing molecular PI chains will cancel in the Boltzmann factor. Accordingly, if
// no target is specified something has gone wrong (wrt to the intended context in which this fxn was meant to be called). 
	assert_all_perturb_targets_exist();
	std::vector<Molecule*> mol;
	for_each(systems.begin(), systems.end(), [&mol](System* s) {
		mol.push_back(s->checkpoint->molecule_altered);
	});
	return  PI_chain_mass_length2(mol);
}
double SimulationControl::PI_chain_mass_length2( std::vector<Molecule*> &molecule ) {
	// 'molecule' should be a vector such that each element points  to the the image of a single molecule, as it appears
	// in each of the P systems. The collection of said molecules, in aggregate, comprises the PI/quantum representation
	// of a single real-world molecule being simulated.
	
	// To compute the mass-weighted length of the "polymer chain", we must compute a harmonic-well-type potential between adjacent 
	// images of congruent molecular species. Systems 1 thru P should all have the same number and type of molecular constiutents.
	// Molecule 1 in System 1 is a partial representation of Molecule 1. Molecule 1 in System 2 is another part of that representation
	// and each system can be thought of as a "many worlds" representation where each member exists in every system but is affected 
	// by slightly different circumstances. The complete representation of Molecule 1 is an average over all the "Molecule 1"'s in 
	// every system. The "polymer chain" of constituent "beads" for each molecule are connected by a harmonic potential between each
	// of them that makes said images form a coherent representation of a single particle and prevents the beads from drifting apart
	// independently of each other. The PI "chain length2" for a single molecule is a loop of harmonic "potentials" between the 
	// COMs of the adjacent molecule images that together constitute the representation of that molecule. E.g.:

	// Mol1COM,Sys1  \___/  Mol1COM,Sys2  \___/  Mol1COM,Sys3  \___/  Mol1COM,Sys1 (and we arrive back at the first system)
	
	// Finally, the PI measure of the distance is weighted by the mass of the molecule in question.

	double PI_chain_mass_length2 = 0;   // "length", for lack of a better word. It is the mass*(distance^2 + distance^2 + ...)
	                                    // for the molecule whose bead representation was perturbed. The distances measure the  
	                                    // separations between adjacent bead COMs on the PI "polymer chain". It is actually some
	                                    // sort of weighted length quantity, but it plays the same role as energy if we were
	                                    // simulating a harmonic potential.


	// record the COM coordinates of each system's version of the target molecule.
	std::vector<Vector3D> COMs;
	for( Molecule *m : molecule ) {
		m->update_COM();
		Vector3D com(m->com[0], m->com[1], m->com[2]);
		COMs.push_back(com);
	}
	
	// Now, 'molecules' points to the respective version of the same molecule in each of the nSys systems. We have the COM
	// coords of each different image of that molecule recorded in x,y,z; so we effectively have the coords of the 'loop'
	// representation of a single molecule in the x,y,z arrays. The loop being formed in the manner, e.g.:

	// x[0] <~> x[1] <~> x[2] <~> x[0]      |  where ( x[0], y[0], z[0] ) is the COM position for the current molecule in
	// y[0] <~> y[1] <~> y[2] <~> y[0]      |  the first system. And where the COM position for the SAME molecule in the 
	// z[0] <~> z[1] <~> z[2] <~> z[0]      |  second system is ( x[1], y[1], z[1] ).

	// Cycle around the COM coordinate loop and sum the squared distance for each adjacent pair, i.e. we will be 
	// finding dist^2 values for each pair of coordinates (x[i],y[i],z[i]) & (x[i+1], y[i+1], z[i+1]), such that the
	// last xyz coords in the list will be paired with the first.

	for (int i = 0; i < nSys; i++) {
		int j = (i + 1) % nSys;
		Vector3D delta = COMs[i] - COMs[j];
		PI_chain_mass_length2 += delta.norm2();
	}
	PI_chain_mass_length2 *= (molecule[0]->mass * AMU2KG) * (ANGSTROM2METER * ANGSTROM2METER);   //  weight chain by the mass and
	                                                                                             //  convert A^2 to m^2.
	return PI_chain_mass_length2;
}




double SimulationControl::PI_orientational_mu_length2_ENTIRE_SYSTEM() {
	return 0.0;
}
double SimulationControl::PI_orientational_mu_length2() {
// If each bond is described as a vector, this function computes the difference vector between bonds on adjacent molecules in the PI
// bead chain. It takes the squared-norm of all these differences and returns the their sum. IT ONLY COMPUTES this difference on the
// the molecule targeted for perturbation (as this quantity will cancel for all the undisturbed PI chains throughout the system. As
// such, if no target is specified, something has gone wrong.	
	
	double PI_orient_diff = 0.0;   // "orientational difference", for lack of a better word. It is the (length^2 + length^2 + ...)
								   // for each of P difference vectors describing the spatial difference between bonds in adjacent 
	                               // images of a PI representation of a diatomic molecule in the system. 

	// To compute the weighted "orientational distance" of the "polymer chain", we must compute the difference vector between adjacent
	// bonds in our diatomic molecule. Systems 1 -> P should all have the same number and type of molecular constiutents. Molecule 1 in
	// System 1 is a partial representation of Molecule 1. Molecule 1 in System 2 is another part of that representation and each system
	// can be thought of as a "many worlds" representation where each member exists in every system but is affected by slightly different
	// circumstances. The complete representation is an average over all the "Molecule 1"'s in every system. But the "orientational 
	// distance" is a harmonic potential of sorts (as a function of the length of the difference vector) between the bond-vectors of
	// each of these images.

	
	std::vector<Vector3D> bond_vectors;
	char  *moleculeID       = systems[rank]->checkpoint->molecule_altered->moleculetype;
	int    orientation_site = SimulationControl::get_orientation_site( moleculeID );
	double bond_length      = SimulationControl::get_bond_length(      moleculeID );
	if (  (orientation_site < 0)   ||   (bond_length <= 0)  )
		return 0.0;

	// form the difference vectors 
	for_each( systems.begin(), systems.end(), [orientation_site, bond_length, & bond_vectors](System * s) {
		
		// Get COM for perturbation target
		Molecule *molecule = s->checkpoint->molecule_altered;
		molecule->update_COM();
		Vector3D rCOM( molecule->com[0], molecule->com[1], molecule->com[2] );

		// Create a vector locating the atomic site that is the designated "handle" for orienting
		// the molecule. This is what we refer to when we reference the "orientation site". 
		int site = 0;
		Atom *aPtr = nullptr;
		for (aPtr = molecule->atoms; site != orientation_site; aPtr = aPtr->next)
			site++;
		Vector3D handle_position( aPtr->pos[0], aPtr->pos[1], aPtr->pos[2]);

		// find the difference between the handle and the COM (the bond direction), and scale it such that the
		// length of this vector equals the bond length of the molecule (thereby representing the bond itself)
		Vector3D bond = handle_position - rCOM;
		bond = bond_length * bond.normalize();
		bond_vectors.push_back( bond );
	});
	
	// Cycle once through the bond vectors and sum the square-length of the difference vector between each adjacent pair, 
	// i.e. we will be finding dist^2 values for each pair of coordinates (x[i],y[i],z[i]) & (x[i+1], y[i+1], z[i+1]), 
	// such that the last difference will be between the last bond vector and the first.
	for (int i = 0; i < nSys; i++) {
		int j = (i + 1) % nSys;
		Vector3D diff_vector = bond_vectors[i] - bond_vectors[j];
		PI_orient_diff += diff_vector.norm2();
	}
	PI_orient_diff *= (ANGSTROM2METER * ANGSTROM2METER); // convert A^2 to m^2

	return PI_orient_diff;
}
double SimulationControl::PI_orientational_mu_length2(std::vector<Vector3D*> &o) {
	return 0;
}




int SimulationControl::PI_pick_NVT_move() {
// PI_pick_NVT_move() determines what move will be made next time make_move() is called and selects
// the molecule to which said move will be applied (and creates pointers to its list location).
// Returns the selected move

	int perturb_target = 0; // array index of molecule to be perturbed by next MC move,
	                        // in the array we are about to create.

	std::vector<Molecule *> perturbableMolecules;  // vector of perturbation-eligible molecules
	Molecule  * molecule_ptr       = nullptr;      // linked list traverser
	Molecule  * prev_molecule_ptr  = nullptr;      // molecule linked prev to "current" (in list)
	
	double spinflip_prob        = sys.spinflip_probability;
	double bead_perturb_prob    = sys.bead_perturb_probability;
	double dice_roll_for_move   = Rando::rand();
	double dice_roll_for_target = Rando::rand();


	for( int s=0; s < nSys; s++ ) {

		// populate array with pointers to elements in the linked list		
		for (molecule_ptr = systems[s]->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
			if (!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target))
				perturbableMolecules.push_back(molecule_ptr);
		}

		// select the target (i.e. pick one of the vector elements) & update checkpoint accordingly
		if (perturbableMolecules.size() == 0)
			throw no_molecules_in_system;
		perturb_target = (int)floor(perturbableMolecules.size() * dice_roll_for_target);
		systems[s]->checkpoint->molecule_altered = perturbableMolecules[perturb_target];


		// pick the move
		if (systems[s]->quantum_rotation && (dice_roll_for_move < spinflip_prob)) {
			systems[s]->checkpoint->movetype = MOVETYPE_SPINFLIP;
		}
		else if (dice_roll_for_move < (bead_perturb_prob + spinflip_prob)) {
			systems[s]->checkpoint->movetype = MOVETYPE_PERTURB_BEADS;
		}
		else {
			systems[s]->checkpoint->movetype = MOVETYPE_DISPLACE;
		}

		// Determine the head and tail of the selected molecule, checkpoint->head will be nullptr if molecule is 1st in list
		prev_molecule_ptr = nullptr;
		for (molecule_ptr = systems[s]->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
			if (molecule_ptr == systems[s]->checkpoint->molecule_altered) {
				systems[s]->checkpoint->head = prev_molecule_ptr;
				systems[s]->checkpoint->tail = molecule_ptr->next;
				break;
			}
			prev_molecule_ptr = molecule_ptr;
		}


		// if we have a molecule already backed up (from a previous accept), go ahead and free it
		if (systems[s]->checkpoint->molecule_backup) {
			delete systems[s]->checkpoint->molecule_backup;
			systems[s]->checkpoint->molecule_backup = nullptr;
		}
		// backup the state that will be altered
		systems[s]->checkpoint->molecule_backup = new Molecule(*systems[s]->checkpoint->molecule_altered);

		perturbableMolecules.clear();

	}

	return systems[rank]->checkpoint->movetype;
}




void SimulationControl::PI_make_move(int movetype) {
// apply the move that was selected in the checkpoint

	// update the cavity grid prior to making a move 
	if (sys.cavity_bias)
		for( System *SYS : systems ) {
			SYS->cavity_update_grid();
			SYS->checkpoint->biased_move = 0;
		}


	switch (movetype) {

		case MOVETYPE_DISPLACE:
			PI_displace();
			break;

		case MOVETYPE_SPINFLIP:
			PI_flip_spin();
			break;

		case MOVETYPE_PERTURB_BEADS:
			PI_perturb_beads();
			break;
		/*
			case MOVETYPE_ADIABATIC:
				// change coords of 'altered'
				displace(checkpoint->molecule_altered, pbc, adiabatic_probability, 1.0);
				break;

			case MOVETYPE_VOLUME:
				volume_change(); // I don't want to contribute to the god damned mess -- kmclaugh
				break;
		*/
		default:
			Output::err("MC_MOVES: invalid mc move\n");
			throw invalid_monte_carlo_move;
	}
}




/*
void SimulationControl::PI_insert() {

// int cavities_array_counter = 0;
// System::cavity_t  * cavities_array = nullptr;
// int random_index = 0;
// double com [3] = { 0 };
// double rand[3] = { 0 };
// Molecule  * molecule_ptr = nullptr;
// Atom      * atom_ptr = nullptr;
// Pair      * pair_ptr = nullptr;

// insert a molecule at a random pos and orientation
			// umbrella sampling
			if (cavity_bias && cavities_open) {
				// doing a biased move - this flag lets mc.c know about it
				checkpoint->biased_move = 1;
				// make an array of possible insertion points
				SafeOps::calloc(cavities_array, cavities_open, sizeof(cavity_t), __LINE__, __FILE__);
				for (int i = 0; i < cavity_grid_size; i++) {
					for (int j = 0; j < cavity_grid_size; j++) {
						for (int k = 0; k < cavity_grid_size; k++) {
							if (!cavity_grid[i][j][k].occupancy) {
								for (int p = 0; p < 3; p++)
									cavities_array[cavities_array_counter].pos[p] = cavity_grid[i][j][k].pos[p];
								++cavities_array_counter;
							}
						} // end k
					} // end j
				} // end i
				// insert randomly at one of the free cavity points
				random_index = (cavities_open - 1) - (int)rint(((double)(cavities_open - 1))*get_rand());
				for (int p = 0; p < 3; p++)
					com[p] = cavities_array[random_index].pos[p];
				// free the insertion array
				SafeOps::free(cavities_array);
			} // end umbrella

			else {
				// insert the molecule to a random location within the unit cell
				for (int p = 0; p < 3; p++)
					rand[p] = 0.5 - get_rand();
				for (int p = 0; p < 3; p++) {
					com[p] = 0;
					for (int q = 0; q < 3; q++)
						com[p] += pbc.basis[q][p] * rand[q];
				}
			}

			// process the inserted molecule
			for (atom_ptr = checkpoint->molecule_backup->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
				// move the molecule back to the origin and then assign it to com
				for (int p = 0; p < 3; p++)
					atom_ptr->pos[p] += com[p] - checkpoint->molecule_backup->com[p];
			}

			// update the molecular com
			for (int p = 0; p < 3; p++)
				checkpoint->molecule_backup->com[p] = com[p];
			// give it a random orientation
			checkpoint->molecule_backup->rotate_rand_pbc(1.0); // , pbc, &mt_rand );

			// insert into the list
			if (num_insertion_molecules) {
				// If inserting a molecule from an insertion list, we will always insert at the end
				checkpoint->head->next = checkpoint->molecule_backup;
				checkpoint->molecule_backup->next = nullptr;
			}
			else {
				if (!checkpoint->head) {      // if we're at the start of the list:
					molecules = checkpoint->molecule_backup;
				}
				else {
					checkpoint->head->next = checkpoint->molecule_backup;
				}
				checkpoint->molecule_backup->next = checkpoint->molecule_altered;
			}

			// set new altered and tail to reflect the insertion
			checkpoint->molecule_altered = checkpoint->molecule_backup;
			checkpoint->tail = checkpoint->molecule_altered->next;
			checkpoint->molecule_backup = nullptr;

			if (num_insertion_molecules) { //multi sorbate
				// Free all pair memory in the list
				for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
					for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
						pair_ptr = atom_ptr->pairs;
						while (pair_ptr) {
							Pair *temp = pair_ptr;
							pair_ptr = pair_ptr->next;
							SafeOps::free(temp);
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
*/
/*
void SimulationControl::PI_remove() {
	case MOVETYPE_REMOVE:

		// remove a randomly chosen molecule
		if (cavity_bias) {
			if (get_rand() < pow((1.0 - avg_observables->cavity_bias_probability), ((double)cavity_grid_size*cavity_grid_size*cavity_grid_size)))
				checkpoint->biased_move = 0;
			else
				checkpoint->biased_move = 1;
		}

		// remove 'altered' from the list
		if (!checkpoint->head) {	// handle the case where we're removing from the start of the list
			checkpoint->molecule_altered = molecules;
			molecules = molecules->next;
		}
		else {
			checkpoint->head->next = checkpoint->tail;
		}
		//free_molecule( system, system->checkpoint->molecule_altered );
		delete checkpoint->molecule_altered;
		update_pairs_remove();

		//reset atom and molecule id's
		enumerate_particles();

		break;
	}
*/




void SimulationControl::PI_flip_spin() {

	const int nSystems = (int)systems.size();

	for( int s=0; s<nSystems; s++ )
		if( systems[s]->checkpoint->molecule_altered->nuclear_spin == NUCLEAR_SPIN_PARA)
			systems[s]->checkpoint->molecule_altered->nuclear_spin = NUCLEAR_SPIN_ORTHO;
		else
			systems[s]->checkpoint->molecule_altered->nuclear_spin = NUCLEAR_SPIN_PARA;

	return;
}




void SimulationControl::PI_displace() {

	double dice_rolls[6];
	for (int i = 0; i < 6; i++)
		dice_rolls[i] = Rando::rand();
		
	int nSystems = (int) systems.size();
	
	if (sys.rd_anharmonic)
		// displace_1D(checkpoint->molecule_altered, sys.move_factor);
		throw unsupported_setting;

	if (sys.spectre)
		// spectre_displace(checkpoint->molecule_altered, move_factor, spectre_max_charge, spectre_max_target);
		throw unsupported_setting;

	if (sys.gwp) {
		// if (checkpoint->molecule_altered->atoms->gwp_spin) {
		//	   displace(checkpoint->molecule_altered, pbc, gwp_probability, rot_factor);
		//     checkpoint->molecule_altered->displace_gwp(gwp_probability, &mt_rand);
		// }
		// else
		//     displace(checkpoint->molecule_altered, pbc, move_factor, rot_factor);
		throw unsupported_setting; // remove once supported
	}

	Vector3D pi_com (0, 0, 0); // Center-of-Mass of the Path Integral representation of the molecule
	std::vector<Molecule *> altered_molecules;
	for (int s = 0; s < nSystems; s++) {
		// The "random translation/rotation" that is applied needs to be the same for each bead
		altered_molecules.push_back( systems[s]->checkpoint->molecule_altered );
		altered_molecules[s]->update_COM();
		altered_molecules[s]->translate_rand_pbc( sys.move_factor, systems[s]->pbc, dice_rolls );
		Vector3D molecule_com(altered_molecules[s]->com[0], altered_molecules[s]->com[1], altered_molecules[s]->com[2]);
		pi_com += molecule_com;
	}
	pi_com /= nSystems;

	// Create a random rotation for the molecule
	double diceX = Rando::rand_normal();
	double diceY = Rando::rand_normal();
	double diceZ = Rando::rand_normal();
	double dice_angle = Rando::rand() * sys.rot_factor;
	Quaternion rotation(diceX, diceY, diceZ, dice_angle, Quaternion::AXIS_ANGLE_DEGREE);

	
	for( Molecule *altered : altered_molecules ) {

		// move the molecule by -pi_com (positioning it as if the PI bead chain's center of mass was at the origin)
		altered->translate( -pi_com );
		
		// rotate each site about the random axis
		Atom *aPtr = altered->atoms;
		while (aPtr) {
			Vector3D atomic_pos( aPtr->pos[0], aPtr->pos[1], aPtr->pos[2] );
			atomic_pos = rotation.rotate(atomic_pos);
			aPtr->pos[0] = atomic_pos.x();
			aPtr->pos[1] = atomic_pos.y();
			aPtr->pos[2] = atomic_pos.z();
			aPtr = aPtr->next;
		}
		
		// move the molecule by  pi_com
		altered->translate( pi_com );
		altered->update_COM();  // remove any rounding errors from COM post trans/rot
	}

}




void SimulationControl::PI_perturb_beads() {
	
	PI_perturb_beads_orientations();
	PI_perturb_bead_COMs();
	 
}




void SimulationControl::PI_perturb_bead_COMs_ENTIRE_SYSTEM() {

	Molecule *molPtr = nullptr;
	std::vector<Molecule *> molecules_ptr;   // this array of ptrs tracks the current molecule in each system 
	std::vector<Molecule *> altered_backup;   
	molecules_ptr.resize(nSys);
	altered_backup.resize(nSys);

	
	// Back up the molecules targeted for perturbation (if any)
	for (int s = 0; s < nSys; s++) {
		molecules_ptr [s] = systems[s]->molecules; 
		altered_backup[s] = systems[s]->checkpoint->molecule_altered;
	}

	// Step through each perturb-able molecule in lockstep (each image of that molecule in all systems). At each step, the
	// current set of molecules are set as the perturb targets for the group of systems (checkpoint->molecule_altered, for
	// each sys) and then the set is perturbed and we advance to the next perturbable molecule in each system.
	while( molPtr = molecules_ptr[0] ) {
	
		// Check to see if the current molecule should be perturbed
		// Corresponding molecules *should* be the same across all systems, and the entries of molecules_ptr should point to 
		// different images of the same molecule, so only the version of the molecule in system 0 (molPtr) is checked for validity
		if ( ! (molPtr->frozen || molPtr->adiabatic || molPtr->target)) {

			// If the molecule in image/sys 0 looked good, we target this molecule for perturbation in all system images...
			for (int s = 0; s < nSys; s++) {
				if (nullptr == molecules_ptr[s]) {
					// ... assuming corresponding molecules exist, that is.
					Output::err("ERROR: System images are not consistent, head system has more molecules than at least one other system image.\n");
					throw internal_error;
				}
				systems[s]->checkpoint->molecule_altered = molecules_ptr[s];
			}
			// Perturb all the beads for the targeted molecule in every image
			PI_perturb_bead_COMs(nSys);
		}
		
		// advance our collection of molecule pointers to the next molecule in their respective lists
		for( int s=0; s<nSys; s++ )
			molecules_ptr[s] = molecules_ptr[s]->next;
	}

	// Restore the targeted molecules 
	for (int s = 0; s < nSys; s++) 
		systems[s]->checkpoint->molecule_altered = altered_backup[s];


}
void SimulationControl::PI_perturb_bead_COMs() {
	PI_perturb_bead_COMs( PI_trial_chain_length ); // number of beads in a "trial chain" -- the number beads to move when generating trial COM configs)
}
void SimulationControl::PI_perturb_bead_COMs(int n) {
// n is the number of beads to move

// The algorithm for center-of-mass bead perturbation methods comes from:
// Coker et al.;  J.Chem.Phys. 86, 5689 (1987); doi: 10.1063/1.452495

// General algorithm:
// 1. Isolate a segment of the PI chain representation for a single molecule. Every molecule will move, but only these
//    will move relative to the rest of the collection.
// 2. We make a copy of the COM coords for every member image of the quantum object and place them all in beads[]
//    We will make various changes to the members of beads[] and apply those changes to the simulator's copies of those molecules at the end.
// 3. We compute the COM of the images themselves (in an average sense, this is the point location of the quantum object): chain_COM
//    If we compute the COM of every atom in every molecule of every image (for THIS 1 molecule) the COM should be the same as chain_COM
// 4. We generate perturbations that we apply to the isolated beads/images.
// 5. We compute the change this would have on the COM of the entire, aggregate quantum object (not just the isolated images): delta_COM 
// 6. We subtract delta_COM from every bead[] such that when the perturbations are applied to the isolated images, the COM of the 
//    targeted quantum object will not be moved. Thus, the molecule cannot drift due to bead perturbations.
// 7. Finally, we apply our perturbations and corrections to the actual molecules in our system. Note that in each System[] this change
//    will affect exactly 1 molecule and recall that together these changes will be the perturbation that is applied to 1 object which is
//    represented by a collection of copies of said object, each in a slightly different position/orientation. The average of these different
//    configurations is the potential of the quantum object.


	double beta = 1.0/(kB*sys.temperature);

	static int starterBead = 0; // This bead will not move, and is the index of the first bead to act as an "initial bead". This function will
	                            // construct shorter and shorter chains of "trial chains" with an anchor bead at the start of each chain, this
	                            // this  is the index of the first bead of the first chain. In lang of the Coker ref, starterBead is r[1] in
	                            // the density matrix K(r[1], r[p+1]; beta). We advance starterBead each time we hit this fxn so that the 
	                            // "anchor beads" are different each time, and the effected area of our perturbations has a tendency to rotate
	                            // around the PI loop. The indices of all other participating beads are computed relative to starterBead.
	
	int prevBead_idx  = starterBead; //  index of the bead that comes directly before the bead that is currently being moved  in the 
	                                 //  PI chain (this bead will remain stationary).
	int bead_idx      = (prevBead_idx+1)   % nSys;   //  index of the bead that is currently being moved
	int finalBead_idx = (prevBead_idx+n+1) % nSys;   //  final bead of "trial chain"
	                                                 //  (this bead is also stationary, and may be the same as prevBead)
		 
	
	// Advance the index for the starter bead to the next in the chain. This will be the starter
	// next time, so that every bead gets to start before any bead can start twice. 
	starterBead = (starterBead + 1 ) % nSys; 
	
	
	double Mass = AMU2KG * systems[0]->checkpoint->molecule_altered->mass; // mass of the molecule whose bead configuration is being perturbed
	

	
	
	// populate a vector of Vecs with the COM data from the selected molecule, and compute the center-of-mass (COM) for the PI bead chain
	std::vector<Vector3D> beads;
	Vector3D chain_COM(0, 0, 0);
	for (int s = 0; s < nSys; s++) {
		
		systems[s]->checkpoint->molecule_altered->update_COM();

		double x = systems[s]->checkpoint->molecule_altered->com[0];
		double y = systems[s]->checkpoint->molecule_altered->com[1];
		double z = systems[s]->checkpoint->molecule_altered->com[2];
		Vector3D image_COM(x, y, z); // center of mass of the  target molecule for bead s,
		                             // which is only 1/nSys of the aggregate PI rep of this molecule
		beads.push_back( image_COM );
		chain_COM += image_COM;
	}
	chain_COM /= nSys;

	// compute the perturbation for bead COM positions
	
						 // n is the qty of beads that will be perturbed
	double tB = n;       // this corresponds to all non-cancelling factors of t[ i ] in reference
	double tA = 1.0 + n; // this corresponds to t[i-1] in the same
	                     // the reference has more factors, but they all cancel such that tB/tA is all that remains.

	for( int j=1; j<=n; j++ ) 
	{
		double init_factor  = tB-- / tA--;          // t[i] / t[i-1]    Eqs. 3.9 & 3.10  (seq is e.g. 4/5 -> 3/4 -> 2/3 -> 1/2)
		double term_factor  = 1.0 - init_factor;    // tau  / t[i-1]    Eq.  3.9
		double sigma_factor = sqrt(  (hBar2*beta*init_factor) / (nSys*Mass)  ) * METER2ANGSTROM;  // Eqs. 3.10/3.12 + conv m -> Angstroms
		
		// create the Vec(tor) along which the target bead will be perturbed out of its 'average' position
		Vector3D perturbation( Rando::rand_normal(), Rando::rand_normal(), Rando::rand_normal() ); 

		//                 |---- weighted average position on line connecting existing beads ----|     |------ perturbation ------|
		beads[bead_idx]  =  (init_factor*beads[prevBead_idx]) + (term_factor*beads[finalBead_idx])  +  (sigma_factor*perturbation); // Eq. 3.12 
		

		prevBead_idx = (prevBead_idx + 1) % nSys;  // advance prevBead index
		bead_idx     = (prevBead_idx + 1) % nSys;  // current bead index will be 1 greater than prevBead
	}

	// Compute the center of mass for the chain, post-perturbation
	Vector3D delta_COM(0, 0, 0);
	for( Vector3D bead : beads )
		delta_COM += bead;
	delta_COM = (delta_COM/nSys) - chain_COM;
	

	// Shift the COM coord of the individual beads, such that the COM of the entire chain will not be changed 
	for( Vector3D bead : beads ) 
		bead -= delta_COM;
	
	// now impose the perturbations we've computed back onto the actual system representations
	for (int s = 0; s < nSys; s++)
		systems[s]->checkpoint->molecule_altered->move_to_(beads[s].x(), beads[s].y(), beads[s].z());
}




void SimulationControl::PI_perturb_beads_orientations() {
	char * moleculeID       = systems[0]->checkpoint->molecule_altered->moleculetype;
	int    orientation_site = SimulationControl::get_orientation_site(moleculeID);
	double bond_length      = SimulationControl::get_bond_length(moleculeID);

	// Exit early if requisite orientation data is not present
	if (  (orientation_site < 0)  ||  (bond_length <= 0)  )
		return;
	
	generate_orientation_configs();
	apply_orientation_configs();
}




void SimulationControl::generate_orientation_configs() {
	
	double sorbate_reduced_mass = SimulationControl::get_reduced_mass(systems[0]->checkpoint->molecule_altered->moleculetype);
	if (sorbate_reduced_mass < 0) {
		char buffer[maxLine];
		sprintf(buffer, "No reduced mass specified for moveable/sorbate molecule \"%s\"\n", systems[0]->checkpoint->molecule_altered->moleculetype);
		Output::err(buffer);
		throw missing_required_datum;
	}
	
	double sorbate_bond_length = SimulationControl::get_bond_length(systems[0]->checkpoint->molecule_altered->moleculetype);
	if (sorbate_bond_length < 0) {
		char buffer[maxLine];
		sprintf(buffer, "No bond length specified for moveable/sorbate molecule \"%s\"\n", systems[0]->checkpoint->molecule_altered->moleculetype);
		Output::err(buffer);
		throw missing_required_datum;
	}
	sorbate_bond_length /= METER2ANGSTROM;
	double b2 = sorbate_bond_length * sorbate_bond_length;

	double u_kB_T = sorbate_reduced_mass * kB * sys.temperature;
	orientations[0].randomize();
	generate_orientation_configs(0, (unsigned int)orientations.size(), 2, (unsigned int)orientations.size(), b2, u_kB_T );
}
void SimulationControl::generate_orientation_configs(unsigned int start, unsigned int end, unsigned int p, unsigned int numBeads, double b2, double ukT ) {
// This algorithm comes from:
// Subramanian et al.;  J. Chem. Phys. 146, 094105 (2017); doi: 10.1063/1.4977597]

	const double two_PI = 2.0 * pi;
	

	if (p <= numBeads) {

		// We are given two vectors, I and K, and we want to place J.   J will be halfway between I and K on the polymer/bead
		// chain. We start with a "spring constant" that varies in a wide range, but is consistent with the orientation range 
		// of the entire quantum particle. At the end of the process the spring constant used for placement exactly describes
		// the distribution, and I & K will be directly adjacent beads (to J) on the bead chain. 
		unsigned int J_idx = (start + end) / 2;
		unsigned int K_idx = (end == numBeads) ? 0 : end;
		Vector3D vec_I = orientations[start];
		Vector3D vec_K = orientations[K_idx];



		// "bisector" is the vector that is the exact middle position between the input vectors, I & K. It will be the starting
		// point for the perturbation---normalized since we intend to use it as an axis of rotation.
		Vector3D bisector((vec_I + vec_K) / 2.0);
		bisector.normalize();


		// vec_IK is a vector that is orthogonal to bisector. In most cases, the vector extending from I to K will be 
		// orthogonal already. However, in the initial case, I and K are the same vector, so we create a vector that is
		// only special in the fact that it is different from bisector. We take a cross product with this vector to find  
		// a vector that is orthogonal to bisector.
		Vector3D vec_IK;
		double psi_IK = 0; // this is the angle (in radians) between the vectors I & K

		if (p > 2) {
			vec_IK = vec_K - vec_I;
			psi_IK = Vector3D::angle(vec_I, vec_K);
		}
		else {
			// In this case, the angle between I & K is 0, so we leave psi_IK as-is.
			vec_IK.set(1, 2, -3);                       // here vec_IK is just acting as a temp/dummy variable
			Vector3D different_vec = vec_IK + bisector; // different_vec is simply a different vector than bisector
			different_vec.normalize();
			vec_IK = different_vec.cross(bisector);     // (different_vec X bisector) is ortho to bisector
		}


		double C = Rando::rand(); 
		double lambda2 = h*h / (two_PI*ukT);           // lambda^2 for reduced mass of molecule (thermal wavelength^2) 
		double kh = pi * b2 / lambda2;                 // Eq (13b) in ref
		double K = 4.0 * kh * p * cos( psi_IK * 0.5 ); // Eq (17b) in ref

		// Alpha in reference (renamed to be consistent with angle_B)
		const double angle_A = acos(1.0 + (1.0/K) * log( 1.0 - C*(1.0 - exp(-2.0*K)))); // Eq. (18) in ref
		const double angle_B = Rando::rand() * two_PI;  // Beta in reference  (renamed to avoid confusion with 1/kT)


		// We rotate our ortho vector, IK, about bisector by Beta, and this will form the axis about which bisector will
		// be rotated by Alpha in order to arrive at our final orientation vector. The plane through which bisector travels
		// is defined by the vectors bisector and (bisector X vec_Beta), that is, a vector that would define the choice of 
		// Beta (per the reference) is ortho to both bisector and IK. Since vec_Beta was chosen at random from [0..2pi] about 
		// bisector, this means (bisector X vec_Beta) is likewise randomly chosen from [0..2pi] about bisector.
		Quaternion betaRotation(bisector.x(), bisector.y(), bisector.z(), angle_B, Quaternion::AXIS_ANGLE_RADIAN);
		Vector3D vec_Beta = betaRotation.rotate(vec_IK);
		
		// Construct the final rotation axis. The bisector vector (the exact average orientation vector), will be rotated in the 
		// direction of the beta vector by an angle alpha (angle_A), to arrive at the unit vector that will define the orientation
		// of the bead we are placing. 
		Quaternion final_rotation(vec_Beta.x(), vec_Beta.y(), vec_Beta.z(), angle_A, Quaternion::AXIS_ANGLE_RADIAN);
		
		// Finally, move the bisector vector--the average orientation given I and K--toward beta by an angle alpha to realize 
		// create the perturbation. It is already a unit vector, so we record the orientation for use later.
		Vector3D vec_J = final_rotation.rotate(bisector);
		orientations[J_idx].set(vec_J.x(), vec_J.y(), vec_J.z());
		
		// Now that the intermediate orientation has been perturbed and set, we repeat the procedure by placing two more
		// perturbed, intermediate, orientation vectors between I and J, and between J and K
		if (p < numBeads) {
			generate_orientation_configs( start, J_idx, p * 2, numBeads, b2, ukT);
			generate_orientation_configs( J_idx,   end, p * 2, numBeads, b2, ukT);
		}
	}
}




void SimulationControl::apply_orientation_configs() {
// Orient the beads representing a molecule according to a set of vectors that specify each bead's orientation in space
// Traverse the beads for a given molecule and orient each molecule according to its respective orientation vector.

	int nSystems = (int) systems.size();

	// Impose the orientational perturbations we've computed onto the actual system representations
	int orientation_site = SimulationControl::get_orientation_site( systems[0]->checkpoint->molecule_altered->moleculetype );
	if (orientation_site < 0)
		return; // Molecule has no "handle" specified, and so is will not be re-oriented. 
	for (int s = 0; s < nSystems; s++) {
		systems[s]->checkpoint->molecule_altered->orient( orientations[s], orientation_site );
	}
}