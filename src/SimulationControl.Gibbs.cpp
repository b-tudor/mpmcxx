#include <cstring>

#include "SafeOps.h"
#include "SimulationControl.h"
#include "Output.h"

extern uidx rank, size;






bool SimulationControl::check_Gibbs_options() {

	// Right now, Gibbs is going through the normal check_mc_options() validity check.
	// Place checks for any additional constraints here.
	// If constraints are violated:
	// return fail;

	// Checks on probabilities are done in setup_Gibbs_Systems()---if volume_probability is
	// not explicitly set, it is set to 1/N and the system is not set up (to be able to
	// compute N) until setup_Gibbs_Systems()

	return ok;
}





void SimulationControl::initialize_Gibbs_systems() {

	char linebuf[maxLine] = { 0 };

	// THIS WILL SET UP TWO SYSTEMS, systems[0] AND systems[1] IN THE USUAL MPMC METHOD
	System **temp_sys;
	SafeOps::calloc(temp_sys, 2, sizeof(System*), __LINE__, __FILE__);
	temp_sys[0] = new System(sys);
	temp_sys[1] = new System(sys);
	systems.push_back(temp_sys[0]);
	systems.push_back(temp_sys[1]);

	// if there is an alternate PQR file for the second system specified, copy the filename
	// to the normal PQR filename spot, so that the alternate system is generated instead.
	if (strlen(systems[1]->pqr_input_B))
		strcpy(systems[1]->pqr_input, systems[1]->pqr_input_B);

	for (uint32_t i = 0; i<2; i++) {

		systems[i]->setup_simulation_box();

		sprintf(linebuf, "SIM_CONTROL, SYS %d: simulation box configured.\n", i);
		Output::out(linebuf);

		systems[i]->allocateStatisticsMem();

		systems[i]->allocate_pair_lists();
		sprintf(linebuf, "SIM_CONTROL, SYS %d: finished allocating pair lists\n", i);
		Output::out(linebuf);

		systems[i]->pairs();          // get all of the pairwise interactions, exclusions, etc.
		systems[i]->flag_all_pairs(); // set all pairs to initially have their energies calculated 
		sprintf(linebuf, "SIM_CONTROL, SYS %d: finished calculating pairwise interactions\n", i);
		Output::out(linebuf);

		if (systems[i]->cavity_bias)
			systems[i]->setup_cavity_grid();

		// if polarization active, allocate the necessary matrices 
		if (systems[i]->polarization && !systems[i]->cuda && !systems[i]->polar_zodid)
			systems[i]->thole_resize_matrices();

		// if histogram calculation flag is set, allocate grid
		if (systems[i]->calc_hist) {
			systems[i]->setup_histogram();
		}

		if (!systems[i]->use_sg || systems[i]->rd_only)
		{
			sprintf(linebuf, "SIM_CONTROL, SYS %d: Ewald gaussian width = %f A\n", i, systems[i]->ewald_alpha);
			Output::out(linebuf);
			sprintf(linebuf, "SIM_CONTROL, SYS %d: Ewald kmax = %d\n", i, systems[i]->ewald_kmax);
			Output::out(linebuf);
		}

		// ensure that all SPECTRE charges lie within the restricted domain 
		if (systems[i]->spectre)
			systems[i]->spectre_wrapall();
	}


	// Set the probability of the volume move to 1/N if it wasn't explicitly set
	if (sys.volume_probability == 0.0) {
		sys.volume_probability =
			systems[0]->volume_probability =
			systems[1]->volume_probability = 1.0 / (double)(systems[0]->observables->N + systems[1]->observables->N);
	}
	sprintf(linebuf, "SIM_CONTROL: volume change probability is %lf.\n", systems[0]->volume_probability);
	Output::out(linebuf);


	// If quantum rotation is turned off, make spinflip probability 0.0, just in case
	if (sys.quantum_rotation) {
		sprintf(linebuf, "SIM_CONTROL:      spinflip probability is %lf.\n", sys.spinflip_probability);
		Output::out(linebuf);
	}
	else
		sys.spinflip_probability = 0.0;


	if (sys.transfer_probability == 0.0) {
		Output::err("(ERROR) SIM_CONTROL: transfer move probability was either not set, or set to 0.0 in a Gibbs NVT simulation. Set with keyword \"transfer_probability\" in input file.");
		throw missing_setting;
	}
	else {
		sprintf(linebuf, "SIM_CONTROL:      transfer probability is %lf.\n", sys.transfer_probability);
		Output::out(linebuf);
	}

	sprintf(linebuf, "SIM_CONTROL:      displace probability is %lf.\n", 1.0 - (sys.spinflip_probability + sys.transfer_probability + sys.volume_probability));
	Output::out(linebuf);

	// Make sure probability of all moves adds up to 1.0 (displacement is the
	// default move so all the other moves must some to a probability < 1.0
	if (sys.spinflip_probability + sys.volume_probability + sys.transfer_probability >= 1.0) {
		// To allow displacement probabilities of 0.0, change the above check to   ... >    1.0
		Output::err("(ERROR) SIM_CONTROL: Invalid probabilities set. The summed frequencies for spinflip, volume, transfer, and displacement moves may not exceed 1.0.");
	}
}





bool SimulationControl::Gibbs_mc() {

	int     move = 0;
	double  initial_energy[2] = { 0 },
	        final_energy  [2] = { 0 },
	        current_energy[2] = { 0 };
	
	System::mpiData mpi[2];


	for (uint32_t i = 0; i<2; i++) {
		if (systems[i]->sorbateCount > 1)
			SafeOps::calloc(systems[i]->sorbateGlobal, (size_t) systems[i]->sorbateCount, sizeof(System::sorbateAverages_t), __LINE__, __FILE__);
		if (systems[i]->cavity_bias)
			systems[i]->cavity_update_grid(); // update the grid for the first time 
		systems[i]->observables->volume = systems[i]->pbc.volume; // set volume observable
		initial_energy[i] = systems[i]->mc_initial_energy();
		mpi[i] = systems[i]->setup_mpi();

	}

	// save the initial state
	backup_observables_ALL_SYSTEMS();
	move = System::pick_Gibbs_move(systems);

	uint32_t s = 1;
	uint32_t max_step = systems[0]->numsteps;

	// main MC loop 
	for (s = 1; s <= max_step; s++) {

		systems[0]->step = systems[1]->step = s;

		// restore the last accepted energy 
		initial_energy[0] = systems[0]->observables->energy;
		initial_energy[1] = systems[1]->observables->energy;

		// perturb the system
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		System::make_move_Gibbs(systems);


		// calculate the energy change/boltzmann factor
		final_energy[0] = systems[0]->energy();
		final_energy[1] = systems[1]->energy();

#ifdef QM_ROTATION
		// solve for the rotational energy levels 
		if (system->quantum_rotation && (system->checkpoint->movetype == MOVETYPE_SPINFLIP))
			quantum_system_rotational_energies(system);
#endif // QM_ROTATION 

		boltzmann_factor_NVT_Gibbs(*systems[0], initial_energy[0], final_energy[0], *systems[1], initial_energy[1], final_energy[1]);




		// ACCEPT -or- REJECT move for independent MC moves
		// Displacement and Spinflips are handled independently between systems
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if ((move == MOVETYPE_DISPLACE) || (move == MOVETYPE_SPINFLIP))
		{
			for (uint32_t i = 0; i < 2; i++)
			{
				// Metropolis function (performed separately for each of the two systems.
				if (( Rando::rand() < systems[i]->nodestats->boltzmann_factor) && !systems[i]->iterator_failed) {

					// ACCEPT ///////////////////////////////////////////////////////////////
					current_energy[i] = final_energy[i];
					systems[i]->register_accept();

					// Simulated Annealing...
					if (systems[i]->simulated_annealing) {
						if( systems[i]->simulated_annealing_linear ) {
							systems[i]->temperature = systems[i]->temperature + (systems[i]->simulated_annealing_target - systems[i]->temperature) / ((size_t)max_step - s);
							if ((max_step - s) == 0)
								systems[i]->temperature = systems[i]->simulated_annealing_target;
						}
						else {
							systems[i]->temperature = systems[i]->simulated_annealing_target + (systems[i]->temperature - systems[i]->simulated_annealing_target)*systems[i]->simulated_annealing_schedule;
						}
					}
				}
				else {
					// REJECT ///////////////////////////////////////////////////////////////
					current_energy[i] = initial_energy[i]; // used in parallel tempering
					systems[i]->iterator_failed = 0;       // reset the polar iterative failure flag and...
					systems[i]->restore();                 // ...restore from last checkpoint
					systems[i]->register_reject();
				}
			}

		}
		else {

			// ACCEPT -or- REJECT move for coordinated MC moves
			// i.e.:   if(  (move == MOVETYPE_INSERT)   ||   (move == MOVETYPE_REMOVE)   ||   (move == MOVETYPE_VOLUME)  )
			// Transfers (insert/remove) and Volume Shifts are accepted/rejected together depending on the net effect on both systems
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			double boltzmann_factor = systems[0]->nodestats->boltzmann_factor; // BF in both systems will be identical in XFER & VOL moves

																			   // Metropolis function (performed separately for each of the two systems.
			if(( Rando::rand() < boltzmann_factor) && !systems[0]->iterator_failed && !systems[1]->iterator_failed)
			{

				// ACCEPT ///////////////////////////////////////////////////////////////

				for (uint32_t i = 0; i<2; i++) {
					current_energy[i] = final_energy[i];

					// backup observables for this system only
					*systems[i]->checkpoint->observables = *systems[i]->observables;
					systems[i]->register_accept();

					// Simulated Annealing...
					if( systems[i]->simulated_annealing )
					{
						if( systems[i]->simulated_annealing_linear )
						{
							systems[i]->temperature = systems[i]->temperature + (systems[i]->simulated_annealing_target - systems[i]->temperature) / ((size_t)max_step - s);
							if ((max_step - s) == 0)
								systems[i]->temperature = systems[i]->simulated_annealing_target;

						}
						else
							systems[i]->temperature = systems[i]->simulated_annealing_target + (systems[i]->temperature - systems[i]->simulated_annealing_target)*systems[i]->simulated_annealing_schedule;
					}
				}

			}
			else {

				// REJECT ///////////////////////////////////////////////////////////////
				for (uint32_t i = 0; i < 2; i++) {
					current_energy[i] = initial_energy[i]; // used in parallel tempering
					systems[i]->iterator_failed = 0;       // reset the polar iterative failure flag and...
					systems[i]->restore();                 // ...restore from last checkpoint
					systems[i]->register_reject();
				}
			}
		}

		backup_observables_SYS_VECTOR();
		move = System::pick_Gibbs_move(systems);



		// perform parallel_tempering
		if ((systems[0]->parallel_tempering) && ((s % systems[0]->ptemp_freq) == 0))
			systems[0]->temper_system(current_energy[0]);
		if ((systems[1]->parallel_tempering) && ((s % systems[1]->ptemp_freq) == 0))
			systems[1]->temper_system(current_energy[1]);

		// track the acceptance_rate 
		systems[0]->track_ar(systems[0]->nodestats);
		systems[1]->track_ar(systems[1]->nodestats);

		// each node calculates its stats 
		systems[0]->update_nodestats(systems[0]->nodestats, systems[0]->avg_nodestats);
		systems[1]->update_nodestats(systems[1]->nodestats, systems[1]->avg_nodestats);

		// In the event that sys1 and sys0 have different corr times, we will only honor that of sys0
		if (!(s % systems[0]->corrtime) || (s == max_step)) {

			// do this every correlation time, and at the very end
			Output::out("\n");
			for (uint32_t i = 0; i < 2; i++) {

				systems[i]->do_corrtime_bookkeeping(mpi[0]);

				// write the averages to stdout 
				if (systems[i]->fp_histogram)
					systems[i]->write_histogram(systems[i]->fp_histogram, systems[i]->grids->avg_histogram->grid);

				if (!i) {
					if (systems[i]->write_performance(s) < 0) {
						Output::err("MC: could not write performance data to stdout\n");
						throw unknown_file_error;
					}
					Output::out("\n");
				}

				if (systems[i]->display_averages((int)i) < 0) {
					Output::err("MC: could not write statistics to stdout\n");
					throw unknown_file_error;
				}
			}
		}

	} // main loop 



	  // output restart files for each node
	for (uint32_t i = 0; i<2; i++) {
		if (systems[i]->write_molecules_wrapper(systems[i]->pqr_output) < 0) {
			Output::err("MC: could not write final state to disk\n");
			throw unknown_file_error;
		}

		if (!rank) {
			SafeOps::free(mpi[i].rcv_strct);
			SafeOps::free(mpi[i].temperature);
		}

		if (systems[i]->sorbateCount > 1)
			SafeOps::free(mpi[i].sinfo);

		SafeOps::free(mpi[i].snd_strct);
		SafeOps::free(mpi[i].observables);
		SafeOps::free(mpi[i].avg_nodestats);
	}

	return ok;
}




// calculate the boltzmann factor for a Gibbs (NVT) system pair
void SimulationControl::boltzmann_factor_NVT_Gibbs(System &sys1, double initEnergy_1, double finalEnergy_1, System &sys2, double initEnergy_2, double finalEnergy_2) {

	System *sys[2] = { nullptr };
	sys[0] = &sys1;
	sys[1] = &sys2;


	double dEnergy[2] = { 0 };
	dEnergy[0] = finalEnergy_1 - initEnergy_1;
	dEnergy[1] = finalEnergy_2 - initEnergy_2;

	int move = sys[0]->checkpoint->movetype;

	// treat a bad contact on a coordinated move as a reject for both systems
	if (!std::isfinite(finalEnergy_1) || !std::isfinite(finalEnergy_2)) {
		if (move == MOVETYPE_INSERT || move == MOVETYPE_REMOVE || move == MOVETYPE_VOLUME) {
			sys[0]->observables->energy = MAXVALUE;
			sys[0]->nodestats->boltzmann_factor = 0;
			sys[1]->observables->energy = MAXVALUE;
			sys[1]->nodestats->boltzmann_factor = 0;
		}
		return;
	}




	// GIBBS DISPLACEMENT MOVE
	////////////////////////////////////////////////////////////////////////////////////

	if ((sys[0]->checkpoint->movetype == MOVETYPE_DISPLACE) || (sys[1]->checkpoint->movetype == MOVETYPE_DISPLACE)) {

		// If both moves are not displacement moves, then the internal state has been corrupted
		if (sys[0]->checkpoint->movetype != sys[1]->checkpoint->movetype)
			throw internal_err_invalid_mc_move;

		// SYSTEM 1 ///////////////////////////////////////
		if (!std::isfinite(finalEnergy_1)) {
			sys[0]->observables->energy = MAXVALUE;
			sys[0]->nodestats->boltzmann_factor = 0;
		}
		else
			sys[0]->nodestats->boltzmann_factor = exp(-dEnergy[0] / sys[0]->temperature);

		// SYSTEM 2 ///////////////////////////////////////
		if (!std::isfinite(finalEnergy_2)) {
			sys[1]->observables->energy = MAXVALUE;
			sys[1]->nodestats->boltzmann_factor = 0;
		}
		else
			sys[1]->nodestats->boltzmann_factor = exp(-dEnergy[1] / sys[1]->temperature);

		return;
	}




	// GIBBS TRANSFER MOVE
	// if (sys[0]->checkpoint->movetype == MOVETYPE_INSERT || sys[0]->checkpoint->movetype == MOVETYPE_REMOVE) {
	if ((sys[0]->checkpoint->movetype == MOVETYPE_INSERT && sys[1]->checkpoint->movetype == MOVETYPE_REMOVE) || (sys[0]->checkpoint->movetype == MOVETYPE_REMOVE && sys[1]->checkpoint->movetype == MOVETYPE_INSERT)) {

		// transfer from box A to box B
		int A, B;
		if (sys[0]->checkpoint->movetype == MOVETYPE_REMOVE) {
			A = 0;
			B = 1;
		}
		else {
			A = 1;
			B = 0;
		}

		double V_A   = sys[A]->observables->volume;
		double N_A  = sys[A]->observables->N;
		double V_B   = sys[B]->observables->volume;
		double N_B   = sys[B]->observables->N;
		double Beta = 1.0 / sys[0]->temperature;

		sys[0]->nodestats->boltzmann_factor =
			sys[1]->nodestats->boltzmann_factor = (N_A / V_A)*(V_B / (N_B + 1)) * exp(-Beta * dEnergy[A] - Beta * dEnergy[B]);

		return;
	}




	// GIBBS VOLUME SHIFT
	if (sys[0]->checkpoint->movetype == MOVETYPE_VOLUME)
	{

		// if both systems are not set to MOVETYPE_VOLUME, then the internal state has been corrupted.
		if (sys[0]->checkpoint->movetype != sys[1]->checkpoint->movetype)
			throw internal_err_invalid_mc_move;


		double     dV = 0;      // Volume  change of system 
		double   beta = 1.0 / sys[0]->temperature;
		double   V[2] = { 0 }; // Current volume of system A/B
		double   N[2] = { 0 }; // Current particle count, system A/B



		for (int i = 0; i < 2; i++) {
			V[i] = sys[i]->observables->volume;
			N[i] = sys[i]->observables->N;
		}
		dV = sys[0]->checkpoint->observables->volume - V[0];
		sys[0]->nodestats->boltzmann_factor =
			sys[1]->nodestats->boltzmann_factor = pow((V[0] + dV) / V[0], N[0])  *  pow((V[1] - dV) / V[1], N[1])   *  exp(-beta * dEnergy[0] - beta * dEnergy[1]);

		return;
	}




	// GIBBS SPINFLIP MOVE
	if (sys[0]->checkpoint->movetype == MOVETYPE_SPINFLIP)
	{
		double g = 0, u = 0, partfunc_ratio = 0;

		// if both systems are not set to SPINFLIP, then the internal state has been corrupted.
		if (sys[0]->checkpoint->movetype != sys[1]->checkpoint->movetype)
			throw internal_err_invalid_mc_move;

		// SYSTEM 1 ////////////////////////////////////////////////////////////////
		if( !std::isfinite( finalEnergy_1 )) {
			sys[0]->observables->energy = MAXVALUE;
			sys[0]->nodestats->boltzmann_factor = 0;
		}
		else {
			g = sys[0]->checkpoint->molecule_altered->rot_partfunc_g;
			u = sys[0]->checkpoint->molecule_altered->rot_partfunc_u;
			if (sys[0]->checkpoint->molecule_altered->nuclear_spin == NUCLEAR_SPIN_PARA)
				partfunc_ratio = g / (g + u);
			else
				partfunc_ratio = u / (g + u);
			// set the boltz factor, including ratio of partfuncs for different symmetry rotational levels
			sys[0]->nodestats->boltzmann_factor = partfunc_ratio;
		}


		// SYSTEM 2 ////////////////////////////////////////////////////////////////
		if (!std::isfinite(finalEnergy_2)) {
			sys[1]->observables->energy = MAXVALUE;
			sys[1]->nodestats->boltzmann_factor = 0;
		}
		else {
			g = sys[1]->checkpoint->molecule_altered->rot_partfunc_g;
			u = sys[1]->checkpoint->molecule_altered->rot_partfunc_u;
			if (sys[1]->checkpoint->molecule_altered->nuclear_spin == NUCLEAR_SPIN_PARA)
				partfunc_ratio = g / (g + u);
			else
				partfunc_ratio = u / (g + u);
			// set the boltz factor, including ratio of partfuncs for different symmetry rotational levels
			sys[1]->nodestats->boltzmann_factor = partfunc_ratio;
		}

		return;
	}


	// if we get here, then an invalid Monte Carlo move was selected
	throw invalid_monte_carlo_move;
}