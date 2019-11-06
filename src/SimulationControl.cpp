#include "SimulationControl.h"
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <time.h>

#include "Fugacity.h"
#include "Output.h"
#include "PeriodicBoundary.h"
#include "SafeOps.h"

#ifdef _MPI
	#include <mpi.h>
#endif
extern int size, rank;
extern bool mpi;



std::map<std::string, int> SimulationControl::sorbate_data_index;
std::vector<SimulationControl::molecular_metadata> SimulationControl::sorbate_data;


SimulationControl::~SimulationControl() {

	for( unsigned int s=0; s<systems.size(); s++ ) {
		free(systems[s]);
		systems[s] = nullptr;
	}

	if( system_energies )
		free( system_energies );
}
SimulationControl::SimulationControl(char *inFilename, bool rAR, bool writeFrames) : report_AR(rAR), write_PI_frames(writeFrames)
{
	char linebuf[maxLine];
	PI_nBeads = 0;
	PI_trial_chain_length = 0;
	

	sprintf(linebuf, "SIM_CONTROL: running parameters found in: %s\n", inFilename);
	Output::out1(linebuf);
 	read_config( inFilename );  // Parse input file
	Output::out1( "SIM_CONTROL: Finished reading config file.\n" );

	if (mpi) {
		if (PI_nBeads) {
			if (PI_nBeads != size) {
				Output::out1("SIM_CONTROL: When computing path integrals with MPI, the Trotter number is set\n");
				Output::out1("             by the MPI process count, not the input file. If the 'trotter_number'\n");
				Output::out1("             option is set, it must match the MPI world size.\n"); 
				throw invalid_setting;
			}
		}
		PI_nBeads = size;
	}

	if( check_system() )             // Check system for compatible parameters.
		Output::out1("SIM_CONTROL: input file validated.\n");
	else {
		Output::err("SIM_CONTROL: input file has failed validation.\n");
		throw invalid_input;
	}

	orientations.resize(size); // allocate space for one orientation vector per bead

}




void SimulationControl::initializeSimulationObjects() {
// Initialize the simulation object (or objects, depending on sim parameters), after which,
// the simulation should be ready to run. This includes allocation of requisite memory,
// initialization of data structures, calculation of certain computed constants, etc.

	char linebuf[maxLine] = {0};

	//  Seed the global Random Number Generator (RNG)  ////////////////////////////////////////////

	uint32_t seed;
	if( sys.preset_seed_on ) 
		seed = sys.preset_seed;
	else {
		seed = (unsigned int) time(0);
		if( mpi )
			if (sys.ensemble == ENSEMBLE_PATH_INTEGRAL_NVT) {
				#ifdef _MPI
					// If using MPI for path integrals, every MPI instance must use the same seed. 
					MPI_Bcast( (void *) &seed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
				#endif
			}
	}
	Rando::seed(seed);

	if (mpi) {

		#ifdef _MPI
			for (int i = 0; i < size; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				if (rank == i) {
					sprintf(linebuf, "SIM_CONTROL->SYSTEM[ %d ]: RNG initialized. Seed = %u\n", rank, seed);
					Output::out(linebuf);
				}
			}
		#endif

	} else {
		sprintf(linebuf, "SIM_CONTROL: RNG initialized. Seed = %u\n", seed);
		Output::out1(linebuf);
	}
	



	//  Allocate and populate required data structures.  //////////////////////////////////////////

	switch( sys.ensemble ) {

		case ENSEMBLE_PATH_INTEGRAL_NVT:
			initialize_PI_NVT_Systems();
			return;
		
		case ENSEMBLE_NVT_GIBBS:
			initialize_Gibbs_systems();
			return;

		default:
 
			//  This is the normal MPMC setup (i.e. non-Gibbs, non-Path Integral)  ////////////////

			// set up the simulation box: pbc and read in molecules
			sys.setup_simulation_box();
			Output::out("SIM_CONTROL: simulation box configured.\n");

			sys.allocateStatisticsMem();

			sys.allocate_pair_lists();
			Output::out( "SIM_CONTROL: finished allocating pair lists\n" );
			sys.pairs();          // get all of the pairwise interactions, exclusions, etc.
			sys.flag_all_pairs(); // set all pairs to initially have their energies calculated 
			Output::out( "SIM_CONTROL: finished calculating pairwise interactions\n" );

			if( sys.cavity_bias ) 
				sys.setup_cavity_grid();

			// if polarization active, allocate the necessary matrices 
			if(  sys.polarization   &&   !sys.cuda   &&   !sys.polar_zodid  )
				sys.thole_resize_matrices();

			// if histogram calculation flag is set, allocate grid
			if( sys.calc_hist ) {
				sys.setup_histogram();
			}

			// seed the rng if neccessary 
			if( sys.ensemble != ENSEMBLE_TE   &&   sys.ensemble != ENSEMBLE_REPLAY ) {
				if (sys.preset_seed_on)
					sys.mt_rand.seed( sys.preset_seed );
				else 
					sys.mt_rand.seed((unsigned int)time(0));
				
			}

			if(  ! sys.use_sg   ||   sys.rd_only  ) 
			{
				sprintf(linebuf, "SIM_CONTROL: Ewald gaussian width = %f A\n", sys.ewald_alpha);
				Output::out(linebuf);
				sprintf(linebuf, "SIM_CONTROL: Ewald kmax = %d\n", sys.ewald_kmax);
				Output::out(linebuf);
			}

			#ifdef OPENCL
				if( sys.opencl) {
					Output::out("SIM_CONTROL: Installing OpenCL kernel...\n");
        				sys.ocl = setup_ocl();
					Output::out("SIM_CONTROL: ...OpenCL kernel installed.\n");
				}
			#endif
	
			// ensure that all SPECTRE charges lie within the restricted domain 
			if( sys.spectre ) 
				sys.spectre_wrapall();

			return;
	}
}




void SimulationControl::read_config(char *inFilename) {
// Parses a simulation input file and populates the System Controller with
// the data found therein. Data read in is largely unvalidated.

	char linebuffer[maxLine], *n;
	char errormsg[maxLine];
	char token[maxTokens][maxLine] = { '\0' };
	FILE *fp;
	int i, linenum;
	
	// open the config file or error 
	fp = SafeOps::openFile(inFilename, "r", __LINE__, __FILE__ );
	
	// loop over each line 
	std::memset( linebuffer, 0, maxLine );
	n = fgets( linebuffer, maxLine, fp);

	linenum = 0;
	while (n) {

		linenum++;
		// grab a line and parse it out 
		for( i=0; i<maxTokens; i++)
			std::memset(token[i], 0, maxLine); //clear a token
		int ck = scanf(linebuffer, "%s %s %s %s %s %s %s %s %s %s",
			token[0], token[1], token[2], token[3], token[4],
			token[5], token[6], token[7], token[8], token[9]);

		//parse and apply a command
		if( !process_command( token )) {
			sprintf(errormsg, "SIM_CONTROL: invalid input on line %d.\n", linenum);
			Output::err(errormsg);
			sprintf(errormsg, "> %s\n", linebuffer);
			// Change a double \n\n ending (common here) to a single \n ending.
			if (errormsg[strlen(errormsg)-2] == '\n')
				errormsg[strlen(errormsg)-1] = 0;
			Output::err(errormsg);
			throw invalid_input;
		}
		
		std::memset(linebuffer, 0, maxLine);
		n = fgets(linebuffer, maxLine, fp);

	}

	// close the config file 
	fclose(fp);

	return;
}




bool SimulationControl::process_command( char token[maxTokens][maxLine] ) {
// Read each input line from the input file and set individual system flags accordingly.
// Minimal error checking performed, mostly in the form of checking for malformed
// command sequences. 


	// check for comment/blanks
	if (!strncmp(token[0], "!", 1)) return ok;
	if (!strncmp(token[0], "#", 1)) return ok;
	if (!strcmp (token[0], ""    )) return ok;


	// set job name (CRC)
	if( SafeOps::iequals(token[0], "job_name") ) {
		strcpy( sys.job_name, token[1] );
		return ok;
	}

	
	// ensemble options
	if( SafeOps::iequals(token[0], "ensemble") ) 
	{
		if( SafeOps::iequals(token[1], "nvt") )
			sys.ensemble = ENSEMBLE_NVT;
		else if( SafeOps::iequals(token[1], "uvt") )
			sys.ensemble = ENSEMBLE_UVT;
		else if( SafeOps::iequals(token[1], "surf") )
			sys.ensemble = ENSEMBLE_SURF;
		else if( SafeOps::iequals(token[1], "surf_fit") )
			sys.ensemble = ENSEMBLE_SURF_FIT;
		else if( SafeOps::iequals(token[1], "nve") )
			sys.ensemble = ENSEMBLE_NVE;
		else if( SafeOps::iequals(token[1], "total_energy") )
			sys.ensemble = ENSEMBLE_TE;
		else if( SafeOps::iequals(token[1], "npt") )
			sys.ensemble = ENSEMBLE_NPT;
		else if( SafeOps::iequals(token[1], "replay") )
			sys.ensemble = ENSEMBLE_REPLAY;
		else if( SafeOps::iequals(token[1], "pi_nvt" ) )
			sys.ensemble = ENSEMBLE_PATH_INTEGRAL_NVT;
		else if( SafeOps::iequals(token[1], "nvt_gibbs" ) )
			sys.ensemble = ENSEMBLE_NVT_GIBBS;
		else
			return fail;

		return ok;
	}

	if (SafeOps::iequals(token[0], "sorbate_orientation_site")) {
		int index;
		if ((strlen(token[1]) <= 0) && (strlen(token[2]) <= 0)) {
			Output::err("ERROR: incorrect syntax for sorbate_orientation_site command. orientation_site <Molecule Type ID> <Site#, such that the site that appears first is '0'>\nE.g.:\n  sorbate_orientation_site H2 1\n");
			return fail;
		}
		SafeOps::atoi(token[2], index);
		
		add_orientation_site_entry(token[1], index);
		return ok;
	}

	if (SafeOps::iequals(token[0], "sorbate_bondlength")) {
		double bondLength;
		if ((strlen(token[1]) <= 0) && (strlen(token[2]) <= 0)) {
			Output::err("ERROR: incorrect syntax for sorbate_bondlength command. sorbate_bondlength <Molecule Type ID> <bond length, in Angstroms>\nE.g.:\n  sorbate_bondlength H2 0.742\n");
			return fail;
		}
		SafeOps::atod(token[2], bondLength );

		add_bond_length_entry(token[1], bondLength);
		return ok;
	}
	if (SafeOps::iequals(token[0], "sorbate_reducedMass")) {
		double reducedMass;
		if ((strlen(token[1]) <= 0) && (strlen(token[2]) <= 0)) {
			Output::err("ERROR: incorrect syntax for sorbate_reducedMass command. sorbate_reducedMass <Molecule Type ID> <reduced mass in kg>\nE.g.:\n  sorbate_reducedMass H2 8.368618e-28\n");
			return fail;
		}
		SafeOps::atod(token[2], reducedMass);

		add_reduced_mass_entry(token[1], reducedMass);
		return ok;
	}
	
		

	// random seed options
	if( SafeOps::iequals(token[0], "preset_seed")  || SafeOps::iequals(token[0], "seed")) {
		if(  !SafeOps::atou(token[1], sys.preset_seed )  ) return fail;
		sys.preset_seed_on = 1;
		return ok;
	}


	// Is this a restart of a parallel job?
	////////////////////////////////////////
	if( SafeOps::iequals(token[0], "parallel_restarts") ) {

		if( SafeOps::iequals(token[1], "on") ) {
			sys.parallel_restarts = 1;
			Output::out1("SIM_CONTROL: parallel restart option selected.\n");
		}
		else if( SafeOps::iequals(token[1], "off") ) {
			sys.parallel_restarts = 0;
		}
		else
			return fail;
		return ok;
	}

	// Is this a path integral run on a single thread? 
	//////////////////////////////////////////////////

	if (SafeOps::iequals(token[0], "trotter_number")) {
		if (!SafeOps::atoi(token[1], PI_nBeads ))
			return fail;
		return ok;
	}


	// surf options
	////////////////////////


	// Option for fitting against arbitrary configurations, VS the default behavior of fitting
	// against a small set of orientations, while only varying their separation distance.
	if( SafeOps::iequals(token[0], "fit_arbitrary_configs") ) {

		if( SafeOps::iequals(token[1], "on") ) 
			sys.surf_fit_arbitrary_configs = 1;
		else if( SafeOps::iequals(token[1], "off") ) 
			sys.surf_fit_arbitrary_configs = 0;
		else
			return fail;
		return ok;
	}

	if( SafeOps::iequals(token[0], "surf_decomp") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.surf_decomp = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.surf_decomp = 0;
		else 
			return fail; //unrecognized argument
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_min") )	{
		if( !SafeOps::atod(token[1], sys.surf_min) ) 
			return fail;
		return ok;
	} 
	if( SafeOps::iequals(token[0], "surf_max") ) {
		if( !SafeOps::atod(token[1], sys.surf_max) ) 
			return fail;
		return ok;
	} 
	if( SafeOps::iequals(token[0], "surf_inc") ) {
		if( !SafeOps::atod(token[1], sys.surf_inc) )
			return fail;
		return ok;
	} 
	if( SafeOps::iequals(token[0], "surf_ang") ) {
		if( !SafeOps::atod(token[1], sys.surf_ang) )
			return fail;
		return ok;
	} 
	if( SafeOps::iequals(token[0], "surf_print_level") ) {
		if( !SafeOps::atoi(token[1], sys.surf_print_level) )
			return fail;
		return ok;
	}
	//allows us to specify the surf-fit scales in the input file
	if( SafeOps::iequals(token[0], "surf_weight_constant") ) {
		if (!SafeOps::atod(token[1], sys.surf_weight_constant))
			return fail;
		sys.surf_weight_constant_on = 1;
		return ok;
	} 
	if( SafeOps::iequals(token[0], "surf_scale_q") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_q) )
			return fail;
		sys.surf_scale_q_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_r") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_r) )
			return fail;
		sys.surf_scale_r_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_epsilon") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_epsilon) )
			return fail; 
		sys.surf_scale_epsilon_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_sigma") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_sigma) )
			return fail; 
		sys.surf_scale_sigma_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_omega") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_omega) ) 
			return fail;
		sys.surf_scale_omega_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_alpha") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_alpha) ) 
			return fail;
		sys.surf_scale_alpha_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_pol") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_pol) ) 
			return fail;
		sys.surf_scale_pol_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_c6") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_c6) ) 
			return fail;
		sys.surf_scale_c6_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_c8") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_c8) )
			return fail;
		sys.surf_scale_c8_on = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_scale_c10") ) {
		if( !SafeOps::atod(token[1], sys.surf_scale_c10) ) 
			return fail;
		sys.surf_scale_c10_on = 1;
		return ok;
	}
	if (SafeOps::iequals(token[0], "surf_qshift")) {
		if (SafeOps::iequals(token[1], "on")) {
			sys.surf_qshift_on = 1;
			Output::out1("SIM_CONTROL: surf_qshift is ON.\n");
			Output::out1("SIM_CONTROL: only use qshift with x-axis aligned linear molecules.\n");
		}
		else if (SafeOps::iequals(token[1], "off")) {
			sys.surf_qshift_on = 0;
			Output::out1("SIM_CONTROL: surf_qshift is OFF.\n");
		}
		else
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_preserve") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.surf_preserve = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.surf_preserve = 0;
		else
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_preserve_rotation") ) {
		if( sys.surf_preserve_rotation_on ) {
			Output::err("SIM_CONTROL: surf_preserve_rotationalready set.\n");
			return fail;
		}
		sys.surf_preserve_rotation_on = 1;
		if( !SafeOps::atod(token[1], sys.surf_preserve_rotation_alpha1) ) return fail; 
		if( !SafeOps::atod(token[2], sys.surf_preserve_rotation_beta1)  ) return fail; 
		if( !SafeOps::atod(token[3], sys.surf_preserve_rotation_gamma1) ) return fail; 
		if( !SafeOps::atod(token[4], sys.surf_preserve_rotation_alpha2) ) return fail; 
		if( !SafeOps::atod(token[5], sys.surf_preserve_rotation_beta2)  ) return fail; 
		if( !SafeOps::atod(token[6], sys.surf_preserve_rotation_gamma2) ) return fail; 
		return ok;
	}
	if (SafeOps::iequals(token[0], "surf_global_axis") ) {
		if( SafeOps::iequals(token[1], "on") ) {
			sys.surf_global_axis_on = 1;
			Output::out1("SIM_CONTROL: surf_global_axis is ON.\n");
		}
		else if( SafeOps::iequals(token[1], "off") ) {
			sys.surf_global_axis_on = 0;
			Output::out1("SIM_CONTROL: surf_global_axis is OFF.\n");
		}
		else 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_descent") ) {
		if( SafeOps::iequals(token[1], "on") ) {
			sys.surf_descent = 1;
			Output::out1("SIM_CONTROL: surf_descent is ON.\n");
		}
		else if( SafeOps::iequals(token[1], "off") ) {
			sys.surf_descent = 0;
			Output::out1("SIM_CONTROL: surf_descent is OFF.\n");
		}
		else 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "ee_local") ) {
		if( SafeOps::iequals(token[1], "on") ) {
			sys.ee_local = 1;
			Output::out1("SIM_CONTROL: Exhaustive enumeration is ON\n");
		}
		else if( SafeOps::iequals(token[1], "off") ) {
			sys.ee_local = 0;
			Output::out1("SIM_CONTROL: Exhaustive enumeration is OFF\n");
		}
		else 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "range_eps") ) {
		if( !SafeOps::atod(token[1], sys.range_eps) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "range_sig") ) {
		if( !SafeOps::atod(token[1], sys.range_sig) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "step_eps") ) {
		if( !SafeOps::atod(token[1], sys.step_eps) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "step_sig") )	{
		if( !SafeOps::atod(token[1], sys.step_sig) )
			return fail;
		return ok;
	}

	//spectre options
	if( SafeOps::iequals(token[0], "spectre") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.spectre = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.spectre = 0;
		else return fail;
		return ok;
	}
	if (SafeOps::iequals(token[0], "spectre_max_charge")) {
		if( !SafeOps::atod(token[1], sys.spectre_max_charge) )
			return fail;
		sys.spectre_max_charge = abs(sys.spectre_max_charge);
		return ok;
	}
	if( SafeOps::iequals(token[0], "spectre_max_target")) {
		if( !SafeOps::atod(token[1], sys.spectre_max_target) )
			return fail;
		sys.spectre_max_target = abs(sys.spectre_max_target);
		return ok;
	}
	
	//cavity options
	if( SafeOps::iequals(token[0], "cavity_bias") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.cavity_bias = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.cavity_bias = 0;
		else return fail; //no match
		return ok;
	}
	if( SafeOps::iequals(token[0], "cavity_grid") ) {
		if ( !SafeOps::atoi(token[1], sys.cavity_grid_size) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "cavity_radius") ) {
		if( !SafeOps::atod(token[1], sys.cavity_radius) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "cavity_autoreject_absolute") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.cavity_autoreject_absolute = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.cavity_autoreject_absolute = 0;
		else return fail; //no match
		return ok;
	}
	if( SafeOps::iequals(token[0], "cavity_autoreject_scale") ) {
		if( !SafeOps::atod(token[1], sys.cavity_autoreject_scale) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "cavity_autoreject_repulsion") ) {
		if( !SafeOps::atod(token[1], sys.cavity_autoreject_repulsion) ) 
			return fail;
		return ok;
	}
	//polar options
	if( SafeOps::iequals(token[0], "polarization") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polarization = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polarization = 0;
		else 
			return fail;
		return ok;
	}
	if(  SafeOps::iequals(token[0], "polarvdw") || SafeOps::iequals(token[0], "cdvdw")  ) {
		if ( SafeOps::iequals(token[1], "on") ) {
			sys.polarvdw = 1;
			sys.polarization = 1;
			sys.polar_iterative = 1; //matrix inversion destroys A_matrix before vdw can use it.
			Output::out1( "SIM_CONTROL: Forcing polar_iterative ON for CP-VdW.\n" );
		}
		else if( SafeOps::iequals(token[1], "evects") ) {
			sys.polarvdw = 2; //calculate eigenvectors
			sys.polarization = 1;
			sys.polar_iterative = 1; //matrix inversion destroys A_matrix before vdw can use it.
			Output::out1( "SIM_CONTROL: Forcing polar_iterative ON for CP-VdW.\n" );
		}
		else if( SafeOps::iequals(token[1], "comp") ) {
			sys.polarvdw = 3; //calculate eigenvectors
			sys.polarization = 1;
			sys.polar_iterative = 1; //matrix inversion destroys A_matrix before vdw can use it.
			Output::out1( "SIM_CONTROL: Forcing polar_iterative ON for CP-VdW.\n" );
		}
		else if( SafeOps::iequals(token[1], "off"))
			sys.polarvdw = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "cdvdw_9th_repulsion") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.cdvdw_9th_repulsion = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.cdvdw_9th_repulsion = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "cdvdw_exp_repulsion") ) {
		if( SafeOps::iequals(token[1], "on"))
			sys.cdvdw_exp_repulsion = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.cdvdw_exp_repulsion = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "cdvdw_sig_repulsion") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.cdvdw_sig_repulsion = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.cdvdw_sig_repulsion = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_ewald_full") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_ewald_full = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_ewald_full = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_ewald") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_ewald = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_ewald = 0;
		else return fail;
		return ok;
	}
	//polar wolf shiz
	if( SafeOps::iequals(token[0], "polar_wolf_full") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_wolf_full = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_wolf_full = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_wolf") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_wolf = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_wolf = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_wolf_alpha_lookup") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_wolf_alpha_lookup = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_wolf_alpha_lookup = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_wolf_damp") ) {
		if( !SafeOps::atod(token[1], sys.polar_wolf_alpha) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_wolf_alpha") ) {  // Same as polar_wolf_damp
		if( !SafeOps::atod(token[1], sys.polar_wolf_alpha) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_wolf_alpha_lookup_cutoff") ) {  // Store data in lookup table
		if( !SafeOps::atod(token[1], sys.polar_wolf_alpha_lookup_cutoff) )
			return fail;
		return ok;
	}
	// replay options
	if( SafeOps::iequals(token[0], "calc_pressure") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.calc_pressure = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.calc_pressure = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "calc_pressure_dv") ) {
		if( !SafeOps::atod(token[1], sys.calc_pressure_dv) ) 
			return fail;
		return ok;
	}
	// set total energy for NVE
	if( SafeOps::iequals(token[0], "total_energy") ) {
		if( !SafeOps::atod(token[1], sys.total_energy) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "numsteps") ) {
		if( !SafeOps::atoi(token[1], sys.numsteps) ) 
			return fail;
		if( sys.numsteps < 1 )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "corrtime") ) {
		if( !SafeOps::atoi(token[1], sys.corrtime) ) 
			return fail;
		return ok;
	}
	// set Monte Carlo options 
	if( SafeOps::iequals(token[0], "move_factor") ) {
		if( !SafeOps::atod(token[1], sys.move_factor) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "move_probability"))	{
		Output::err("move_probability is no longer supported as this is not a probability, but a maximum factor by which to scale the length of random moves. Use move_factor instead.\n");
		return fail;
	}
	if( SafeOps::iequals(token[0], "rot_probability") ) {
		Output::out1("rot_probability is no longer supported as this is not a probability, but the maximum rotation that can occur as a Monte Carlo rotational move. Use rot_factor instead.\n"); 
		return fail;
	}
	if( SafeOps::iequals(token[0], "rot_factor") ) {
		if( !SafeOps::atod(token[1], sys.rot_factor) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "gwp_probability") ) {
		if( !SafeOps::atod(token[1], sys.gwp_probability) )
			return fail;
		return ok;
	}
	if (SafeOps::iequals( token[0], "insert_probability") ) {
		if ( !SafeOps::atod(token[1], sys.insert_probability) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "adiabatic_probability") ) {
		if( !SafeOps::atod(token[1], sys.adiabatic_probability) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "spinflip_probability") ) {
		if( !SafeOps::atod(token[1], sys.spinflip_probability) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals( token[0], "volume_probability") ) {
		if( !SafeOps::atod(token[1], sys.volume_probability) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "volume_change_factor") ) {
		if( !SafeOps::atod(token[1], sys.volume_change_factor) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals( token[0], "transfer_probability") ) {
		if( !SafeOps::atod(token[1], sys.transfer_probability) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals( token[0], "bead_perturb_probability" )) {
		if( !SafeOps::atod(token[1], sys.bead_perturb_probability ))
			return fail;
		return ok;
	}
	if (SafeOps::iequals(token[0], "bead_perturb_factor")) {
		if (!SafeOps::atod(token[1], sys.PI_bead_perturb_factor ))
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "PI_trial_chain_length" )) {
		if( !SafeOps::atoi( token[1], PI_trial_chain_length ))
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "ptemp_freq") ) {
		if( !SafeOps::atoi(token[1], sys.ptemp_freq) )
			return fail;
		return ok;
	}
#ifndef QM_ROTATION
	sys.spinflip_probability = 0;
#endif // !QM_ROTATION

	
	// parallel tempering options 
	if( SafeOps::iequals(token[0], "parallel_tempering") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.parallel_tempering = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.parallel_tempering = 0;
		else return fail; //no match
		return ok;
	}
	if( SafeOps::iequals(token[0], "max_temperature") ) {
		if( !SafeOps::atod(token[1], sys.max_temperature) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "temperature") ) {
		if( !SafeOps::atod(token[1], sys.temperature)) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "simulated_annealing") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.simulated_annealing = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.simulated_annealing = 0;
		else return fail; //no match
		return ok;
	}

	if( SafeOps::iequals(token[0], "simulated_annealing_linear") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.simulated_annealing_linear = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.simulated_annealing_linear = 0;
		else return fail; //no match
		return ok;
	}

	if( SafeOps::iequals(token[0], "simulated_annealing_schedule") ) {
		if( !SafeOps::atod( token[1], sys.simulated_annealing_schedule) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "simulated_annealing_target") ) {
		if( !SafeOps::atod(token[1], sys.simulated_annealing_target) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "pressure") ) {
		if( !SafeOps::atod(token[1], sys.pressure )) 
			return fail;
		return ok;
	}
	
	// fugacity settings
	if( SafeOps::iequals(token[0], "h2_fugacity") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.h2_fugacity = 1;
		else if(SafeOps::iequals(token[1], "off") )
			sys.h2_fugacity = 0;
		else return fail;
		return ok;
	}
	if(SafeOps::iequals(token[0], "co2_fugacity") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.co2_fugacity = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.co2_fugacity = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "ch4_fugacity") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.ch4_fugacity = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.ch4_fugacity = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "n2_fugacity") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.n2_fugacity = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.n2_fugacity = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "user_fugacities") ) { //read in a list of user-defined fugacities
		sys.user_fugacities = 1;
		sys.fugacitiesCount = 0;
		for(int i=1; i<maxTokens; i++)
			if( strlen(token[i]) )
				sys.fugacitiesCount++;
		if( sys.fugacitiesCount == 0 ) 
			return fail;

		for( int i=1; i<=sys.fugacitiesCount; i++)
			if( !SafeOps::atod(token[i], sys.fugacities[i-1]) )
				return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "free_volume") ) {
		if( !SafeOps::atod(token[1], sys.free_volume) ) 
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "rd_only") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.rd_only = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.rd_only = 0;
		else return fail;
		return ok;
	}

	if( SafeOps::iequals(token[0], "gwp") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.gwp = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.gwp = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "wolf") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.wolf = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.wolf = 0;
		else return fail;
		return ok;
	}
	// rd options
	if( SafeOps::iequals(token[0], "rd_lrc") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.rd_lrc = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.rd_lrc = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "rd_crystal") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.rd_crystal = 1;
		else if( SafeOps::iequals(token[1], "off"))
			sys.rd_crystal = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "rd_crystal_order") ) {
		if( !SafeOps::atoi(token[1], sys.rd_crystal_order) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "rd_anharmonic") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.rd_anharmonic = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.rd_anharmonic = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "rd_anharmonic_k") ) {
		if( !SafeOps::atod(token[1], sys.rd_anharmonic_k) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "rd_anharmonic_g") ) {
		if( !SafeOps::atod(token[1], sys.rd_anharmonic_g) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "feynman_hibbs") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.feynman_hibbs = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.feynman_hibbs = 0;
		else return fail;
		return ok;
	}
	// auxillary feynman-hibbs correction for polarvdw (the default correction is better)
	if( SafeOps::iequals(token[0], "vdw_fh_2be") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.vdw_fh_2be = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.vdw_fh_2be = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "feynman_kleinert") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.feynman_kleinert = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.feynman_kleinert = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "feynman_hibbs_order") ) {
		if( !SafeOps::atoi(token[1], sys.feynman_hibbs_order) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "sg") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.use_sg = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.use_sg = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "waldmanhagler") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.waldmanhagler = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.waldmanhagler = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "halgren_mixing") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.halgren_mixing = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.halgren_mixing = 0;
		else return fail;
		return ok;
	}

	if( SafeOps::iequals(token[0], "dreiding") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.use_dreiding = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.use_dreiding = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "lj_buffered_14_7") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.using_lj_buffered_14_7 = true;
		else if( SafeOps::iequals(token[1], "off") )
			sys.using_lj_buffered_14_7 = false;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "disp_expansion") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.using_disp_expansion = true;
		else if( SafeOps::iequals(token[1], "off") )
			sys.using_disp_expansion = false;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "extrapolate_disp_coeffs") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.extrapolate_disp_coeffs = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.extrapolate_disp_coeffs = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "damp_dispersion") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.damp_dispersion = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.damp_dispersion = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "disp_expansion_mbvdw") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.disp_expansion_mbvdw = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.disp_expansion_mbvdw = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "axilrod_teller") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.using_axilrod_teller = true;
		else if( SafeOps::iequals(token[1], "off"))
			sys.using_axilrod_teller = false;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "midzuno_kihara_approx") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.midzuno_kihara_approx = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.midzuno_kihara_approx = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "schmidt_mixing") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.schmidt_mixing = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.schmidt_mixing = 0;
		else return fail;
		return ok;
	}
	if (SafeOps::iequals(token[0], "force_mixing")) {
		if (SafeOps::iequals(token[1], "on"))
			sys.force_mixing = 1;
		else if (SafeOps::iequals(token[1], "off"))
			sys.force_mixing = 0;
		else return fail;
		return ok;
	}
	if (SafeOps::iequals(token[0], "bohm_ahlrichs_mixing")) {
		if (SafeOps::iequals(token[1], "on"))
			sys.bohm_ahlrichs_mixing = 1;
		else if (SafeOps::iequals(token[1], "off"))
			sys.bohm_ahlrichs_mixing = 0;
		else return fail;
		return ok;
	}
	if (SafeOps::iequals(token[0], "wilson_popelier_mixing")) {
		if (SafeOps::iequals(token[1], "on"))
			sys.wilson_popelier_mixing = 1;
		else if (SafeOps::iequals(token[1], "off"))
			sys.wilson_popelier_mixing = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "c6_mixing") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.c6_mixing = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.c6_mixing = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "wrapall") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.wrapall = 1;
		else if( SafeOps::iequals( token[1], "off") )
			sys.wrapall = 0;
		else return fail; // no match
		return ok;
	}
	if( SafeOps::iequals(token[0], "scale_charge") ) {
		if( !SafeOps::atod(token[1], sys.scale_charge) )
			return fail;
		return ok;
	}

	if( SafeOps::iequals(token[0], "ewald_alpha") ) {
		if( !SafeOps::atod(token[1], sys.ewald_alpha) )
			return fail;
		sys.ewald_alpha_set = 1;
		return ok;
	}

	if( SafeOps::iequals(token[0], "ewald_kmax") ) {
		if( !SafeOps::atoi(token[1], sys.ewald_kmax) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "pbc_cutoff") ) {
		if(  ! SafeOps::atod(token[1], sys.pbc.cutoff )  )
			return fail;
		return ok;
	}
	//polar options
	if( SafeOps::iequals(token[0], "polar_ewald") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_ewald = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_ewald = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_ewald_alpha") ) {
		if( !SafeOps::atod(token[1], sys.polar_ewald_alpha) )
			return fail;
		sys.polar_ewald_alpha_set = 1;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polarizability_tensor") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polarizability_tensor = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polarizability_tensor = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_zodid") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_zodid = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_zodid = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_iterative") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_iterative = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_iterative = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_palmo") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_palmo = 1;
		else if (SafeOps::iequals(token[1], "off") )
			sys.polar_palmo = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_gs") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_gs = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_gs = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_gs_ranked") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_gs_ranked = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_gs_ranked = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_sor") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_sor = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_sor = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_esor") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_esor = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_esor = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_gamma") ) {
		if( !SafeOps::atod(token[1], sys.polar_gamma) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_damp") ) {
		if( !SafeOps::atod(token[1], sys.polar_damp) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_precision") ) {
		if( !SafeOps::atod(token[1], sys.polar_precision) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_max_iter") ) {
		if( !SafeOps::atoi(token[1], sys.polar_max_iter) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_damp_type") ) {
		if( SafeOps::iequals(token[1], "none") )
			sys.damp_type = DAMPING_OFF;
		else if( SafeOps::iequals(token[1], "off") )
			sys.damp_type = DAMPING_OFF;
		else if( SafeOps::iequals(token[1], "linear") )
			sys.damp_type = DAMPING_LINEAR;
		else if( SafeOps::iequals(token[1], "exponential") )
			sys.damp_type = DAMPING_EXPONENTIAL;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "polar_rrms") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.polar_rrms = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.polar_rrms = 0;
		else return fail;
		return ok;
	}

	if( SafeOps::iequals(token[0], "cuda") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.cuda = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.cuda = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "independent_particle") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.independent_particle = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.independent_particle = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals( token[0], "pqr_input") ) {
		if( strlen(token[1]) )
			strcpy(sys.pqr_input, token[1]);
		else return fail;
		return ok;
	}
	if (SafeOps::iequals(token[0], "pqr_input_B")) {
		if (strlen(token[1]))
			strcpy(sys.pqr_input_B, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "pqr_output") ) {
		if( strlen(token[1]) )
			strcpy(sys.pqr_output, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "pqr_restart") ) {
		if( strlen(token[1]) )
			strcpy(sys.pqr_restart, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "traj_output") ) {
		if( strlen(token[1]) )
			strcpy(sys.traj_output, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "traj_input")) {
		if( strlen(token[1]))
			strcpy(sys.traj_input, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "energy_output") ) {
		if( strlen(token[1]) )
			strcpy(sys.energy_output, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "energy_output_csv") ) {
		if( strlen(token[1]) )
			strcpy(sys.energy_output_csv, token[1]);
		else return fail;
		return ok;
	}
	if (SafeOps::iequals(token[0], "xyz_output")) {
		if (strlen(token[1]))
			strcpy(sys.xyz_output, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "pop_histogram_output") ) {
		if( strlen(token[1]) )
			strcpy(sys.histogram_output, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "dipole_output") ) {
		if( strlen(token[1]) )
			strcpy(sys.dipole_output, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "field_output") ) {
		if( strlen(token[1]) )
			strcpy(sys.field_output, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "frozen_output") ) {
		if( strlen(token[1]) )
			strcpy(sys.frozen_output, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "insert_input") ) {
		if( strlen(token[1]) )
			strcpy(sys.insert_input, token[1]);
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "surf_output") ) {
		if( strlen(token[1]) ) 
			strcpy(sys.surf_output, token[1]);
		else return fail;
		return ok;
	}
	// Force long (non-PDB compliant, %11.6f) output of coordinates
	if( SafeOps::iequals(token[0], "long_output") ) {
		if( SafeOps::iequals(token[1], "on") ) {
			sys.long_output = 1;
			Output::out1("SIM_CONTROL: Long coordinate output requested.\n");
		}
		else if( SafeOps::iequals(token[1], "off") ) 
			sys.long_output = 0;
		else return fail;
		return ok;
	}
	// read box limits from pqr input
	if( SafeOps::iequals(token[0], "read_pqr_box") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.read_pqr_box_on = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.read_pqr_box_on = 0;
		else return fail;
		return ok;
	}
	// surface fit input parameters
	if( SafeOps::iequals(token[0], "fit_schedule") ) {
		if( !SafeOps::atod(token[1], sys.fit_schedule) )
			return fail;
		if( sys.fit_schedule <= 0.0 || sys.fit_schedule >= 1.0 ) {
			Output::err("SIM_CONTROL: Invalid fit_schedule. Acceptable values are (0,1]\n");
			return fail;
		}
		return ok;
	}
	if( SafeOps::iequals(token[0], "fit_max_energy") ) {
		if( !SafeOps::atod(token[1], sys.fit_max_energy) )
			return fail;
		if( sys.fit_max_energy <= 0.0 ) {
			Output::err("SIM_CONTROL: fit_max_energy parameter must be greater than zero.\n");
			return fail;
		}
		return ok;
	}
	if( SafeOps::iequals(token[0], "fit_start_temp")) {
		if( !SafeOps::atod(token[1], sys.fit_start_temp) )
			return fail;
		if( sys.fit_start_temp <= 0.0 ) {
			Output::err("SIM_CONTROL: fit_start_temp parameter must be greater than zero.\n");
			return fail;
		}
		return ok;
	}
	if( SafeOps::iequals(token[0], "fit_boltzmann_weight") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.fit_boltzmann_weight = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.fit_boltzmann_weight = 0;
		else return fail;
		return ok;
	}
	
	if( SafeOps::iequals(token[0], "max_bondlength") ) {
		if( !SafeOps::atod(token[1], sys.max_bondlength) )
			return fail;
		return ok;
	}
	
	// set basis
	// normal way
	if( SafeOps::iequals(token[0], "basis1") ) {
		if( !SafeOps::atod(token[1], sys.pbc.basis[0][0]) ) return fail;
		if( !SafeOps::atod(token[2], sys.pbc.basis[0][1]) ) return fail; 
		if( !SafeOps::atod(token[3], sys.pbc.basis[0][2]) ) return fail; 
		return ok;
	}
	if( SafeOps::iequals(token[0], "basis2")) {
		if( !SafeOps::atod(token[1], sys.pbc.basis[1][0]) ) return fail; 
		if( !SafeOps::atod(token[2], sys.pbc.basis[1][1]) ) return fail; 
		if( !SafeOps::atod(token[3], sys.pbc.basis[1][2]) ) return fail; 
		return ok;
	}
	if( SafeOps::iequals(token[0], "basis3") ) {
		if( !SafeOps::atod(token[1], sys.pbc.basis[2][0]) ) return fail;
		if( !SafeOps::atod(token[2], sys.pbc.basis[2][1]) ) return fail;
		if( !SafeOps::atod(token[3], sys.pbc.basis[2][2]) ) return fail;
		return ok;
	}
	// .car file way
	// so if both carbasis and basis1/2/3 are in the input file, the last one will overwrite
	if (SafeOps::iequals(token[0], "carbasis")) {
		sys.car2basis(atof(token[1]), atof(token[2]), atof(token[3]), atof(token[4]), atof(token[5]), atof(token[6]));
	}
	if( SafeOps::iequals(token[0], "pop_histogram") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.calc_hist = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.calc_hist = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "pop_hist_resolution") ) {
		if( !SafeOps::atod(token[1], sys.hist_resolution) ) 
			return fail;
		return ok;
	}

	
#ifdef QM_ROTATION
	//quantum rotation stuff
	if( SafeOps::iequals(token[0], "quantum_rotation") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.quantum_rotation = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.quantum_rotation = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "quantum_rotation_hindered") ) {
		if( SafeOps::iequals(token[1], "on") )
			sys.quantum_rotation_hindered = 1;
		else if( SafeOps::iequals(token[1], "off") )
			sys.quantum_rotation_hindered = 0;
		else return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "quantum_rotation_hindered_barrier") ) {
		if( !SafeOps::atod(token[1], sys.quantum_rotation_hindered_barrier) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "quantum_rotation_B") ) {
		if( !SafeOps::atod(token[1], sys.quantum_rotation_B) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "quantum_rotation_level_max") ) {
		if (!SafeOps::atoi(token[1], sys.quantum_rotation_level_max))
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "quantum_rotation_l_max") ) {
		if( SafeOps::atoi(token[1], sys.quantum_rotation_l_max) )
			return fail;
		return ok;
	}
	if( SafeOps::iequals(token[0], "quantum_rotation_sum") ) {
		if( SafeOps::atoi(token[1], sys.quantum_rotation_sum) )
			return fail;
		return ok;
	}
#endif //end QM rotation

	/* #ifdef XXX
	if( SafeOps::iequals(token[0], "quantum_vibration") ) {
		if( SafeOps::iequals(token[1],"on") )
			sys.quantum_vibration = 1;
		else if( SafeOps::iequals(token[1],"off") )
			sys.quantum_vibration = 0;
		else return fail;
		return ok;
	}
	#endif */
	
	/********************************************************************************************
	 *  THIS IS A PROBLEM FOR FUTURE-BRANT---2016.03.09                                         *
	 ********************************************************************************************

	if( SafeOps::iequals(token[0], "fit_input") ) {
		// navigate to the end of the input file linked list
		fileNode_t *node = &(system->fit_input_list);
		while (node->next)
			node = node->next;
		// allocate a new node
		if (!(node->next = malloc(sizeof(fileNode_t)))) {
			error("SIM_CONTROL: Exhausted memory during input file node allocation.\n");
			return (1);
		}
		// advance to new node and initialize
		node = node->next;
		node->next = 0; // terminate list
		if (!(node->data.filename = calloc(MAXLINE, sizeof(char)))) {
			error("SIM_CONTROL: Exhausted memory during string allocation for fit input filename.\n");
			return (1);
		}
		// copy filename to node and increment list count
		strcpy(node->data.filename, token[1]);
		system->fit_input_list.data.count++;
	}
	*/

	
	Output::err("SIM_CONTROL: Unknown Command.\n");
	return fail;
}




bool SimulationControl::check_system() {
// Validates data that was read in from the system file. Checks for invalid or 
// contradictory options, etc. 

	char linebuf[maxLine];

	switch( sys.ensemble ) {
		case ENSEMBLE_UVT:
			Output::out1("SIM_CONTROL: Grand canonical ensemble\n");
			if( ! check_mc_options() )
				return fail;
			break;
		case ENSEMBLE_NVT:
			Output::out1("SIM_CONTROL: Canonical ensemble\n");
			if( ! check_mc_options() )
				return fail;
			break;
		case ENSEMBLE_PATH_INTEGRAL_NVT:
			Output::out1( "SIM_CONTROL: Canonical ensemble for Path Integrals\n" );
			if(  ! check_mc_options()   ||   ! check_PI_options()  )
				return fail;
			break;
		case ENSEMBLE_NVT_GIBBS:
			Output::out1( "SIM_CONTROL: Gibbs ensemble\n" );
			if(  ! check_mc_options() || ! check_Gibbs_options() )
				return fail;
			break;
		case ENSEMBLE_SURF:
			Output::out1("SIM_CONTROL: Potential energy surface\n");
			//ensemble_surf_options(system);/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			break;
		case ENSEMBLE_SURF_FIT:
			Output::out1("SIM_CONTROL: Potential energy surface fitting\n");
			// ensemble_surf_fit_options(system);/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			break;
		case ENSEMBLE_NVE:
			Output::out1("SIM_CONTROL: Microcanonical ensemble\n");
			if( ! check_mc_options() )
				return fail;
			break;
		case ENSEMBLE_TE:
			Output::out1("SIM_CONTROL: Single-point energy calculation\n");
			sys.numsteps = 0;
			sys.corrtime = 0;
			// ensemble_te_options(system);/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			break;
		case ENSEMBLE_NPT:
			Output::out1("SIM_CONTROL: Isobaric-Isothermal ensemble\n");
			if( ! check_mc_options() )
				return fail;
			break;
		case ENSEMBLE_REPLAY:
			Output::out1("SIM_CONTROL: Replaying trajectory\n");
			// ensemble_replay_options(system); /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			break;
		default:
			Output::err("SIM_CONTROL: improper ensemble specified.\n");
			return fail;		
	}

	if( sys.spectre        ) 
		if( ! check_spectre_options() )
			return fail;
	if( sys.rd_only )
		Output::out1("SIM_CONTROL: calculating repulsion/dispersion only\n");
	if( sys.wolf )
		Output::out1("SIM_CONTROL: ES Wolf summation active\n");
	if( sys.rd_lrc ) 
		Output::out1("SIM_CONTROL: rd long-range corrections are ON\n");
	else 
		Output::out1("SIM_CONTROL: rd long-range corrections are OFF\n");
	if( sys.rd_crystal ) {
		if( sys.rd_crystal_order <= 0 ) {
			Output::err("SIM_CONTROL: rd crystal order must be positive\n");
			return fail;
		}
		else {
			sprintf(linebuf,"SIM_CONTROL: rd crystal order set to %d.\n", sys.rd_crystal_order);
			Output::out1(linebuf);
		}
	}
	if( sys.use_sg )
		Output::out1("SIM_CONTROL: Molecular potential is Silvera-Goldman\n");
	if( sys.waldmanhagler )
		Output::out1("SIM_CONTROL: Using Waldman-Hagler mixing rules for LJ-interactions.\n");
	if( sys.halgren_mixing )
		Output::out1("SIM_CONTROL: Using Halgren mixing rules for LJ-interactions.\n");
	if( sys.c6_mixing )
		Output::out1("SIM_CONTROL: Using C6 mixing rules for LJ-interactions.\n");
	if(
		   ( sys.waldmanhagler   &&  sys.halgren_mixing )
		|| ( sys.waldmanhagler   &&  sys.c6_mixing      )
		|| ( sys.halgren_mixing  &&  sys.c6_mixing      ) 
	) { 
		Output::err("SIM_CONTROL: more than one mixing rule specified\n");
		return fail;
	}
	if( sys.use_dreiding )
		Output::out1("SIM_CONTROL: Molecular potential is DREIDING\n");
	if( sys.using_lj_buffered_14_7 )
		Output::out1("SIM_CONTROL: Molecular potential is lj_buffered_14_7\n");
	if( sys.using_lj_buffered_14_7 )
		Output::out1("SIM_CONTROL: Using Halgren mixing rules for LJ-interactions.\n");
	if( sys.using_disp_expansion )
		Output::out1("SIM_CONTROL: Using the dispersion coefficient expansion and exponential repulsion for LJ-interactions.\n");
	if( sys.extrapolate_disp_coeffs ) 
		Output::out1("SIM_CONTROL: Extrapolating the C10 coefficient from the C6 and C8 coefficients with disp_expansion.\n");
	if( sys.damp_dispersion )
		Output::out1("SIM_CONTROL: Using Tang-Toennies damping for dispersion interactions with disp_expansion.\n");
	if( sys.schmidt_mixing )
		Output::out1("SIM_CONTROL: Using the Schmidt mixing rule for exponential repulsions with disp_expansion.\n");
	if (sys.force_mixing)
		Output::out1("SIM_CONTROL: Using the force matching mixing rule for exponential repulsions with disp_expansion.\n");
	if (sys.bohm_ahlrichs_mixing)
		Output::out1("SIM_CONTROL: Using the Bohm Ahlrichs mixing rule for exponential repulsions with disp_expansion.\n");
	if (sys.wilson_popelier_mixing)
		Output::out1("SIM_CONTROL: Using the Wilson Popelier mixing rule for exponential repulsions with disp_expansion.\n");
	if( sys.feynman_hibbs   &&   ! check_feynman_hibbs_options() ) 
		return fail;
	if( sys.simulated_annealing   &&   ! check_simulated_annealing_options() ) 
		return fail;
	if( sys.calc_hist   &&   ! check_hist_options() ) 
		return fail;
	if( sys.polarization   &&   ! check_polarization_options() ) 
		return fail;
	#ifdef QM_ROTATION
		if( sys.quantum_rotation && ! check_qrot_options() )
			return fail;
	#endif 
	#ifdef XXX
		if( quantum_vibration ) 
			Output::out1( "SIM_CONTROL: Quantum vibrational eigenspectrum calculation enabled\n" );
	#endif // XXX

	// Require a job name (CRC)
	if( sys.job_name[0] == 0) {
		Output::err("SIM_CONTROL: must specify a job name\n");
		return fail;
	} else {
		sprintf(linebuf, "SIM_CONTROL: Job Name: %s\n", sys.job_name);
		Output::out1(linebuf);
		check_io_files_options();
	}

	//miscellaneous options

	if( sys.gwp ) {
		Output::out1("SIM_CONTROL: Gaussian wavepacket code active\n");
		if( sys.gwp_probability == 0.) {
			Output::out1("SIM_CONTROL: GWP move scaling not input - setting equal to move_factor\n");
			sys.gwp_probability = sys.move_factor;
		}
	}

	if( sys.scale_charge != 1.0 ) {
		sprintf(linebuf, "SIM_CONTROL: frozen atom charges scaled by %.2f\n", sys.scale_charge );
		Output::out1(linebuf);
	}
	
	if( sys.cuda ) {
		#ifndef CUDA
			Output::err("SIM_CONTROL: cuda keyword enabled, but not compiled into this code\n");
			return fail;
		#else
			Output::out1("SIM_CONTROL: CUDA GPU acceleration activated\n");
		#endif // CUDA 
	}

	if( sys.rd_anharmonic ) {
		if( !sys.rd_only ) {
			Output::err("SIM_CONTROL: rd_anharmonic being set requires rd_only\n");
			return fail;
		} else {
			sprintf(linebuf, "SIM_CONTROL: rd_anharmonic_k = %.3f K/A^2\n", sys.rd_anharmonic_k);
			Output::out1(linebuf);
			sprintf(linebuf, "SIM_CONTROL: rd_anharmonic_g = %.3f K/A^4\n", sys.rd_anharmonic_g);
			Output::out1(linebuf);
		}
	}

	return ok;
}




bool SimulationControl::check_mc_options( ) {
// Validates options specifically dealing with Monte Carlo 


	char linebuf[maxLine];
	

	// Check for valid simulation steps, correlation time, free volume and temperature
	// (temp for all ensembles except NVE)                  
	
	if( sys.numsteps>=1 ) {
		sprintf(linebuf, "SIM_CONTROL: Each core performing %d simulation steps.\n", sys.numsteps);
		Output::out1(linebuf);
	}
	else {
		Output::err("SIM_CONTROL: Improper number of steps specified.\n");
		return fail;
	}


	if( sys.corrtime >= 1 )  {
		sprintf(linebuf, "SIM_CONTROL: System correlation time is %d steps.\n", sys.corrtime);
		Output::out1(linebuf); 
	} 
	else {
		Output::err("SIM_CONTROL: Improper correlation time specified.\n");
		return fail;
	}


	if( sys.free_volume > 0.0 ) {
		sprintf(linebuf, "SIM_CONTROL: system free_volume is %.3f A^3\n", sys.free_volume);
		Output::out1(linebuf);
	}
	//else?? shouldn't sys.free_volume always be > 0.0?


	if( sys.ensemble != ENSEMBLE_NVE) {

		if( sys.temperature > 0.0 ) {
			sprintf(linebuf, "SIM_CONTROL: system temperature is %.3f K\n", sys.temperature);
			Output::out1(linebuf); 
		}
		else {
			Output::err("SIM_CONTROL: Invalid temperature specified.\n");
			return fail;
		}
	}
	
	
	
	

	 //\   PARALLEL TEMPERING OPTIONS HAVE NOT BEEN TESTED for PI CODE
	//  \_______________________________________________________________________________________________________________________________________

	/*
	if( sys.parallel_tempering ) {
		if( size > 1 ) {
			if( !sys.ptemp_freq )
				sys.ptemp_freq = ptemp_freq_default;

			SafeOps::calloc( sys.ptemp, 1, __LINE__, __FILE__ );
			SafeOps::calloc( sys.ptemp->templist, size, __LINE__, __FILE__ );
		
			if ( sys.max_temperature > sys.temperature ) {
				Output::out1("SIM_CONTROL: Parallel tempering activated\n");
				sprintf(linebuf, "SIM_CONTROL: Parallel tempering frequency set to %d steps.\n", sys.ptemp_freq);
				Output::out1(linebuf);
			}
			else {
				Output::err("SIM_CONTROL: Parallel tempering requires max_temperature > temperature.\n");
				return fail;
			}
			if( sys.ensemble == ENSEMBLE_NVE ) {
				Output::err("SIM_CONTROL: Parallel tempering is not implemented for NVE.\n");
				return fail;
			}
			if( sys.simulated_annealing ) {
				Output::err("SIM_CONTROL: Parallel tempering is incompatible with simulated annealing.\n");
				return fail;
			}
			if( size < 2 ) {
				Output::err("SIM_CONTROL: Multiple cores are required for parallel tempering.\n");
				return fail;
			}
			if( sys.feynman_hibbs ) {
				Output::err("SIM_CONTROL: Parallel tempering is not compatible with temperature-dependent potentials (Feynman Hibbs).\n");
				return fail;
			}
			// set temperature of baths 
			SafeOps::calloc( sys.ptemp->index, size, __LINE__, __FILE__ );
		
			for( int j=0; j<size; j++ ) {
				sys.ptemp->templist[j] = sys.temperature * pow(pow(sys.max_temperature/sys.temperature,1.0/(size-1)), j);
				sys.ptemp->index[j] = j;
			}
			//set local temperature
			sys.temperature = sys.ptemp->templist[rank];

		} 
			
		else { 
			// Size == 1
			Output::err("SIM_CONTROL: Parallel tempering can only be used when running in parallel.\n");
			return fail;
		}
	}
	*/

	


	
	//    ENSEMBLE SPECIFIC CHECKS
	//\________________________________________________________________________________________________________________________

	if( sys.ensemble == ENSEMBLE_NVE) {

		sprintf(linebuf, "SIM_CONTROL: NVE energy is %.3f K\n", sys.total_energy);
		Output::out1(linebuf);
	}



	if(  (sys.ensemble == ENSEMBLE_NVE)  ||  (sys.ensemble == ENSEMBLE_NVT)  ) {
		
		if( sys.spinflip_probability > 1.0) {
			Output::out1("The requested spinflip probabilities is greater than 1.0.\n");
			return fail;
		}

		sprintf(linebuf, "SIM_CONTROL: spinflip probability is %lf.\n", sys.spinflip_probability);
		Output::out1(linebuf);

		sprintf(linebuf, "SIM_CONTROL: displace probability is %lf.\n", 1.0 - sys.spinflip_probability);
		Output::out1(linebuf);
	}



	if( sys.ensemble == ENSEMBLE_PATH_INTEGRAL_NVT ) {
		
		// Check and report relevant probability settings 
		if((sys.spinflip_probability + sys.bead_perturb_probability) > 1.0) {
			Output::err("The requested probabilities for all MC moves sum to a value greater than 1.0.\n");
			return fail;
		}
		sprintf(linebuf, "SIM_CONTROL: spinflip probability is %lf.\n", sys.spinflip_probability);
		Output::out1(linebuf);
		sprintf(linebuf, "SIM_CONTROL: bead perturbation probability is %lf.\n", sys.bead_perturb_probability);
		Output::out1(linebuf);
		sprintf(linebuf, "SIM_CONTROL: displace probability is %lf.\n", 1.0 - sys.spinflip_probability - sys.bead_perturb_probability);
		Output::out1(linebuf);
	}

	

	if( sys.ensemble == ENSEMBLE_NPT ) {

		if( sys.pressure <= 0.0 ) {
			Output::err("SIM_CONTROL: invalid pressure set for NPT\n");
			return fail;

		} else {

			sprintf(linebuf, "SIM_CONTROL: reservoir pressure is %.3f atm\n", sys.pressure);
			Output::out1(linebuf);
		}

		if (sys.volume_probability == 0.0)
		{
			sprintf(linebuf, "SIM_CONTROL: volume change probability is 1/N_molecules.\n");
			Output::out1(linebuf);

			sprintf(linebuf, "SIM_CONTROL: displace probability is 1-1/N_molecules.\n");
			Output::out1(linebuf);
		}
		else
		{
			sprintf(linebuf, "SIM_CONTROL: volume change probability is %.3f\n", sys.volume_probability);
			Output::out1(linebuf);

			sprintf(linebuf, "SIM_CONTROL: displace probability is %.3f\n", 1.0 - sys.volume_probability);
			Output::out1(linebuf);
		}

		sprintf(linebuf, "SIM_CONTROL: volume change factor is %lf.\n", sys.volume_change_factor);
		Output::out1(linebuf);
	}



	if( sys.ensemble == ENSEMBLE_UVT ) {

		if( sys.user_fugacities) {

			sprintf(linebuf, "SIM_CONTROL: user defined fugacities are in use.\n");
			Output::out1(linebuf);
			for( int i=0; i<sys.fugacitiesCount; i++ ) {
				sprintf(linebuf, "SIM_CONTROL: fugacity[%d] is set to %.3f atm\n", i, sys.fugacities[i]);
				Output::out1(linebuf);
			}

			if( sys.pressure != 0.0 ) {
				sprintf(linebuf, "SIM_CONTROL: User defined fugacities are not compatible with pressure specification.\n");
				Output::out1(linebuf);
				return fail;
			}
		}

		else if( sys.pressure <= 0.0) {

			Output::err("SIM_CONTROL: invalid pressure set for GCMC\n");
			return fail;

		}
		else {
			
			sprintf(linebuf, "SIM_CONTROL: reservoir pressure is %.3f atm\n", sys.pressure );
			Output::out1(linebuf);
			
			if( sys.h2_fugacity ) {

				if( sys.fugacities[0] != 0.0 ) {
					Output::err( "SIM_CONTROL: h2_fugacity called, but fugacities are already set.\n" );
					return fail;
				}
				
				sys.fugacities[0] = Fugacity::h2_fugacity( sys.temperature, sys.pressure);
				if( sys.h2_fugacity <= 0.0) {
					Output::err("SIM_CONTROL: error in H2 fugacity assignment\n");
					return fail;
				}
				sprintf(linebuf, "SIM_CONTROL: H2 fugacity = %.3f atm\n", sys.fugacities[0]);
				Output::out1(linebuf);
			}

			if( sys.co2_fugacity ) {

				if( sys.fugacities[0] != 0.0 ) {
					Output::err("SIM_CONTROL: co2_fugacity called, but fugacities are already set.\n");
					return fail;
				}
				sys.fugacities[0] = Fugacity::get_peng_robinson_fugacity(sys.temperature, sys.pressure, "co2");
				if( sys.co2_fugacity <= 0.0 ) {
					Output::err("SIM_CONTROL: error in CO2 fugacity assignment\n");
					return fail;
				}

				sprintf(linebuf, "SIM_CONTROL: CO2 fugacity = %.3f atm\n", sys.fugacities[0] );
				Output::out1(linebuf);
			}

			if( sys.ch4_fugacity ) {

				if( sys.fugacities[0] != 0.0 ) {
					Output::err("SIM_CONTROL: ch4_fugacity called, but fugacities are already set.\n");
					return fail;
				}
				sys.fugacities[0] = Fugacity::ch4_fugacity(sys.temperature, sys.pressure);
				if( sys.ch4_fugacity <= 0.0 ) {
					Output::err("SIM_CONTROL: error in CH4 fugacity assignment\n");
					return fail;
				}

				sprintf(linebuf, "SIM_CONTROL: CH4 fugacity = %.3f atm\n", sys.fugacities[0]);
				Output::out1(linebuf);
			}

			if( sys.n2_fugacity ) {

				if( sys.fugacities[0] != 0.0 ) {
					Output::err("SIM_CONTROL: n2_fugacity called, but fugacities are already set.\n");
					return fail;
				}
				
				sys.fugacities[0] = Fugacity::n2_fugacity(sys.temperature, sys.pressure);
				if( sys.n2_fugacity <= 0.0 ) {
					Output::err("SIM_CONTROL: error in N2 fugacity assignment\n");
					return fail;
				}

				sprintf( linebuf, "SIM_CONTROL: N2 fugacity = %.3f atm\n", sys.fugacities[0] );
				Output::out1(linebuf);
			}

		} //calculated fugacities

		sprintf(linebuf, "SIM_CONTROL: insert/delete probability is %lf.\n", sys.insert_probability);
		Output::out1(linebuf);

		if (sys.quantum_rotation) {
			sprintf(linebuf, "SIM_CONTROL: spinflip probability is %lf.\n", sys.spinflip_probability * (1.0 - sys.insert_probability));
			Output::out1(linebuf);

			sprintf(linebuf, "SIM_CONTROL: displace probability is %lf.\n", (1.0 - sys.spinflip_probability)*(1.0 - sys.insert_probability));
			Output::out1(linebuf);
		}
		else
		{
			sprintf(linebuf, "SIM_CONTROL: displace probability is %lf.\n", 1.0 - sys.insert_probability);
			Output::out1(linebuf);
		}

	} //ensemble uVT
	


	if (sys.ensemble == ENSEMBLE_NVT_GIBBS)
	{
		// Probabilities are established and reported in setup_Gibbs_Systems()
		// because volume_probability defaults to 1/N if not explicitly set, 
		// and N cannot be computed until setup_Gibbs_Systems() is called.
	}
	




	//   GENERAL CHECKS/PARAMETER REPORTING
	//\________________________________________________________________________________________________________________________


	sprintf(linebuf, "SIM_CONTROL: translation change factor is %.5f\n", sys.move_factor);
	Output::out1(linebuf);

	sprintf(linebuf, "SIM_CONTROL: rotation change factor is %.5f\n", sys.rot_factor);
	Output::out1(linebuf);

	if (ENSEMBLE_PATH_INTEGRAL_NVT == sys.ensemble) {
		sprintf( linebuf,"SIM_CONTROL: bead perturbation scaling factor is %.5f\n", sys.PI_bead_perturb_factor);
		Output::out1(linebuf);
	}

	if( sys.gwp ) {
		sprintf(linebuf, "SIM_CONTROL: gwp change factor is %.3f\n", sys.gwp_probability);
		Output::out1(linebuf);
	}


	if( sys.cavity_autoreject_absolute ) { 
		Output::out1("SIM_CONTROL: cavity autoreject absolute activated\n");
		if ((sys.cavity_autoreject_scale <= 0.0) || (sys.cavity_autoreject_scale > 1.78)) {
			Output::err("SIM_CONTROL: cavity_autoreject_scale either not set or out of range\n");
			return fail;
		}
	}
		

	if( sys.cavity_bias ) {

		if(   (sys.cavity_grid_size <= 0)  ||  (sys.cavity_radius <= 0.0)   ) {
			Output::err("SIM_CONTROL: invalid cavity grid or radius specified\n"); 
			return fail;
		
		} else {
			Output::out1("SIM_CONTROL: cavity-biased umbrella sampling activated\n");
			sprintf(linebuf, "SIM_CONTROL: cavity grid size is %dx%dx%d points with a sphere radius of %.3f A\n", sys.cavity_grid_size, sys.cavity_grid_size, sys.cavity_grid_size, sys.cavity_radius );
			Output::out1(linebuf);
		}
	}

	return ok;
}




bool SimulationControl::check_spectre_options() {
	char linebuf[maxLine];
	if( sys.ensemble != ENSEMBLE_NVT ) {
		Output::err("SIM_CONTROL: SPECTRE algorithm requires canonical ensemble\n");
		return fail;

	} else {

		Output::out1("SIM_CONTROL: SPECTRE algorithm activated\n");
		sprintf(linebuf, "SIM_CONTROL: SPECTRE max charge = %.3f\n", sys.spectre_max_charge);
		Output::out1(linebuf);
		sprintf(linebuf, "SIM_CONTROL: SPECTRE max target = %.3f\n", sys.spectre_max_target);
		Output::out1(linebuf);
		return ok;
	}
}




bool SimulationControl::check_io_files_options() {

	char linebuf[maxLine];
	int file_count = (PI_nBeads) ? PI_nBeads : size;


	if (SafeOps::iequals(sys.pqr_restart, "off")) { // Optionally turn off restart configuration output
		Output::err("SIM_CONTROL: **Warning**: PQR restart file option disabled; writing restart configuration to /dev/null\n");
		strcpy(sys.pqr_restart, "/dev/null");

	} else {

		if (sys.pqr_restart[0] == 0) {	// (CRC)
			// If not filename template specified, generate the default
			strcpy(sys.pqr_restart, sys.job_name);
			strcat(sys.pqr_restart, ".restart.pqr");
		}

		
		if (file_count > 1) {
			for (int j = 0; j < file_count; j++) {
				if (mpi) {
					#ifdef _MPI
						MPI_Barrier(MPI_COMM_WORLD);
						std::string filename = Output::make_filename(sys.pqr_restart, rank);
						pqr_restart_filenames.push_back(filename);
						if (j == rank) {
							sprintf(linebuf, "SIM_CONTROL: Thread/SYSTEM %d will be writing restart configuration to ./%s\n", rank, pqr_restart_filenames[rank].c_str());
							Output::out(linebuf);
						}
					#endif
				}
				else {
					std::string filename = Output::make_filename(sys.pqr_restart, j);
					pqr_restart_filenames.push_back(filename);
					sprintf(linebuf, "SIM_CONTROL: SYSTEM %d will be writing restart configuration to ./%s\n", j, pqr_restart_filenames[j].c_str());
					Output::out1(linebuf);
				}
			}
		}
		else {
			pqr_restart_filenames.push_back(sys.pqr_restart);
			sprintf(linebuf, "SIM_CONTROL: will be writing restart configuration to ./%s\n", pqr_restart_filenames[0].c_str());
			Output::out(linebuf);
		}
	}
	


	if (SafeOps::iequals(sys.pqr_output, "off")) { // Optionally turn off final configuration output
		Output::err("SIM_CONTROL: **Warning: PQR final configuration file disabled; writing to /dev/null\n");
		strcpy(sys.pqr_output, "/dev/null");

	} else {

		if (sys.pqr_output[0] == 0) {
			// If no filename template specified, generate the default
			strcpy(sys.pqr_output, sys.job_name);
			strcat(sys.pqr_output, ".final.pqr");
		}

		if (file_count > 1) {

			for (int j = 0; j < file_count; j++) {
					
				if( mpi ) {
					#ifdef _MPI
					MPI_Barrier(MPI_COMM_WORLD);
					
					std::string filename = Output::make_filename(sys.pqr_output, rank);
					pqr_final_filenames.push_back(filename);
					if (j == rank) {
						sprintf(linebuf, "SIM_CONTROL: Thread/SYSTEM %d will be writing final configuration to ./%s\n", rank, pqr_final_filenames[rank].c_str());
						Output::out(linebuf);
					}
					#endif

				} else {
					std::string filename = Output::make_filename(sys.pqr_output, j);
					pqr_final_filenames.push_back(filename);
					sprintf(linebuf, "SIM_CONTROL: SYSTEM %d will be writing final configuration to ./%s\n", j, pqr_final_filenames[j].c_str());
					Output::out(linebuf);
				}
			}
		} 
		else {
			pqr_final_filenames.push_back(sys.pqr_output);
			sprintf(linebuf, "SIM_CONTROL: will be writing final configuration to ./%s\n", pqr_final_filenames[0].c_str());
			Output::out(linebuf);
		}
	}



	
	if( sys.parallel_restarts ) {

		

		if (file_count > 1) {
			for (int j = 0; j < file_count; j++) {

				FILE *test;

				// Try to open the plain restart file.
				int id = mpi ? rank : j;
				std::string filename = Output::make_filename(sys.pqr_restart, id);
				test = fopen(filename.c_str(), "r");
				if (test) {
					// File was successfully opened so this is the one we will use
					fclose(test);
				
				} else {
					// No restart files available, try to open a "last" file
					std::string basename = Output::make_filename(sys.pqr_restart, rank);
					filename = basename + ".last";
					//SafeOps::calloc(filename, (int)strlen(basename) + 16, sizeof(char), __LINE__, __FILE__);
					//sprintf(filename, "%s.last", basename);
					//free(basename);
					test = fopen(filename.c_str(), "r");
					if (test) {
						fclose(test);

					} else {
						// if pqr_input is non-NULL, the user specified a default input file 
						// and we will use that, exiting the if/else tree here
						if (sys.pqr_input[0] == 0) {
							// if they did not specify a default input geometry file (and evidently none of the parallel options
							// were available) then we will construct the default geometry input file name, based on the job name.
							filename = std::string(sys.job_name) + ".initial.pqr";
						} 
					} // default
				} // pqr.last


				if (mpi ) {
					#ifdef _MPI
					MPI_Barrier(MPI_COMM_WORLD);
					pqr_input_filenames.push_back(filename);
					strcpy(sys.pqr_input, filename.c_str());
					if (j == rank) {
						sprintf(linebuf, "SIM_CONTROL: Thread/SYSTEM %d will be reading coordinates from file: /%s\n", rank, sys.pqr_input);
						Output::out(linebuf);
					}
					#endif
				}
				else {
					
					pqr_input_filenames.push_back(filename);
					sprintf(linebuf, "INPUT (system %d): Reading coordinates from file: %s\n", j, pqr_input_filenames[j].c_str());
					Output::out1(linebuf);
				}
			}
		}
	} else {

		if (sys.pqr_input[0] == 0) {
			strcpy(sys.pqr_input, sys.job_name);
			strcat(sys.pqr_input, ".initial.pqr");
			pqr_input_filenames.push_back(sys.pqr_input);
			sprintf(linebuf, "SIM_CONTROL: input PQR file not specified...will try to read coordinates from ./%s\n", pqr_input_filenames[0].c_str());
			Output::out(linebuf);
		}
		else {
			pqr_input_filenames.push_back(sys.pqr_input);
			sprintf(linebuf, "SIM_CONTROL: reading initial molecular coordinates from: %s\n", pqr_input_filenames[0].c_str());
			Output::out(linebuf);
		}
	}





	// Energy output will default to on if not specified

	if( sys.energy_output[0] == 0 ) {  // (CRC)
		strcpy( sys.energy_output, sys.job_name );
		strcat( sys.energy_output, ".energy.dat"    );
	}

	if( SafeOps::iequals( sys.energy_output, "off") ) { // Optionally turn off energy printing
		Output::out( "SIM_CONTROL: energy file output disabled; writing to /dev/null\n" );
		strcpy( sys.energy_output, "/dev/null" );
	} else {
		sprintf(linebuf, "SIM_CONTROL: will be writing energy output to ./%s\n", sys.energy_output);
		Output::out(linebuf);
	}




	if( sys.surf_virial) { //if we are doing virial
		if (sys.virial_output[0] == 0) {
			strcpy(sys.virial_output, sys.job_name);
			strcat(sys.virial_output, ".virial.dat");
		}
		
		if( SafeOps::iequals( sys.virial_output, "off") )  {   //optionally turn off virial printing
			Output::out("SIM_CONTROL: virial file output disabled; writing to /dev/null\n");
			strcpy( sys.virial_output, "/dev/null" );

		} else {
			sprintf(linebuf, "SIM_CONTROL: will be writing virial output to ./%s\n",  sys.virial_output);
			Output::out(linebuf);
		}
	}



	// Trajectory file will default to on if not specified 
	if( sys.traj_output[0] == 0 ) {	// (CRC)
		strcpy( sys.traj_output, sys.job_name );
		strcat( sys.traj_output, ".traj.pqr"      );
	} 

	if( SafeOps::iequals( sys.traj_output, "off")) { // Optionally turn off trajectory printing
		Output::out("SIM_CONTROL: trajectory file output disabled; writing to /dev/null\n");
		strcpy( sys.traj_output, "/dev/null" ); 
	} else {
		sprintf(linebuf, "SIM_CONTROL: will be writing trajectory to ./%s\n", sys.traj_output);
		Output::out(linebuf);
	}




	if( sys.insert_input[0] ) {
		sprintf( linebuf, "SIM_CONTROL: inserted molecules will be selected from ./%s\n", sys.insert_input );
		Output::out( linebuf );
	} 



	if (sys.polarization) {

		if (!sys.dipole_output[0]) {  // (CRC)
			strcpy(sys.dipole_output, sys.job_name);
			strcat(sys.dipole_output, ".dipole.dat");
		}

		if (SafeOps::iequals(sys.dipole_output, "off")) {
			Output::out("SIM_CONTROL: dipole file output disabled; writing to /dev/null\n");
			strcpy(sys.dipole_output, "/dev/null");
		} else {
			sprintf(linebuf, "SIM_CONTROL: dipole field will be written to ./%s\n", sys.dipole_output);
			Output::out(linebuf);
		}



		if (!sys.field_output[0]) { // (CRC)
			strcpy(sys.field_output, sys.job_name);
			strcat(sys.field_output, ".field.dat");
		}

		if (SafeOps::iequals(sys.field_output, "off")) {
			Output::out("SIM_CONTROL: field file output disabled; writing to /dev/null\n");
			strcpy(sys.field_output, "/dev/null");
		} else {
			sprintf(linebuf, "SIM_CONTROL: field field will be written to ./%s\n", sys.field_output);
			Output::out(linebuf);
		}
	}

	return ok;
}




bool SimulationControl::check_feynman_hibbs_options( ) {

	char linebuf[maxLine];
	Output::out("SIM_CONTROL: Feynman-Hibbs effective potential activated\n");

	if( sys.feynman_kleinert ) {
		Output::out("SIM_CONTROL: Feynman-Kleinert iteration method activated\n");

		if( ! sys.rd_anharmonic ) {
			Output::err("SIM_CONTROL: Feynman-Kleinert iteration only implemented for anharmonic oscillator\n");
			return fail;
		}
	} 
	else {
		switch( sys.feynman_hibbs_order ) {
			case 2:
				sprintf(linebuf, "SIM_CONTROL: Feynman-Hibbs second-order quantum correction activated\n");
				Output::out(linebuf);
				break;
			case 4:
				sprintf(linebuf, "SIM_CONTROL: Feynman-Hibbs fourth-order quantum correction activated\n");
				Output::out(linebuf);
				break;
			default:
				Output::out("SIM_CONTROL: Feynman-Hibbs order unspecified or specified with unsupported value--defaulting to h^2\n");
				sys.feynman_hibbs_order = 2;
				break;
		}
	}
	//if using polarvdw and FH, cavity_autoreject_absolute must be on (otherwise energies blow up)
	if(   sys.polarvdw   &&   !(sys.cavity_autoreject_absolute)   &&   ( sys.ensemble != ENSEMBLE_REPLAY)   ) {
		Output::err("SIM_CONTROL: cavity_autoreject_absolute must be used with polarvdw + Feynman Hibbs.\n");
		return fail;
	}

	if( sys.temperature <= 0 ) {
		Output::err("SIM_CONTROL: feynman_hibbs requires positive temperature.\n");
		return fail;
	}

	return ok;
}




bool SimulationControl::check_simulated_annealing_options()
{
	char linebuf[maxLine];
	Output::out("SIM_CONTROL: Simulated annealing active\n");

	if(   (sys.simulated_annealing_schedule < 0.0)  ||  (sys.simulated_annealing_schedule > 1.0)   ) {
		Output::err("SIM_CONTROL: invalid simulated annealing temperature schedule specified\n");
		return fail;
	} else {
		sprintf(linebuf, "SIM_CONTROL: Simulated annealing temperature schedule = %.3f\n", sys.simulated_annealing_schedule);
		Output::out(linebuf);
	}

	if( sys.simulated_annealing_target < 0.0 ) {
		Output::err("SIM_CONTROL: invalid simulated annealing target specified\n");
		return fail;
	} else {
		sprintf(linebuf, "SIM_CONTROL: Simulated annealing target %lfK.", sys.simulated_annealing_target );
		Output::out(linebuf);
	}
	
	if( sys.simulated_annealing_linear ) {
		sprintf(linebuf, "SIM_CONTROL: Simulated annealing using a linear ramp.");
		Output::out(linebuf);
	}

	return ok;
}




bool SimulationControl::check_hist_options()
{
	char linebuf[maxLine];

	static const double defaultMaxBondLength = 1.8;
	static const double defaultHistResolution = 0.7;
	static const double histLowerBound = 0.01;
	static const double histUpperBound = 5.00;
	

	Output::out("SIM_CONTROL: Histogram calculation will be performed.\n");
	if( sys.hist_resolution == 0.0 ){
		Output::out("SIM_CONTROL: No histogram resolution set but histogram calculation requested\n");
		sprintf( linebuf, "SIM_CONTROL: Setting hist_resolution to default value of %0.3fA\n", defaultHistResolution );
		Output::out( linebuf );
		sys.hist_resolution = defaultHistResolution;
	}
	else if(   sys.hist_resolution < histLowerBound   ||   sys.hist_resolution > histUpperBound   ) {
		sprintf( linebuf, "SIM_CONTROL: Requested histogram resolution out of range. Valid range is [%0.3f,%0.3f]\n", histLowerBound, histUpperBound );
		Output::out( linebuf );
		sprintf( linebuf, "SIM_CONTROL: Setting hist_resolution to default value of %0.3fA.\n", defaultHistResolution );
		Output::out( linebuf );
		sys.hist_resolution = defaultHistResolution;
	}
	else if( ! sys.histogram_output[0] ){
		Output::out("SIM_CONTROL: No histogram outputfile selected, defaulting to histogram.dx\n");
		strcpy( sys.histogram_output, "histogram.dx" );
	}
	else{
		sprintf(linebuf,"SIM_CONTROL: histogram resolution set to %.3f A\n", sys.hist_resolution);
		Output::out(linebuf);
	}

	if(  sys.max_bondlength < .5  ){
		Output::out("SIM_CONTROL: max_bondlength either not set or out of range. max_bondlength must be >= 0.5\n");
		sprintf( linebuf, "SIM_CONTROL: setting max_bondlength to default value of %0.3fA\n", defaultMaxBondLength );
		Output::out( linebuf );
		sys.max_bondlength = defaultMaxBondLength;
	}	

	if(  ! sys.frozen_output[0]  ){
		Output::out("SIM_CONTROL: no frozen_output set! setting frozen coordinate output file to frozen.dx\n");
		strcpy( sys.frozen_output, "frozen.dx" );
	} else {
		sprintf(linebuf, "SIM_CONTROL: will be writing frozen coordinates to %s\n", sys.frozen_output);
		Output::out(linebuf);
	}

	return ok;
}




bool SimulationControl::check_polarization_options()
{
	char linebuf[maxLine];

	Output::out("SIM_CONTROL: Thole polarization activated\n");

	if( sys.cuda ) {
		if( ! sys.polar_iterative ) {
			Output::err("SIM_CONTROL: CUDA GPU acceleration available for iterative Thole only, enable polar_iterative\n");
			return fail;
		} 
		else if( sys.damp_type != DAMPING_EXPONENTIAL ) {
			Output::err( "SIM_CONTROL: CUDA GPU accleration available for exponential Thole damping only\n" );
			return fail;
		} else if( ! sys.polar_max_iter ) {
			// XXX disable for 1 iter testing 
			Output::err("SIM_CONTROL: Must set polar_max_iter for CUDA GPU acceleration.\n");
			return fail;
		}
		else
			Output::out("SIM_CONTROL: CUDA GPU Thole SCF solver activated\n");
	}

	if(  sys.polar_iterative   &&   sys.polarizability_tensor  ) {
		Output::err("SIM_CONTROL: iterative polarizability tensor method not implemented\n");
		return fail;
	}

	if(  !(sys.polar_iterative)  &&  sys.polar_zodid  ) {
		Output::err("SIM_CONTROL: ZODID and matrix inversion cannot both be set!\n");
		return fail;
	}

	if(  sys.polar_wolf  ||  sys.polar_wolf_full  ) {
		if ( sys.polar_wolf )
			Output::out("SIM_CONTROL: Polar wolf activated. Thole field calculated using wolf method.\n");
		if( sys.polar_wolf_full ) 
			Output::out("SIM_CONTROL: Full polar wolf treatment activated.\n");
		if( sys.polar_wolf_alpha_lookup ) {
			if( sys.polar_wolf_alpha_lookup_cutoff <= 0 ) {
				sprintf(linebuf,"SIM_CONTROL: invalid polar_wolf_alpha_lookup_cutoff\n");
				return fail;
			}
			sprintf(linebuf,"SIM_CONTROL: Polar wolf alpha will be performed via lookup table with cutoff %lf Ang.\n", sys.polar_wolf_alpha_lookup_cutoff );
			Output::out(linebuf);
		}	
		if(  (sys.polar_wolf_alpha > 0)   &&   (sys.polar_wolf_alpha <= 1)  ) {
			sprintf( linebuf, "SIM_CONTROL: Polar wolf damping set to %lf. (0 is default)\n", sys.polar_wolf_alpha);
			Output::out(linebuf);
		} 
		else if( sys.polar_wolf_alpha != 0 )
		{
			Output::err("SIM_CONTROL: 1 >= polar_wolf_alpha >= 0 is required.\n");
			return fail;
		}
	}

	if( sys.polar_ewald )
		Output::out("SIM_CONTROL: Polar ewald activated. Thole field calculated using ewald method.\n");
	
	if( sys.polar_ewald_full ) 
		Output::out("SIM_CONTROL: Full ewald polarization activated.\n");
	
	if( sys.damp_type == DAMPING_LINEAR )
		Output::out("SIM_CONTROL: Thole linear damping activated\n");
	else if( sys.damp_type == DAMPING_OFF )
		Output::out("SIM_CONTROL: Thole linear damping is OFF\n");
	else if( sys.damp_type == DAMPING_EXPONENTIAL )
		Output::out("SIM_CONTROL: Thole exponential damping activated\n");
	else {
		Output::err("SIM_CONTROL: Thole damping method not specified\n");
		return fail;
	}

	if(   sys.polar_damp <= 0.0    &&    sys.damp_type != DAMPING_OFF   ) {
		Output::err("SIM_CONTROL: damping factor must be specified\n");
		return fail;
	} else if ( sys.polar_damp > 0.0 ) {
		sprintf( linebuf, "SIM_CONTROL: Thole damping parameter is %.4f\n", sys.polar_damp );
		Output::out( linebuf );
	}

	if( sys.polar_iterative ) {

		Output::out( "SIM_CONTROL: Thole iterative solver activated\n" );
		if( sys.polar_zodid ) {
			Output::out( "SIM_CONTROL: ZODID polarization enabled\n" );
		}

		if(  (sys.polar_precision > 0.0)   &&   (sys.polar_max_iter > 0)  ) {
			Output::err( "SIM_CONTROL: cannot specify both polar_precision and polar_max_iter, must pick one\n" );
			return fail;
		}
	
		if( sys.polar_precision < 0.0 ) {
			Output::err( "SIM_CONTROL: invalid polarization iterative precision specified\n" );
			return fail;
		} else if( sys.polar_precision > 0.0 ) {
			sprintf(linebuf, "SIM_CONTROL: Thole iterative precision is %e A*sqrt(KA) (%e D)\n", sys.polar_precision, sys.polar_precision/DEBYE2SKA );
			Output::out( linebuf );
		} else {
			sprintf(linebuf, "SIM_CONTROL: using polar max SCF iterations = %d\n", sys.polar_max_iter );
			Output::out(linebuf);
		}

		if( sys.polar_precision > 0.0  ||  sys.polar_rrms ) 
			Output::out("SIM_CONTROL: polar_rrms activated. Dipole rrms will be reported.\n");
		
		if( sys.polar_sor  &&  sys.polar_esor ) {
			Output::err( "SIM_CONTROL: cannot specify both SOR and ESOR SCF methods\n" );
			return fail;
		}

		if( sys.polar_sor )
			Output::out( "SIM_CONTROL: SOR SCF scheme active\n" );
		else if( sys.polar_esor )
			Output::out( "SIM_CONTROL: ESOR SCF scheme active\n" );

		if( sys.polar_gamma < 0.0 ) {
			Output::err( "SIM_CONTROL: invalid Pre-cond/SOR/ESOR gamma set\n" );
			return fail;
		} else {
			sprintf( linebuf, "SIM_CONTROL: Pre-cond/SOR/ESOR gamma = %.3f\n", sys.polar_gamma );
			Output::out( linebuf );
		}

		if( sys.polar_gs   &&   sys.polar_gs_ranked ) {
			Output::err( "SIM_CONTROL: both polar_gs and polar_gs_ranked cannot be set\n" );
			return fail;
		}

		if( sys.polar_gs )
			Output::out( "SIM_CONTROL: Gauss-Seidel iteration scheme active\n" );
		else if( sys.polar_gs_ranked )
			Output::out( "SIM_CONTROL: Gauss-Seidel Ranked iteration scheme active\n" );

		if( sys.polar_palmo )
			Output::out( "SIM_CONTROL: Polarization energy of Palmo and Krimm enabled\n" );


	} else {

		Output::out( "SIM_CONTROL: Matrix polarization activated\n" );
		if( sys.polarizability_tensor )
			Output::out( "SIM_CONTROL: Polarizability tensor calculation activated\n" );
	}


	if( sys.polarvdw ) {
		Output::out( "SIM_CONTROL: polarvdw (coupled-dipole van der Waals) activated\n" );
		if( sys.feynman_hibbs ) {
			if( sys.vdw_fh_2be )
				Output::out( "SIM_CONTROL: two-body-expansion feynman-hibbs for polarvdw is active\n" );
			else
				Output::out( "SIM_CONTROL: polarvdw feynman-hibbs will be calculated using MPFD\n" );
		}
		if( sys.cdvdw_exp_repulsion )
			Output::out( "SIM_CONTROL: exponential repulsion activated\n" );
		if( sys.cdvdw_sig_repulsion )
			Output::out( "SIM_CONTROL: C_6*sig^6 repulsion activated\n" );
		if( sys.cdvdw_9th_repulsion )
			Output::out( "SIM_CONTROL: 9th power repulsion mixing activated\n" );
		if(  (sys.cdvdw_exp_repulsion + sys.cdvdw_sig_repulsion + sys.cdvdw_9th_repulsion + sys.waldmanhagler + sys.halgren_mixing)   > 1   ) {
			Output::err( "SIM_CONTROL: more than one mixing rules specified" );
			return fail;
		}
	}
	if ( ! sys.polarvdw ) {
		if( sys.cdvdw_exp_repulsion ) {
			Output::err( "SIM_CONTROL: exponential repulsion must be used in conjunction with polarvdw\n" );
			return fail;
		}
		if( sys.cdvdw_sig_repulsion) {
			Output::err("SIM_CONTROL: sig repulsion is used in conjunction with polarvdw\n");
			return fail;
		}
	}

	return ok;
}




bool SimulationControl::check_qrot_options() {

	char linebuf[maxLine];

	Output::out( "SIM_CONTROL: Quantum rotational eigenspectrum calculation enabled\n" );
	if( sys.quantum_rotation_B <= 0.0 ) {
		Output::err( "SIM_CONTROL: invalid quantum rotational constant B specified\n" );
		return fail;
	} else {
		sprintf(linebuf, "SIM_CONTROL: Quantum rotational constant B = %.3f K (%.3f cm^-1)\n", sys.quantum_rotation_B, sys.quantum_rotation_B*kB/(100.0*h*c) );
		Output::out( linebuf );
	}

	if( sys.quantum_rotation_level_max <= 0 ) {
		Output::err( "SIM_CONTROL: invalid quantum rotation level max\n" );
		return fail;
	} else {
		sprintf( linebuf, "SIM_CONTROL: Quantum rotation level max = %d\n", sys.quantum_rotation_level_max );
		Output::out( linebuf );
	}

	if( sys.quantum_rotation_l_max <= 0 ) {
		Output::err( "SIM_CONTROL: invalid quantum rotation l_max\n" );
		return fail;
	} else {
		sprintf( linebuf, "SIM_CONTROL: Quantum rotation l_max = %d\n", sys.quantum_rotation_l_max );
		Output::out( linebuf );
	}

	if(    sys.quantum_rotation_level_max   >   (sys.quantum_rotation_l_max + 1)*(sys.quantum_rotation_l_max + 1)    ) {
		Output::err( "SIM_CONTROL: quantum rotational levels cannot exceed l_max + 1 X l_max +1\n" );
		return fail;
	}

	if(  sys.quantum_rotation_sum <= 0   ||   sys.quantum_rotation_sum > sys.quantum_rotation_level_max  ) {
		Output::err( "SIM_CONTROL: quantum rotational sum for partition function invalid\n" );
		return fail;
	} else {
		sprintf( linebuf, "SIM_CONTROL: Quantum rotation sum = %d\n", sys.quantum_rotation_sum );
		Output::out( linebuf );
	}

	return ok;
}




bool SimulationControl::runSimulation() {
// Starts main run loop, according to simulation type

	char start_up_msg[maxLine] = { 0 };   // Message to display when sim starts
	char err_exit_msg[maxLine] = { 0 };   // Message to display if sim errors out
	bool (SimulationControl::*run_it)();  // run_it is a pointer to the function that we will run in order to execute the 
	                                      // simulation. The syntax is weird because this is a pointer to a member function
	                                      // rather than a C or static function. Member fxns invoked via: (this->*run_it)()

	if (mpi) {
		#ifdef _MPI
			MPI_Barrier(MPI_COMM_WORLD);
			sprintf(start_up_msg, "SIM_CONTROL: all %d cores are in sync\n", size);
			Output::out1(start_up_msg);
		#endif
	}
	


	// Set which simulation function to run, set the start-up message, and the message to display upon failure
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	switch (sys.ensemble) {
	case ENSEMBLE_PATH_INTEGRAL_NVT:
		run_it = &SimulationControl::PI_nvt_mc;
		strcat(start_up_msg, "SIM_CONTROL: *********************************************************\n");
		strcat(start_up_msg, "SIM_CONTROL: *** starting Path-Integral Monte Carlo NVT simulation ***\n");
		strcat(start_up_msg, "SIM_CONTROL: *********************************************************\n");
		strcat(err_exit_msg, "SIM_CONTROL: PIMC failed\n");
		break;

	case ENSEMBLE_NVT_GIBBS:
		run_it = &SimulationControl::Gibbs_mc;
		strcat(start_up_msg, "SIM_CONTROL: *************************************************\n");
		strcat(start_up_msg, "SIM_CONTROL: *** starting Gibbs NVT Monte Carlo simulation ***\n");
		strcat(start_up_msg, "SIM_CONTROL: *************************************************\n");
		strcat(err_exit_msg, "SIM_CONTROL: Gibbs failed\n");	
		break;

	case ENSEMBLE_SURF:
		run_it = &SimulationControl::surface;
		strcat(start_up_msg, "SIM_CONTROL: *****************************************************\n");
		strcat(start_up_msg, "SIM_CONTROL: *** starting potential energy surface calculation ***\n");
		strcat(start_up_msg, "SIM_CONTROL: *****************************************************\n");
		strcat(err_exit_msg, "SIM_CONTROL: surface module failed on error, exiting\n");
		break;

	case ENSEMBLE_SURF_FIT:
		if (sys.surf_fit_arbitrary_configs) {
			run_it = &SimulationControl::surface_fit_arbitrary;
			strcat(err_exit_msg, "SIM_CONTROL: surface fitting module (for arbitrary configurations) failed on error, exiting\n");
		}
		else {
			run_it = &SimulationControl::surface_fit;
			strcat(err_exit_msg, "SIM_CONTROL: surface fitting module failed on error, exiting\n");
		}
		strcat(start_up_msg, "SIM_CONTROL: *************************************************************\n");
		strcat(start_up_msg, "SIM_CONTROL: *** starting potential energy surface fitting calculation ***\n");
		strcat(start_up_msg, "SIM_CONTROL: *************************************************************\n");
		strcat(err_exit_msg, "SIM_CONTROL: surface fitting module (for arbitrary configurations) failed on error, exiting\n");
		break;

	case ENSEMBLE_REPLAY:  // replay trajectory and recalc energies, etc. 
		run_it = &SimulationControl::replay_trajectory;
		strcat(start_up_msg, "SIM_CONTROL: **********************************\n");
		strcat(start_up_msg, "SIM_CONTROL: *** starting trajectory replay ***\n");
		strcat(start_up_msg, "SIM_CONTROL: **********************************\n");
		strcat(err_exit_msg, "SIM_CONTROL: trajectory replay failed, exiting\n");
		break;

	case ENSEMBLE_TE:
		run_it = &SimulationControl::calculate_te;
		strcat(start_up_msg, "SIM_CONTROL: *************************************************\n");
		strcat(start_up_msg, "SIM_CONTROL: *** starting single-point energy calculation  ***\n");
		strcat(start_up_msg, "SIM_CONTROL: *************************************************\n");
		strcat(err_exit_msg, "SIM_CONTROL: single-point energy calculation failed, exiting\n");
		break;

	
	// case ENSEMBLE_UVT | ENSEMBLE_NVT | ENSEMBLE_NVE:
	default:

		run_it = &SimulationControl::mc;

		if (sys.ensemble == ENSEMBLE_UVT) {
			strcat(start_up_msg, "SIM_CONTROL: *******************************************************\n");
			strcat(start_up_msg, "SIM_CONTROL: *** starting Grand Canonical Monte Carlo simulation ***\n");
			strcat(start_up_msg, "SIM_CONTROL: *******************************************************\n");
		}
		else if (sys.ensemble == ENSEMBLE_NVT) {
			strcat(start_up_msg, "SIM_CONTROL: *******************************************************\n");
			strcat(start_up_msg, "SIM_CONTROL: ***  starting  Canonical  Monte  Carlo  simulation  ***\n");
			strcat(start_up_msg, "SIM_CONTROL: *******************************************************\n");
		}
		else if (sys.ensemble == ENSEMBLE_NVE) {
			strcat(start_up_msg, "SIM_CONTROL: *******************************************************\n");
			strcat(start_up_msg, "SIM_CONTROL: *** starting Microcanonical  Monte Carlo simulation ***\n");
			strcat(start_up_msg, "SIM_CONTROL: *******************************************************\n");
		}

		strcat(err_exit_msg, "SIM_CONTROL: MC failed on error, exiting\n");
		break;
	}



	// Run the simulation
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Output::out1(start_up_msg);
	bool sim_completed = (this->*run_it)();
	if( ! sim_completed ) {
		Output::err(err_exit_msg);
		return fail;
	}
	Output::out1("SIM_CONTROL: Simulation complete!\n");
	Output::out1("Cleaning up & exit. ");
	return ok;
}




void SimulationControl::add_orientation_site_entry(const char *id, const int site_idx ) {

	//Check to see if entry already exists
	std::map<std::string, int>::iterator it;
	it = sorbate_data_index.find(id);
	if (it == sorbate_data_index.end()) {

		// no orientation site found, create new entry
		molecular_metadata data;
		data.orientation_site = site_idx;
		data.bond_length = 0;
		data.reduced_mass = 0;
		// add data to collection and log its index
		sorbate_data.push_back(data);
		int metadata_idx = (int) sorbate_data.size() - 1;
		sorbate_data_index.insert(std::make_pair(id, metadata_idx));
	
	} else {
		sorbate_data[it->second].orientation_site = site_idx;
	}
}
int SimulationControl::get_orientation_site( std::string molecule_id) {

	std::map<std::string, int>::iterator it;
	it = sorbate_data_index.find(molecule_id);
	if (it == sorbate_data_index.end()) {
		return -1; // no orientation site found
	}

	return sorbate_data[it->second].orientation_site;
}



void SimulationControl::add_bond_length_entry(const char *id, const double bond_length) {

	//Check to see if entry already exists
	std::map<std::string, int>::iterator it;
	it = sorbate_data_index.find(id);
	if (it == sorbate_data_index.end()) {

		// no orientation site found, create new entry
		molecular_metadata data;
		data.orientation_site = -1;
		data.bond_length = bond_length;
		data.reduced_mass = 0;
		// add data to collection and log its index
		sorbate_data.push_back(data);
		int index = (int) sorbate_data.size() - 1;
		sorbate_data_index.insert(std::make_pair(id, index));

	} else {
		sorbate_data[it->second].bond_length = bond_length;
	}
}
double SimulationControl::get_bond_length(std::string molecule_id) {

	std::map<std::string, int>::iterator it;
	it = sorbate_data_index.find(molecule_id);
	if (it == sorbate_data_index.end()) {
		return 0; // no orientation site found
	}

	return sorbate_data[it->second].bond_length;
}



void SimulationControl::add_reduced_mass_entry(const char *id, const double reduced_mass) {

	//Check to see if entry already exists
	std::map<std::string, int>::iterator it;
	it = sorbate_data_index.find(id);
	if (it == sorbate_data_index.end()) {

		// no orientation site found, create new entry
		molecular_metadata data;
		data.orientation_site = -1;
		data.bond_length = 0;
		data.reduced_mass = reduced_mass;
		// add data to collection and log its index
		sorbate_data.push_back(data);
		int index = (int) sorbate_data.size() - 1;
		sorbate_data_index.insert(std::make_pair(id, index));

	}
	else {
		sorbate_data[it->second].reduced_mass = reduced_mass;
	}
}
double SimulationControl::get_reduced_mass(std::string molecule_id) {

	std::map<std::string, int>::iterator it;
	it = sorbate_data_index.find(molecule_id);
	if (it == sorbate_data_index.end()) {
		return -1.0; // no orientation site found
	}

	return sorbate_data[it->second].reduced_mass;
}