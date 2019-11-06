#include <cstring>
#include <stdio.h>

#include "Atom.h"
#include "Output.h"
#include "Pair.h"
#include "PeriodicBoundary.h"
#include "System.h"
#include "SafeOps.h"
#include "UsefulMath.h"

extern int size, rank;


static const double  ewald_alpha_default              = 0.5;
static const int     ewald_kmax_default               = 7;
static const int     ptemp_freq_default               = 20;   // default frequency for parallel tempering bath swaps
static const double  wolf_alpha_lookup_cutoff_default = 30.0; //angstroms





System::~System() {

	// ***cavity_grid
	// *ptemp
	// *polar_wolf_alpha_table
	// **A_matrix
	// **B_matrix
	// *vdw_eiso_info
	// *insertion_molecules
	// **insertion_molecules_array
	// **atom_array
	// **molecule_array
	// * molecules
	// * nodestats 
	// * avg_nodestats
	// * observables
	// * avg_observables
	// * sorbateInfo 
	// * sorbateGlobal 
	// * checkpoint;
	if( sorbateCount > 1 )
		free( sorbateGlobal );

	if( !rank ) {
		//if(fp_energy)        fclose(fp_energy);
		//if(fp_energy_csv)    fclose(fp_energy_csv);
		//if(fp_field)         fclose(fp_field);
		//if(fp_histogram)     fclose(fp_field);
		//if(fp_frozen)        fclose(fp_frozen);
		//if(fp_traj_replay)   fclose(fp_traj_replay);
		//if(fp_surf)          fclose(fp_surf);
	}

	fp_energy      = nullptr;
	fp_energy_csv  = nullptr;
	fp_field       = nullptr;
	fp_histogram   = nullptr;
	fp_frozen      = nullptr;
	fp_traj_replay = nullptr;
	fp_surf        = nullptr;

	if(nodestats)        free( nodestats );
	if(avg_nodestats)    free( avg_nodestats );
	if(observables)      free( observables );
    if(avg_observables)  free( avg_observables ); 
    if(grids) {
		if(grids->histogram)
			free( grids->histogram );
		if(grids->avg_histogram)
			free( grids->avg_histogram );
		free( grids );
	}
	if( checkpoint ) {
		if(checkpoint->observables)
			free( checkpoint->observables );
		free( checkpoint );
	}
	if(mpi_data.rcv_strct ) 
		free(mpi_data.rcv_strct);
	if( mpi_data.temperature )
		free(mpi_data.temperature);
	if( mpi_data.snd_strct)
		free(mpi_data.snd_strct);
	if( mpi_data.observables )
		free(mpi_data.observables);
	if( mpi_data.avg_nodestats)
		free(mpi_data.avg_nodestats);
	if (mpi_data.sinfo && (sorbateCount > 1))
		free(mpi_data.sinfo);
};




System::System() {
	
	// a priori system defaults

	cuda       = 0;
	opencl     = 0;

	ensemble   = 0;	
	numsteps   = 0;
	step       = 0;
	corrtime   = 0;
	ptemp_freq = 0;


	

	move_factor               = 1.0; 
	rot_factor                = 1.0; 
	scale_charge              = 1.0;
	last_volume               = 0.0;  // NPT
	volume_change_factor      = 0.25; // set default volume change factor (for NPT) to 0.25 
	adiabatic_probability     = 0.0;
	gwp_probability           = 0.0;
	insert_probability        = 0.0;
	spinflip_probability      = 0.0; 
	volume_probability        = 0.0;
	transfer_probability      = 0.0;
	bead_perturb_probability  = 0.0;
	PI_bead_perturb_factor    = 1.0;


	// io filenames
	dipole_output     [0] = 0;
	energy_output     [0] = 0;
	energy_output_csv [0] = 0;
	field_output      [0] = 0;
	frozen_output     [0] = 0;
	histogram_output  [0] = 0;
	insert_input      [0] = 0;
	pqr_input         [0] = 0;
	pqr_input_B       [0] = 0;
	pqr_output        [0] = 0;
	pqr_restart       [0] = 0;
	surf_output       [0] = 0;
	traj_input        [0] = 0;
	traj_output       [0] = 0;
	virial_output     [0] = 0;
	

	// Observables
	temperature     = 0.0; 
	pressure        = 0.0;
	free_volume     = 0.0;
	total_energy    = 0.0;	
	N               = 0.0;

	// Simulation Flags & Data
	gwp               = 0;
	long_output       = 0; // prints extended (%11.6f) coordinates
	max_bondlength    = 0; // threshold to bond (re:output files)
	parallel_restarts = 0; // is this a restart of a parallel job?


	// Cavity Stuff
	cavity_grid                  = nullptr;
	cavity_bias                  = 0; 
	cavity_grid_size             = 0;
	cavities_open                = 0;
	cavity_autoreject_repulsion  = 0.0;
	cavity_autoreject_scale      = 0.0; 
	cavity_radius                = 0.0; 
	cavity_volume                = 0.0;


	// Auto-reject Options
	//first is in terms of sigma and only applies to LJ; latter is in Angstroms and applies to all pairs
	cavity_autoreject            = 0, 
	cavity_autoreject_absolute   = 0;


	// Parallel Tempering Options
	parallel_tempering           = 0;
	max_temperature              = 0.0;
	ptemp                        = nullptr;

	
	// Info regarding  periodic boundary and unit cell geometry
	wrapall                      = 1;
	


	// RNG
	//preset_seed_on               = 0;  //for manual specification of random seed
	//preset_seed                  = 0;
	
	
	read_pqr_box_on              = 0; //read box basis from pqr
	
	// Simulated Annealing
	simulated_annealing          = 0; 
	simulated_annealing_linear   = 0;
	simulated_annealing_schedule = 0.0;
	simulated_annealing_target   = 0.0;

	// Spectre
	spectre             = 0;
	spectre_max_charge  = 0.0; 
	spectre_max_target  = 0.0;


	// Energy-corrections
	feynman_hibbs       = 0; 
	feynman_kleinert    = 0;
	feynman_hibbs_order = 0;
	vdw_fh_2be          = 0; //2BE method for polarvdw 
	rd_lrc              = 0; 
	rd_crystal          = 0;
	rd_crystal_order    = 0;

	// uVT Fugacity Functions
	h2_fugacity         = 0; 
	co2_fugacity        = 0; 
	ch4_fugacity        = 0;
	n2_fugacity         = 0;
	user_fugacities     = 0;
	fugacitiesCount     = 0;
	for( int i=0; i<maxTokens; i++ )
		fugacities[i]   = 0.0;


	// Force-field Options
	rd_only                 = 0; 
	rd_anharmonic           = 0;
	rd_anharmonic_k         = 0.0;
	rd_anharmonic_g         = 0.0;
	using_axilrod_teller    = false; 
	c6_mixing               = 0; 
	damp_dispersion         = 0; 
	using_disp_expansion    = false;
	disp_expansion_mbvdw    = 0;
	use_dreiding            = 0; 
	extrapolate_disp_coeffs = 0; 
	halgren_mixing          = 0;
	using_lj_buffered_14_7  = false;
	midzuno_kihara_approx   = 0;
	schmidt_mixing          = 0;
	force_mixing			= 0;
	bohm_ahlrichs_mixing	= 0;
	wilson_popelier_mixing	= 0;
	use_sg                  = false;
	waldmanhagler           = 0;


	// ES Options
	wolf                    = 0;
	ewald_alpha_set         = 0; 
	ewald_kmax              = 0; 
	polar_ewald_alpha_set   = 0;
	ewald_alpha             = 0;
	polar_ewald_alpha       = 0;
	
	
	// Thole Options
	polarization            = 0; 
	polarvdw                = 0;
	polarizability_tensor   = 0;
	cdvdw_exp_repulsion     = 0; 
	cdvdw_sig_repulsion     = 0;
	cdvdw_9th_repulsion     = 0;
	iterator_failed         = 0; //flag set when iterative solver fails to converge (when polar_precision is used)
	polar_iterative         = 0; 
	polar_ewald             = 0;
	polar_ewald_full        = 0;
	polar_zodid             = 0;
	polar_palmo             = 0;
	polar_rrms              = 0;
	polar_gs                = 0;
	polar_gs_ranked         = 0;
	polar_sor               = 0;
	polar_esor              = 0;
	polar_max_iter          = 0;
	polar_wolf              = 0;
	polar_wolf_full         = 0;
	polar_wolf_alpha_lookup = 0;
	polar_wolf_alpha        = 0.0; 
	polar_gamma             = 0.0;
	polar_damp              = 0.0;
	field_damp              = 0.0;
	polar_precision         = 0.0;
	polar_wolf_alpha_table         = nullptr;
	polar_wolf_alpha_lookup_cutoff = 0.0;
	polar_wolf_alpha_table_max     = 0;       //stores the total size of the array
	damp_type                      = 0;
	A_matrix                       = nullptr; // A matrix, B matrix and polarizability tensor 
	B_matrix                       = nullptr;
	for( int i=0; i<3; i++)
		for( int j=0; j<3; j++ ) {
			C_matrix[i][j]  = 0;
		}

	vdw_eiso_info = nullptr; //keeps track of molecule vdw self energies
	

	//misc
	independent_particle = 0;

	// insertions from a separate linked list
	num_insertion_molecules   = 0;         // the number of elements found in both lists below:
	insertion_molecules       = nullptr;   // linked list of molecules to be randomly inserted
	insertion_molecules_array = nullptr;   // array providing direct access to elements of above list


	// quantum rotation stuff
	quantum_rotation                  = 0; 
	quantum_rotation_hindered         = 0;
	quantum_rotation_l_max            = 0;  
	quantum_rotation_level_max        = 0; 
	quantum_rotation_phi_max          = 0;
	quantum_rotation_sum              = 0;
	quantum_rotation_theta_max        = 0;
	quantum_vibration                 = 0;
	quantum_rotation_B                = 0.0;
	quantum_rotation_hindered_barrier = 0.0;
	
	// histogram stuff
	grids            = nullptr; // Initialize Grid fields
	calc_hist        = 0;       // flag to calculate a 3D histogram 
	hist_resolution  = 0.0;
	n_histogram_bins = 0;

	//atom array	
	natoms         = 0;
	molecules      = nullptr;
	atom_array     = nullptr;
	molecule_array = nullptr;

	//replay option
	calc_pressure    = 0;
	calc_pressure_dv = 0.0;

	// Linked list head that will keep track of separate average-observables for each sorbate in the system.
	sorbateCount     = 0;        // number of sorbates in the system.
	sorbateInfo      = nullptr;  // stores an array of sorbate Info
	sorbateInsert    = 0;        // which sorbate was last inserted
	sorbateGlobal    = nullptr;  // where the global average is stored

	fp_energy        = nullptr;
	fp_energy_csv    = nullptr;
	fp_field         = nullptr;
	fp_histogram     = nullptr;
	fp_frozen        = nullptr;
	fp_traj_replay   = nullptr;
	fp_surf          = nullptr;


	ee_local                      = 0;   // Exhaustive Enumeration option
	fit_best_square_error         = 0.0;
	fit_max_energy                = 0.0;
	fit_schedule                  = 0.0;
	fit_start_temp                = 0.0;
	fit_boltzmann_weight          = 0;
	surf_decomp                   = 0;
	surf_fit_arbitrary_configs    = 0;
	surf_qshift_on                = 0;
	surf_weight_constant_on       = 0;
	surf_global_axis_on           = 0;
	surf_ang                      = 0.0;
	surf_inc                      = 0.0;
	surf_min                      = 0.0;
	surf_max                      = 0.0;
	surf_descent                  = 0;
	surf_preserve                 = 0;
	surf_preserve_rotation_on     = 0;
	surf_preserve_rotation_alpha1 = 0.0;
	surf_preserve_rotation_alpha2 = 0.0;
	surf_preserve_rotation_beta1  = 0.0;
	surf_preserve_rotation_beta2  = 0.0;
	surf_preserve_rotation_gamma1 = 0.0;
	surf_preserve_rotation_gamma2 = 0.0;
	range_eps                     = 0.0; 
	range_sig                     = 0.0; 
	step_eps                      = 0.0;  
	step_sig                      = 0.0;
	surf_print_level              = 0;   // sets the amount of output (1-6) that correspond to the nested loops in surface.c
	surf_scale_alpha_on           = 0; 
	surf_scale_epsilon_on         = 0;
	surf_scale_omega_on           = 0;
	surf_scale_sigma_on           = 0;
	surf_scale_pol_on             = 0;
	surf_scale_q_on               = 0;
	surf_scale_r_on               = 0;
	surf_scale_c6_on              = 0;
	surf_scale_c8_on              = 0;
	surf_scale_c10_on             = 0;
	surf_scale_epsilon            = 0.0; 
	surf_scale_r                  = 0.0;
	surf_scale_omega              = 0.0;
	surf_scale_sigma              = 0.0;
	surf_scale_q                  = 0.0;
	surf_scale_alpha              = 0.0;
	surf_scale_pol                = 0.0;
	surf_scale_c6                 = 0.0;
	surf_scale_c8                 = 0.0;
	surf_scale_c10                = 0.0;
	surf_virial                   = 0;
	surf_quadrupole               = 0.0;
	surf_weight_constant          = 0.0;

	// default ewald parameters 
	ewald_alpha                    = ewald_alpha_default;
	ewald_kmax                     = ewald_kmax_default;
	polar_ewald_alpha              = ewald_alpha_default;
	polar_wolf_alpha_lookup_cutoff = wolf_alpha_lookup_cutoff_default; 

	// default polarization parameters 
	polar_gamma = 1.0;
	polar_max_iter = 10;

	// default rd LRC flag 
	rd_lrc = 1;

	// Initialize fit_input_list to reflect an empty list
	fit_input_list.next = nullptr;
	fit_input_list.data.count = 0;

	//set default jobname
	sprintf( job_name, "untitled");

	// Allocate space for statistics
	nodestats       = nullptr;
	avg_nodestats   = nullptr;
	observables     = nullptr;
    avg_observables = nullptr;
    checkpoint      = nullptr;

	// Mpi structs
	mpi_data.msgsize = 0;
	mpi_data.observables   = nullptr;
	mpi_data.avg_nodestats = nullptr;
	mpi_data.sinfo         = nullptr;
	mpi_data.temperature   = nullptr;
	mpi_data.snd_strct     = nullptr;
	mpi_data.rcv_strct     = nullptr;

    
}




System::System( const System &sd ) {
// This creates a shallow copy, if you need a deep copy you will have to do it manually
// io filenames

	for( int i=0; i<maxLine; i++ ) {
		job_name          [ i ] = sd.job_name          [ i ];
		dipole_output     [ i ] = sd.dipole_output     [ i ];
		energy_output     [ i ] = sd.energy_output     [ i ];
		energy_output_csv [ i ] = sd.energy_output_csv [ i ];
		field_output      [ i ] = sd.field_output      [ i ];
		frozen_output     [ i ] = sd.frozen_output     [ i ];
		histogram_output  [ i ] = sd.histogram_output  [ i ];
		insert_input      [ i ] = sd.insert_input      [ i ];
		pqr_input         [ i ] = sd.pqr_input         [ i ];
		pqr_input_B       [ i ] = sd.pqr_input_B       [ i ];
		pqr_output        [ i ] = sd.pqr_output        [ i ];
		pqr_restart       [ i ] = sd.pqr_restart       [ i ];
		surf_output       [ i ] = sd.surf_output       [ i ];
		traj_input        [ i ] = sd.traj_input        [ i ];
		traj_output       [ i ] = sd.traj_output       [ i ];
		virial_output     [ i ] = sd.virial_output     [ i ];
	}

	// Compilation type flags and data
	cuda                          = sd.cuda;
	opencl                        = sd.opencl;
	#ifdef OPENCL
		ocl                       = sd.ocl;
	#endif


	// Monte Carlo Controls
	ensemble                      = sd.ensemble;
	numsteps                      = sd.numsteps;
	step                          = sd.step;
	corrtime                      = sd.corrtime;
	ptemp_freq                    = sd.ptemp_freq;
	move_factor                   = sd.move_factor;
	rot_factor                    = sd.rot_factor;
	last_volume                   = sd.last_volume; // NPT
	volume_change_factor          = sd.volume_change_factor; // NPT
	adiabatic_probability         = sd.adiabatic_probability;
	gwp_probability               = sd.gwp_probability;
	insert_probability            = sd.insert_probability;
	spinflip_probability          = sd.spinflip_probability;
	volume_probability            = sd.volume_probability;
	transfer_probability          = sd.transfer_probability;
	bead_perturb_probability      = sd.bead_perturb_probability;
	PI_bead_perturb_factor        = sd.PI_bead_perturb_factor;

	// Observables
	temperature                   = sd.temperature;
	pressure                      = sd.pressure;
	free_volume                   = sd.free_volume;
	total_energy                  = sd.total_energy;
	N                             = sd.N;
	fugacitiesCount               = sd.fugacitiesCount;
	for(int i=0;i<maxTokens;i++)
		fugacities[ i ] = sd.fugacities[ i ];
	
	// Simulation Flags & Data
	gwp                           = sd.gwp;
	long_output                   = sd.long_output; // prints extended (%11.6f) coordinates
	max_bondlength                = sd.max_bondlength; // threshold to bond (re:output files)
	parallel_restarts             = sd.parallel_restarts; // is this a restart of a parallel job?

	// Cavity Stuff
//	cavity_t ***cavity_grid;
	cavity_bias                   = sd.cavity_bias;
	cavity_grid_size              = sd.cavity_grid_size;
	cavities_open                 = sd.cavities_open;
	cavity_autoreject_repulsion   = sd.cavity_autoreject_repulsion;
	cavity_autoreject_scale       = sd.cavity_autoreject_scale;
	cavity_radius                 = sd.cavity_radius;
	cavity_volume                 = sd.cavity_volume;

	// Auto-reject Options
	//first is in terms of sigma and only applies to LJ; latter is in Angstroms and applies to all pairs
	cavity_autoreject             = sd.cavity_autoreject; 
	cavity_autoreject_absolute    = sd.cavity_autoreject_absolute;

	// Parallel Tempering Options
	parallel_tempering            = sd.parallel_tempering;
	max_temperature               = sd.max_temperature;
//	ptemp_t *ptemp;
	
	// Info regarding  periodic boundary and unit cell geometry
	wrapall                       = sd.wrapall;
	pbc.cutoff                    = sd.pbc.cutoff; // radial cutoff (A)
	pbc.volume                    = sd.pbc.volume; // unit cell volume (A^3) 

	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) {
			C_matrix            [ i ][ j ] = sd.C_matrix            [ i ][ j ];
			pbc.basis           [ i ][ j ] = sd.pbc.basis           [ i ][ j ]; // unit cell lattice (A)
			pbc.reciprocal_basis[ i ][ j ] = sd.pbc.reciprocal_basis[ i ][ j ]; // reciprocal space lattice (1/A)
		}
	}
	
	// RNG
	preset_seed_on                = sd.preset_seed_on; //for manual specification of random seeds
	preset_seed                   = sd.preset_seed;
	

	read_pqr_box_on               = sd.read_pqr_box_on; //read box basis from pqr
	
	// Simulated Annealing
	simulated_annealing           = sd.simulated_annealing;
	simulated_annealing_linear    = sd.simulated_annealing_linear;
	simulated_annealing_schedule  = sd.simulated_annealing_schedule;
	simulated_annealing_target    = sd.simulated_annealing_target;

	// Spectre
	spectre                       = sd.spectre;
	spectre_max_charge            = sd.spectre_max_charge;
	spectre_max_target            = sd.spectre_max_target;

	// Energy-corrections
	feynman_hibbs                 = sd.feynman_hibbs; 
	feynman_kleinert              = sd.feynman_kleinert;
	feynman_hibbs_order           = sd.feynman_hibbs_order;
	vdw_fh_2be                    = sd.vdw_fh_2be; //2BE method for polarvdw 
	rd_lrc                        = sd.rd_lrc;
	rd_crystal                    = sd.rd_crystal;
	rd_crystal_order              = sd.rd_crystal_order;

	// uVT Fugacity Functions
	h2_fugacity                   = sd.h2_fugacity; 
	co2_fugacity                  = sd.co2_fugacity;
	ch4_fugacity                  = sd.ch4_fugacity;
	n2_fugacity                   = sd.n2_fugacity;
	user_fugacities               = sd.user_fugacities;

	// Force-field Options
	rd_only                       = sd.rd_only; 
	rd_anharmonic                 = sd.rd_anharmonic;
	rd_anharmonic_k               = sd.rd_anharmonic_k; 
	rd_anharmonic_g               = sd.rd_anharmonic_g;
	using_axilrod_teller          = sd.using_axilrod_teller;
	c6_mixing                     = sd.c6_mixing; 
	damp_dispersion               = sd.damp_dispersion;
	using_disp_expansion          = sd.using_disp_expansion;
	disp_expansion_mbvdw          = sd.disp_expansion_mbvdw;
	use_dreiding                  = sd.use_dreiding;
	extrapolate_disp_coeffs       = sd.extrapolate_disp_coeffs;
	halgren_mixing                = sd.halgren_mixing;
	using_lj_buffered_14_7        = sd.using_lj_buffered_14_7;
	midzuno_kihara_approx         = sd.midzuno_kihara_approx;
	schmidt_mixing                = sd.schmidt_mixing;
	force_mixing				  = sd.force_mixing;
	bohm_ahlrichs_mixing		  = sd.bohm_ahlrichs_mixing;
	wilson_popelier_mixing		  = sd.wilson_popelier_mixing;
	use_sg                        = sd.use_sg;
	waldmanhagler                 = sd.waldmanhagler;

	// ES Options
	wolf                          = sd.wolf;
	ewald_alpha_set               = sd.ewald_alpha_set;
	ewald_kmax                    = sd.ewald_kmax;
	polar_ewald_alpha_set         = sd.polar_ewald_alpha_set;
	ewald_alpha                   = sd.ewald_alpha;
	polar_ewald_alpha             = sd.polar_ewald_alpha;
	
	// Thole Options
	polarization                  = sd.polarization; 
	polarvdw                      = sd.polarvdw;
	polarizability_tensor         = sd.polarizability_tensor;
	cdvdw_exp_repulsion           = sd.cdvdw_exp_repulsion;
	cdvdw_sig_repulsion           = sd.cdvdw_sig_repulsion;
	cdvdw_9th_repulsion           = sd.cdvdw_9th_repulsion;
	iterator_failed               = sd.iterator_failed; 
	polar_iterative               = sd.polar_iterative;
	polar_ewald                   = sd.polar_ewald;
	polar_ewald_full              = sd.polar_ewald_full;
	polar_zodid                   = sd.polar_zodid;
	polar_palmo                   = sd.polar_palmo;
	polar_rrms                    = sd.polar_rrms;
	polar_gs                      = sd.polar_gs;
	polar_gs_ranked               = sd.polar_gs_ranked;
	polar_sor                     = sd.polar_sor;
	polar_esor                    = sd.polar_esor;
	polar_max_iter                = sd.polar_max_iter;
	polar_wolf                    = sd.polar_wolf;
	polar_wolf_full               = sd.polar_wolf_full;
	polar_wolf_alpha_lookup       = sd.polar_wolf_alpha_lookup;
	polar_wolf_alpha              = sd.polar_wolf_alpha;
	polar_gamma                   = sd.polar_gamma;
	polar_damp                    = sd.polar_damp;
	field_damp                    = sd.field_damp;
	polar_precision               = sd.polar_precision;
	polar_wolf_alpha_lookup_cutoff= sd.polar_wolf_alpha_lookup_cutoff;
	polar_wolf_alpha_table_max    = sd.polar_wolf_alpha_table_max; //stores the total size of the array
	damp_type                     = sd.damp_type;
	

	//misc
	scale_charge                  = sd.scale_charge;
	independent_particle          = sd.independent_particle;

	// insertions from a separate linked list
	num_insertion_molecules       = sd.num_insertion_molecules; // the number of elements found in both lists below:

	// quantum rotation stuff
	quantum_rotation                  = sd.quantum_rotation;
	quantum_rotation_hindered         = sd.quantum_rotation_hindered;
	quantum_rotation_l_max            = sd.quantum_rotation_l_max;
	quantum_rotation_level_max        = sd.quantum_rotation_level_max;
	quantum_rotation_phi_max          = sd.quantum_rotation_phi_max;
	quantum_rotation_sum              = sd.quantum_rotation_sum;
	quantum_rotation_theta_max        = sd.quantum_rotation_theta_max;
	quantum_vibration                 = sd.quantum_vibration;
	quantum_rotation_B                = sd.quantum_rotation_B;
	quantum_rotation_hindered_barrier = sd.quantum_rotation_hindered_barrier;
	
	// histogram stuff
	grids                         = nullptr;
	calc_hist                     = sd.calc_hist; // flag to calculate a 3D histogram 
	hist_resolution               = sd.hist_resolution;
	n_histogram_bins              = sd.n_histogram_bins;

	//atom array	
	natoms                        = sd.natoms;
	
	//replay option
	calc_pressure                 = sd.calc_pressure;
	calc_pressure_dv              = sd.calc_pressure_dv;
	
	// Linked list head that will keep track of separate average-observables
	// for each sorbate in the system.
	sorbateCount                  = sorbateCount;  // Number of sorbates in the system.
	sorbateInsert                 = sd.sorbateInsert; //which sorbate was last inserted


	checkpoint                    = sd.checkpoint;
	fit_input_list.data.count     = 0;
	fit_input_list.next           = nullptr;
	ee_local                      = sd.ee_local;              // Exhaustive Enumeration option
	fit_best_square_error         = sd.fit_best_square_error;
	fit_max_energy                = sd.fit_max_energy;	       
	fit_schedule                  = sd.fit_schedule;
	fit_start_temp                = sd.fit_start_temp;
	fit_boltzmann_weight          = sd.fit_boltzmann_weight;
	surf_decomp                   = sd.surf_decomp;
	surf_fit_arbitrary_configs    = sd.surf_fit_arbitrary_configs;
	surf_qshift_on                = sd.surf_qshift_on;
	surf_weight_constant_on       = sd.surf_weight_constant_on;
	surf_global_axis_on           = sd.surf_global_axis_on;
	surf_ang                      = sd.surf_ang;
	surf_inc                      = sd.surf_inc;
	surf_min                      = sd.surf_min;
	surf_max                      = sd.surf_max;
	surf_descent                  = sd.surf_descent;
	surf_preserve                 = sd.surf_preserve;
	surf_preserve_rotation_on     = sd.surf_preserve_rotation_on;
	surf_preserve_rotation_alpha1 = sd.surf_preserve_rotation_alpha1;
	surf_preserve_rotation_alpha2 = sd.surf_preserve_rotation_alpha2;
	surf_preserve_rotation_beta1  = sd.surf_preserve_rotation_beta1;
	surf_preserve_rotation_beta2  = sd.surf_preserve_rotation_beta2;
	surf_preserve_rotation_gamma1 = sd.surf_preserve_rotation_gamma1;
	surf_preserve_rotation_gamma2 = sd.surf_preserve_rotation_gamma2;
	range_eps                     = sd.range_eps;
	range_sig                     = sd.range_sig;
	step_eps                      = sd.step_eps;
	step_sig                      = sd.step_sig;
	surf_print_level              = sd.surf_print_level; // sets the amount of output (1-6) that correspond to the nested loops in surface.c
	surf_scale_alpha_on           = sd.surf_scale_alpha_on; 
	surf_scale_epsilon_on         = sd.surf_scale_epsilon_on;
	surf_scale_omega_on           = sd.surf_scale_omega_on;
	surf_scale_sigma_on           = sd.surf_scale_sigma_on;
	surf_scale_pol_on             = sd.surf_scale_pol_on;
	surf_scale_q_on               = sd.surf_scale_q_on;
	surf_scale_r_on               = sd.surf_scale_r_on;
	surf_scale_c6_on              = sd.surf_scale_c6_on;
	surf_scale_c8_on              = sd.surf_scale_c8_on;
	surf_scale_c10_on             = sd.surf_scale_c10_on;
	surf_scale_epsilon            = sd.surf_scale_epsilon; 
	surf_scale_r                  = sd.surf_scale_r;
	surf_scale_omega              = sd.surf_scale_omega;
	surf_scale_sigma              = sd.surf_scale_sigma;
	surf_scale_q                  = sd.surf_scale_q;
	surf_scale_alpha              = sd.surf_scale_alpha;
	surf_scale_pol                = sd.surf_scale_pol;
	surf_scale_c6                 = sd.surf_scale_c6;
	surf_scale_c8                 = sd.surf_scale_c8;
	surf_scale_c10                = sd.surf_scale_c10;
	surf_virial                   = sd.surf_virial;
	surf_quadrupole               = sd.surf_quadrupole;
	surf_weight_constant          = sd.surf_weight_constant;

	// The following structures are not copied
	atom_array                    = nullptr;
	molecule_array                = nullptr;
	molecules                     = nullptr;
	polar_wolf_alpha_table        = nullptr;
	A_matrix                      = nullptr;
	B_matrix                      = nullptr;
	vdw_eiso_info                 = nullptr;
	insertion_molecules           = nullptr;
	insertion_molecules_array     = nullptr;

	nodestats                     = nullptr;
	avg_nodestats                 = nullptr;
	observables                   = nullptr;
	avg_observables               = nullptr;
	sorbateInfo                   = nullptr;
	sorbateGlobal                 = nullptr;

	fp_energy                     = nullptr;
	fp_energy_csv                 = nullptr;
	fp_field                      = nullptr;
	fp_histogram                  = nullptr;
	fp_frozen                     = nullptr;
	fp_traj_replay                = nullptr;
	fp_surf                       = nullptr;

	// Mpi structs
	mpi_data.msgsize = 0;
	mpi_data.observables          = nullptr;
	mpi_data.avg_nodestats        = nullptr;
	mpi_data.sinfo                = nullptr;
	mpi_data.temperature          = nullptr;
	mpi_data.snd_strct            = nullptr;
	mpi_data.rcv_strct            = nullptr;

}




bool System::setup_simulation_box() {
// reads molecular input (geometry) file and computes periodic boundary conditions
// If a molecule insertion list is provided, it is created

	char *input_file;

	if( ensemble == ENSEMBLE_REPLAY ) 
		input_file = traj_input;		
	else
		input_file = pqr_input;
	
	
		
	// read in input.pqr molecules
	read_molecules( input_file );
	
	// Were molecules read into the simulation?
	if(   ! molecules  ) {
		if( ensemble == ENSEMBLE_REPLAY ) {
			Output::out( "SYSTEM: end of trajectory file\n" );
			return ok;
		} else {
			Output::err( "SYSTEM: error reading in input molecules. No molecules found in system.\n" );
			return fail;
		}
	}
	else 
		Output::out( "SYSTEM: finished reading in molecules\n" );


	// read in pqr box and calculate periodic boundary conditions if neccessary
	if( read_pqr_box_on ) 
		read_pqr_box( input_file );

	
	update_pbc();
	if ( ensemble != ENSEMBLE_SURF && ensemble != ENSEMBLE_SURF_FIT ) {
		if(  (pbc.volume <= 0.0)   ||   (pbc.cutoff <= 0.0)  ) {
			Output::err("SYSTEM: invalid simulation box dimensions.\n");
			throw invalid_box_dimensions;
		}
	}
	if( pbc.volume > 0 )
		pbc.printboxdim();

	
	
	
	// read in the insertion molecules
	if( insert_input[0] ) {
		/*
		// Multi-sorbate stuffs - problem for later-brant ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		insertion_molecules = read_insertion_molecules(system);
		if( ! insertion_molecules ) {
			Output::err("SYSTEM: error read in insertion molecules\n");
			return fail;
		} else
			Output::out("SYSTEM: finished reading in insertion molecules\n");
		*/
	}
	else //else only 1 sorbate type
		sorbateCount = 1;

	// now that we've read in the sorbates, we can check that user_fugacities is properly set (if used)
	if( user_fugacities ) {
		if( fugacitiesCount != sorbateCount ) {
			Output::err("SYSTEM: number of fugacities set via user_fugacities does not match the number of sorbates.\n");
			return fail;
		}
	}
	
	return ok;
}




void System::read_molecules( char *input_file ) {
// reads the molecular geometry file and allocates/populates the data structures
// that will represent the system. 
	FILE   *fp = SafeOps::openFile( input_file, "r", __LINE__, __FILE__ );
	read_molecules(fp);
	fclose(fp);

}
void System::read_molecules( FILE *fp ) {
// reads the molecular geometry file and allocates/populates the data structures
// that will represent the system. 

	fpos_t file_pos;
	fgetpos( fp, &file_pos ); //get file pointer position, we will restore this when done
	
	Atom *atom_ptr         = nullptr;
	Atom *prev_atom_ptr    = nullptr;
	Molecule *molecule_ptr = nullptr;
	
	char   linebuf           [maxLine] = {0},   err_msg     [maxLine] = {0};
	char   token_atom        [maxLine] = {0},   token_atomid[maxLine] = {0},   token_atomtype  [maxLine] = {0},
	       token_moleculetype[maxLine] = {0},   token_frozen[maxLine] = {0},   token_moleculeid[maxLine] = {0},
	       token_x           [maxLine] = {0},   token_y     [maxLine] = {0},   token_z         [maxLine] = {0},
	       token_mass        [maxLine] = {0},   token_charge[maxLine] = {0},   token_alpha     [maxLine] = {0},
	       token_epsilon     [maxLine] = {0},   token_sigma [maxLine] = {0},   token_omega     [maxLine] = {0},
	       token_gwp_alpha   [maxLine] = {0},   token_c6    [maxLine] = {0},   token_c8        [maxLine] = {0},
	       token_c10         [maxLine] = {0},   token_c9    [maxLine] = {0};

	int    current_frozen     = 0,  current_adiabatic = 0,  current_spectre       = 0,  current_target = 0,
	       current_moleculeid = 0,  current_atomid    = 0,  current_site_neighbor = 0,   
	       moveable           = 0,  spectres          = 0,  targets               = 0,  atom_counter   = 0;

	double current_x       = 0,  current_y      = 0,  current_z     = 0,
	       current_mass    = 0,  current_charge = 0,  current_alpha = 0,
	       current_epsilon = 0,  current_sigma  = 0,  current_omega = 0,  current_gwp_alpha = 0,
	       current_c6      = 0,  current_c8     = 0,  current_c10   = 0,  current_c9        = 0;

	
	// allocate the start of the list 
	SafeOps::calloc(molecules, 1, sizeof(Molecule), __LINE__, __FILE__ );
	molecule_ptr = molecules;

	molecule_ptr->id = 1;
	SafeOps::calloc( molecule_ptr->atoms, 1, sizeof(Atom), __LINE__, __FILE__ );
	atom_ptr = molecule_ptr->atoms;
	prev_atom_ptr = atom_ptr;

	// clear the linebuffer and read the tokens in 
	atom_counter = 0;
	std::memset(linebuf, 0, maxLine);

	while(fgets(linebuf,maxLine,fp)) {
		// clear the tokens 
		std::memset(token_atom,         0, maxLine);
		std::memset(token_atomid,       0, maxLine);
		std::memset(token_atomtype,     0, maxLine);
		std::memset(token_moleculetype, 0, maxLine);
		std::memset(token_frozen,       0, maxLine);
		std::memset(token_moleculeid,   0, maxLine);
		std::memset(token_x,            0, maxLine);
		std::memset(token_y,            0, maxLine);
		std::memset(token_z,            0, maxLine);
		std::memset(token_mass,         0, maxLine);
		std::memset(token_charge,       0, maxLine);
		std::memset(token_alpha,        0, maxLine);
		std::memset(token_epsilon,      0, maxLine);
		std::memset(token_sigma,        0, maxLine);
		std::memset(token_omega,        0, maxLine);
		std::memset(token_gwp_alpha,    0, maxLine);
		std::memset(token_c6,           0, maxLine);
		std::memset(token_c8,           0, maxLine);
		std::memset(token_c10,          0, maxLine);
		std::memset(token_c9,           0, maxLine);

		// parse the line 
		sscanf(linebuf, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", 
			token_atom, token_atomid, token_atomtype, token_moleculetype, token_frozen, 
			token_moleculeid, token_x, token_y, token_z, token_mass, token_charge, 
			token_alpha, token_epsilon, token_sigma, token_omega, token_gwp_alpha,
			token_c6, token_c8, token_c10, token_c9);

		if( SafeOps::strncasecmp(token_atom, "END", 3)) break; //we've reached the end of the current molecule, quit looping

		if(   SafeOps::iequals(token_atom,"ATOM")   &&   ! SafeOps::iequals(token_moleculetype,"BOX")   ) {

			current_frozen    = 0; 
			current_adiabatic = 0;
			current_spectre   = 0;
			current_target    = 0;

			if(SafeOps::iequals(token_frozen, "F"))
				current_frozen = 1;
			if(SafeOps::iequals(token_frozen, "A"))
				current_adiabatic = 1;
			if(SafeOps::iequals(token_frozen, "S"))
				current_spectre = 1;
			if(SafeOps::iequals(token_frozen, "T"))
				current_target = 1;

			bool flag =  SafeOps::atoi( token_moleculeid, current_moleculeid ) 
			          && SafeOps::atoi( token_atomid,     current_atomid     )
			          && SafeOps::atod( token_x,          current_x          ) 
			          && SafeOps::atod( token_y,          current_y          ) 
			          && SafeOps::atod( token_z,          current_z          )
			          && SafeOps::atod( token_mass,       current_mass       )	// mass in amu
			          && SafeOps::atod( token_alpha,      current_alpha      )
			          && SafeOps::atod( token_epsilon,    current_epsilon    )
			          && SafeOps::atod( token_sigma,      current_sigma      )
			          && SafeOps::atod( token_omega,      current_omega      )
			          && SafeOps::atod( token_gwp_alpha,  current_gwp_alpha  )
			          && SafeOps::atod( token_charge,     current_charge     );
			if(  ! flag  ) {
				Output::out(linebuf);
				throw invalid_datum;
			}
			current_charge *= E2REDUCED;  // convert charge into reduced units

			
			

			// if the following tokens were picked up, attempt to convert them:
			if(   token_c6 [0]   &&   ! SafeOps::atod( token_c6,  current_c6 )   ) {
				sprintf( err_msg, "ERROR: Invalid data token \"%s\" found in geometry file:\nERROR: %s\n", token_c6, linebuf );
				Output::err(err_msg);
				throw invalid_datum;
			}

			if(   token_c8 [0]   &&   ! SafeOps::atod( token_c8,  current_c8 )   ) {
				sprintf( err_msg, "ERROR: Invalid data token \"%s\" found in geometry file:\nERROR: %s\n", token_c8, linebuf );
				Output::err(err_msg);
				throw invalid_datum;
			}

			if(   token_c10[0]   &&   ! SafeOps::atod( token_c10, current_c10)   ) {
				sprintf( err_msg, "ERROR: Invalid data token \"%s\" found in geometry file:\nERROR: %s\n", token_c10, linebuf );
				Output::err(err_msg);
				throw invalid_datum;
			}

			if(   token_c9 [0]   &&   ! SafeOps::atod( token_c9,  current_c9 )   ) {
				sprintf( err_msg, "ERROR: Invalid data token \"%s\" found in geometry file:\nERROR: %s\n", token_c9, linebuf );
				Output::err(err_msg);
				throw invalid_datum;
			}
			


			if ( cdvdw_sig_repulsion ) {
				if ( current_epsilon != 1.0 ) {
					Output::err("warning: setting epsilon to 1.0 (due to sig_repulsion)\n");
					current_epsilon = 1.0;
				}
			}
			else if ( polarvdw && !cdvdw_exp_repulsion ) {
				if ( current_sigma != 1.0 ) {
					Output::err("warning: setting sigma to 1.0 (due to polarvdw)\n");
					current_sigma = 1.0;
				}
			}
			// Functionality of site_neighbor disabled in favor of omega/gwp_alpha parameters
			// Current behavior is to default to atom 0, typically the center of mass for
			// the molecule.
			current_site_neighbor = 0; //atoi( token_site_neighbor );
                        
			if( current_frozen )
				current_charge *= scale_charge;

			if( molecule_ptr->id != current_moleculeid ) {
				SafeOps::calloc( molecule_ptr->next, 1, sizeof(Molecule), __LINE__, __FILE__ );
				molecule_ptr = molecule_ptr->next;
				SafeOps::calloc(molecule_ptr->atoms, 1, sizeof(Atom), __LINE__, __FILE__ );
				prev_atom_ptr->next = nullptr;
				free(atom_ptr);
				atom_ptr = molecule_ptr->atoms;
			}

			strcpy(molecule_ptr->moleculetype, token_moleculetype);

			molecule_ptr->id        = current_moleculeid;
			molecule_ptr->frozen    = current_frozen;
			molecule_ptr->adiabatic = current_adiabatic;
			molecule_ptr->spectre   = current_spectre;
			molecule_ptr->target    = current_target;
			molecule_ptr->mass     += current_mass;

			#ifdef QM_ROTATION
				// if quantum rot calc. enabled, allocate the necessary structures
				if( quantum_rotation && !molecule_ptr->frozen && !molecule_ptr->quantum_rotational_energies )
					allocqmrotation(system,molecule_ptr);
			#endif // QM_ROTATION

			#ifdef XXX
				// if quantum vib calc. enabled, allocate the necessary structures 
				if( quantum_vibration && !molecule_ptr->frozen )
					allocqmvibration(system,molecule_ptr);
			#endif // XXX

			++atom_counter;
			atom_ptr->id             = atom_counter;
			atom_ptr->bond_id        = current_atomid;
			atom_ptr->frozen         = current_frozen;
			atom_ptr->adiabatic      = current_adiabatic;
			atom_ptr->spectre        = current_spectre;
			atom_ptr->target         = current_target;
			atom_ptr->pos[0]         = current_x;
			atom_ptr->pos[1]         = current_y;
			atom_ptr->pos[2]         = current_z;
			atom_ptr->mass           = current_mass;
			atom_ptr->charge         = current_charge;
			atom_ptr->polarizability = current_alpha;
			atom_ptr->epsilon        = current_epsilon;
			atom_ptr->sigma          = current_sigma;
			atom_ptr->omega          = current_omega;
			atom_ptr->gwp_alpha      = current_gwp_alpha;
			atom_ptr->c6             = current_c6;
			atom_ptr->c8             = current_c8;
			atom_ptr->c10            = current_c10;
			atom_ptr->c9             = current_c9;
			std::memset( atom_ptr->atomtype, 0, maxLine );
			strcpy( atom_ptr->atomtype, token_atomtype );
			if(current_gwp_alpha != 0.)
				atom_ptr->gwp_spin = 1;
			else
				atom_ptr->gwp_spin = 0;

			atom_ptr->site_neighbor_id = current_site_neighbor;
			SafeOps::calloc( atom_ptr->next, 1, sizeof(Atom), __LINE__, __FILE__ );
			prev_atom_ptr = atom_ptr;
			atom_ptr      = atom_ptr->next;
		}

		std::memset(linebuf, 0, maxLine);
	}

	// terminate the atom list 
	prev_atom_ptr->next = nullptr;
	free( atom_ptr );

	// scan the list, make sure that there is at least one moveable molecule 
	moveable = 0;
	spectres = 0;
	targets  = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		if(!molecule_ptr->frozen ) ++moveable;
		if(molecule_ptr ->target ) ++targets;
		if(molecule_ptr ->spectre) ++spectres;
	}

	if(spectre) {
		if(!spectres || !targets) {
			Output::err("SYSTEM: either no targets or spectres found\n");
			throw missing_required_datum;
		}
	} else {
		if(!moveable) {
			Output::err("SYSTEM: no moveable molecules found, there must be at least one in your PQR file\n");
			throw missing_required_datum;
		}
	}
	if (!atom_counter) {
		free(molecules);
		free(molecule_ptr);
		throw molecule_wo_atoms;
	}

	// restore file position--some routines require this
	fsetpos( fp, &file_pos );
}




void System::read_pqr_box( char * input_file ) {

	FILE *fp = SafeOps::openFile( input_file, "r", __LINE__, __FILE__ );

	char buffer[maxLine], token[7][maxLine],
	     msg   [maxLine];
	int  basis_set[3] = { 0, 0, 0 };  //flags to make sure we set all basis vectors

	Output::out("SYSTEM: (read_pqr_box) checking input pqr for basis info\n");



	while(  fgets(buffer, maxLine, fp)  !=  NULL  ) 
	{

		if( basis_set[0]  &&  basis_set[1]  &&  basis_set[2] )
			break; // All basis vectors set

		sscanf(buffer, "%s %s %s %s %s %s %s", 
			token[0], token[1], token[2], token[3], token[4], token[5], token[6]);

		if( (!strncmp(token[0],"END",3)) ) break; //if end of molecule, then stop searching

		if(  (!strcmp(token[0],"REMARK"))   &&   (!strcmp(token[1],"BOX"))   &&   (!strcmp(token[3],"="))  )
		{
			if (!strcmp(token[2],"BASIS[0]")) { //set basis[0]
				if(  ! SafeOps::atod( token[4], pbc.basis[0][0] )  )    continue; //make sure each conversion is successful
				if(  ! SafeOps::atod( token[5], pbc.basis[0][1] )  )    continue;
				if(  ! SafeOps::atod( token[6], pbc.basis[0][2] )  )    continue;
				//if we get this far, then we've successfully read in the basis vector
				basis_set[0] = 1;
			}
			if(!strcmp(token[2],"BASIS[1]")) { //set basis[0]
				if(  ! SafeOps::atod( token[4], pbc.basis[1][0] )  )    continue;  //make sure each conversion is successful
				if(  ! SafeOps::atod( token[5], pbc.basis[1][1] )  )    continue; 
				if(  ! SafeOps::atod( token[6], pbc.basis[1][2] )  )    continue; 
				//if we get this far, then we've successfully read in the basis vector
				basis_set[1] = 1;
			}
			if(!strcmp(token[2],"BASIS[2]")) { //set basis[0]
				if(  ! SafeOps::atod( token[4], pbc.basis[2][0] )  )    continue;  //make sure each conversion is successful
				if(  ! SafeOps::atod( token[5], pbc.basis[2][1] )  )    continue;
				if(  ! SafeOps::atod( token[6], pbc.basis[2][2] )  )    continue;
				//if we get this far, then we've successfully read in the basis vector
				basis_set[2] = 1;
			}
			else continue;
		}
		else continue;
	}


	if( basis_set[0] ) {
		sprintf( msg, "SYSTEM: basis[0] successfully read from pqr {%.5lf %.5lf %.5lf}\n", pbc.basis[0][0], pbc.basis[0][1], pbc.basis[0][2]);
		Output::out(msg);
	} else {
		sprintf( msg, "SYSTEM: unable to read basis[0] from pqr file.\n");
		Output::err(msg);
	}


	if( basis_set[1] ) {
		sprintf( msg, "SYSTEM: basis[1] successfully read from pqr {%.5lf %.5lf %.5lf}\n", pbc.basis[1][0], pbc.basis[1][1], pbc.basis[1][2]);
		Output::out(msg);
	} else {
		sprintf(msg, "SYSTEM: unable to read basis[1] from pqr file.\n");
		Output::err(msg);
	}


	if( basis_set[2] ) {
		sprintf( msg, "SYSTEM: basis[2] successfully read from pqr {%.5lf %.5lf %.5lf}\n", pbc.basis[2][0], pbc.basis[2][1], pbc.basis[2][2]);
		Output::out(msg);
	} else {
		sprintf(msg,"SYSTEM: unable to read basis[2] from pqr file.\n");
		Output::err(msg);
	}

}




void System::update_pbc() {
/*
	if(  (pbc.cutoff != 0.0)  &&  (checkpoint->movetype != MOVETYPE_VOLUME)  &&  (ensemble != ENSEMBLE_REPLAY)  )
		// If the cutoff is unset, the movetype is MOVETYPE_VOLUME or it is a replay,
		// then this computation needs to happen. Otherwise 
		return; // nothing to do.
*/

	// compute the unit cell volume, cutoff and reciprocal space lattice vectors
	pbc.compute_volume;
	if (pbc.cutoff == 0.) pbc.compute_cutoff;
	pbc.compute_reciprocal;

	// calculate ewald_alpha and polar_ewald_alpha unless manually set
	if (ewald_alpha_set != 1)
		ewald_alpha = 3.5 / pbc.cutoff;
	if (polar_ewald_alpha_set != 1)
		polar_ewald_alpha = 3.5 / pbc.cutoff;
	
}




void System::rebuild_arrays () {
// makes an array to access system's linked list of atoms (and corresponding molecule)

	Molecule * molecule_ptr;
	Atom * atom_ptr;

	free( atom_array );
	free( molecule_array );

	//allocate the arrays
	SafeOps::calloc( molecule_array, natoms, sizeof( Molecule *), __LINE__, __FILE__ );
	SafeOps::calloc( atom_array,     natoms, sizeof( Atom *    ), __LINE__, __FILE__ );

	int n=0;
	//build the arrays
	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			molecule_array[n] = molecule_ptr;
			atom_array[n] = atom_ptr;
			n++;
		}
	}
	natoms = n;
}




void System::countN() {
// count the number of molecules currently in the system excluding frozen, adiabatic, etc.

	Molecule * molecule_ptr;

	observables->N = 0;
	observables->spin_ratio = 0;
	for(molecule_ptr = molecules, observables->N = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		if(!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target)) {
			// update the molecule counter
			observables->N += 1.0;
			// update the nuclear spin ratio
			if(molecule_ptr->nuclear_spin == NUCLEAR_SPIN_PARA)
				observables->spin_ratio += 1.0;
		}

		if(ensemble == ENSEMBLE_NVE)
			N = observables->N;
	}

	return;
}




int System::countNatoms() {
// Counts the number of atoms in the system

	Molecule * m;
	Atom * a;
	int NAtoms = 0;

	for ( m = molecules;  m;  m = m->next )
		for ( a = m->atoms;  a;  a = a->next )
			NAtoms++;

	return NAtoms;
}




void System::pairs() {
// Update everything necessary to describe the complete pairwise system 

	int n;
	Pair *pair_ptr;
	
	// needed for GS ranking metric
	double rmin;

	// get array of atom ptrs
	rebuild_arrays();
	n = natoms;

	// loop over all atoms and pair
	for( int i = 0; i < (n - 1); i++) {
		pair_ptr = atom_array[i]->pairs;
		for( int j = (i + 1); j < n; j++ ) {

			// set the link 
			pair_ptr->atom     = atom_array[j];
			pair_ptr->molecule = molecule_array[j];

			// this:
			//
			//     if ( !pair_ptr->frozen && !(pair_ptr->rd_excluded && pair_ptr->es_excluded) ) { pair_exculsions(...); }
			//
			// is dangerous and has already been responsible for numerous bugs, most recently in UVT runs.
			// After an insert/remove move there is no guarantee that pair_ptr->rd_excluded is properly set.
			
			pair_exclusions( molecule_array[i], molecule_array[j], atom_array[i], atom_array[j], pair_ptr );

			// recalc min image (the induced-induced interaction is needed for frozen atoms)
			if( !pair_ptr->frozen || polarization )
				minimum_image( atom_array[i], atom_array[j], pair_ptr );

			pair_ptr = pair_ptr->next;

		} // for j
	} // for i


	// update the Center-of-Masa of each molecule
	update_com();

	// store wrapped coords
	wrap_all();

	// rank metric 
	if( polar_iterative && polar_gs_ranked ) {
		// determine smallest polarizable separation
		rmin = MAXVALUE;
		for( int i = 0; i < n; i++ ) {
			if ( atom_array[i]->polarizability == 0.0 )	continue;
			for ( pair_ptr = atom_array[i]->pairs; pair_ptr; pair_ptr=pair_ptr->next ) {
				if ( pair_ptr->atom->polarizability == 0.0 ) continue; 
				if ( pair_ptr->rimg < rmin ) rmin = pair_ptr->rimg;
			}
		}

		//calculate rank shits
		for( int i = 0; i < n; i++ )
			atom_array[i]->rank_metric = 0;	

		for( int i = 0; i < n; i++) {
			if ( atom_array[i]->polarizability == 0.0 )
				continue;
			for ( pair_ptr = atom_array[i]->pairs; pair_ptr; pair_ptr=pair_ptr->next ) {
				if ( pair_ptr->atom->polarizability == 0.0 )
					continue;
				if ( pair_ptr->r <= rmin*1.5 ) {
					atom_array[i]->rank_metric += 1.0;
					pair_ptr->atom->rank_metric += 1.0;
				}
			} // for pair_ptr
		} // for i

	}
}




void System::pair_exclusions( Molecule *molecule_i, Molecule *molecule_j, Atom *atom_i, Atom *atom_j, Pair *pair_ptr) {
// set the exclusions and LJ mixing for relevant pairs 

	double si3, sj3, si6, sj6;
	double repul1, repul2, repulmix;

	// recalculate exclusions
	if(   (molecule_i == molecule_j)  &&  !gwp   )   // if both on same molecule, exclude all interactions
	{
		pair_ptr->rd_excluded = 1;
		pair_ptr->es_excluded = 1;

	} else {

		// exclude null repulsion/dispersion interactions 
		if( ((atom_i->epsilon == 0.0) || (atom_i->sigma == 0.0) || (atom_j->epsilon == 0.0) || (atom_j->sigma == 0.0))    
		    &&
		    ( atom_i->c6  == 0.0 &&  atom_i->c8  == 0.0 && atom_i->c10 == 0.0 &&  atom_j->c6  == 0.0 &&  atom_j->c8  == 0.0 &&  atom_j->c10 == 0.0 )
		)
			pair_ptr->rd_excluded = 1;
		else
			pair_ptr->rd_excluded = 0;

		// exclude null electrostatic interactions
		if(  (atom_i->charge == 0.0)  ||  (atom_j->charge == 0.0)  )
			pair_ptr->es_excluded = 1;
		else
			pair_ptr->es_excluded = 0;

	}

	// determine the frozen interactions
	pair_ptr->frozen = atom_i->frozen && atom_j->frozen;

	// get the mixed LJ parameters
	if( !use_sg ) {

		if( waldmanhagler && !cdvdw_sig_repulsion ) { //wh mixing rule
			si3  = atom_i->sigma;
			si3 *= si3*si3;
			si6  = si3*si3;

			sj3  = atom_j->sigma;
			sj3 *= sj3*sj3;
			sj6  = sj3*sj3;

			if(  (atom_i->sigma < 0.0)  ||  (atom_j->sigma < 0.0)  ) {
				pair_ptr->attractive_only = 1;
				pair_ptr->sigma = pow(0.5*(si6+sj6),1./6.);
			} else if ((atom_i->sigma == 0 || atom_j->sigma == 0 )) {
				pair_ptr->sigma = 0;
				pair_ptr->epsilon = sqrt(atom_i->epsilon*atom_j->epsilon); //can't use sigma weights -> div by 0
			} else {
				pair_ptr->sigma = pow(0.5*(si6+sj6),1./6.);
				pair_ptr->epsilon = sqrt(atom_i->epsilon*atom_j->epsilon) * 2.0*si3*sj3/(si6+sj6);
			}
		} 

		else if( halgren_mixing ) { //halgren mixing rules

			if(   atom_i->sigma > 0.0   &&   atom_j->sigma > 0.0   ) {
				pair_ptr->sigma = (atom_i->sigma*atom_i->sigma*atom_i->sigma+atom_j->sigma*atom_j->sigma*atom_j->sigma)/(atom_i->sigma*atom_i->sigma+atom_j->sigma*atom_j->sigma);
			} else {
				pair_ptr->sigma = 0;
			}

			if (atom_i->epsilon>0.0&&atom_j->epsilon>0.0) 	{
				pair_ptr->epsilon = 4*atom_i->epsilon*atom_j->epsilon/pow(sqrt(atom_i->epsilon)+sqrt(atom_j->epsilon),2);
			} else {
				pair_ptr->epsilon = 0;
			}
		}

		else if( cdvdw_9th_repulsion ) {  //9th power mixing for repulsion
			si3      = atom_i->sigma;
			si3     *= si3 * si3;
			si6      = si3 * si3;
			sj3      = atom_j->sigma;
			sj3     *= sj3 * sj3;
			sj6      = sj3 * sj3;
			repul1   = 4.0 * si6 * si6 * atom_i->epsilon;
			repul2   = 4.0 * sj6 * sj6 * atom_j->epsilon;
			repulmix = pow(0.5*(pow(repul1,1./9.) + pow(repul2,1./9.)),9);
			pair_ptr->sigma   = 1.0;
			pair_ptr->epsilon = repulmix/4.0;
		}

		else if ( cdvdw_sig_repulsion ) {  //sigma repulsion for coupled-dipole vdw
			si3  = atom_i->sigma;
			si3 *= si3 * si3; 
			si6  = si3 * si3;
			sj3  = atom_j->sigma;
			sj3 *= sj3 * sj3;
			sj6  = sj3 * sj3;
			pair_ptr->sigma  = pow(0.5*(si6+sj6),1./6.);
			pair_ptr->sigrep = 1.5 * hBar/kB * au2invseconds * atom_i->omega * atom_j->omega * atom_i->polarizability 
			                   * atom_j->polarizability  /  (atom_i->omega + atom_j->omega)  /  pow(pair_ptr->sigma,6);
		}
		else if(  polarvdw  &&  cdvdw_exp_repulsion  ) {  // mix for buckingham repulsion
			// sigma == C, epsilon == rho
			// U = C exp(-R/(2*rho))
			pair_ptr->sigma = pow(pow(atom_i->sigma,atom_i->epsilon) * pow(atom_j->sigma,atom_j->epsilon),1.0/((atom_i->epsilon+atom_j->epsilon)));
			pair_ptr->epsilon = 0.5*(atom_i->epsilon + atom_j->epsilon);
		}
		else if( using_disp_expansion ) {  // mix for buckingham repulsion
			// sigma == r, epsilon == alpha, C ~= 316 K
			// U = C exp(-alpha(R-r))
			
			// forumlas for these (except JR Schmidt's mixing rule) is here http://pubs.acs.org/doi/pdf/10.1021/acs.jpca.6b10295
			if (schmidt_mixing)
			{
				pair_ptr->sigma = 0.5 * (atom_i->sigma + atom_j->sigma);
				pair_ptr->epsilon = (atom_i->epsilon + atom_j->epsilon) * atom_i->epsilon * atom_j->epsilon / (atom_i->epsilon * atom_i->epsilon + atom_j->epsilon * atom_j->epsilon);
			}

			else if (force_mixing)
			{
				double c = 315.7750382111558307123944638;
				double Aii = c * exp(atom_i->epsilon * atom_i->sigma);
				double Ajj = c * exp(atom_j->epsilon * atom_j->sigma);
				double Bii = atom_i->epsilon;
				double Bjj = atom_j->epsilon;
				pair_ptr->epsilon = 2.0 * atom_i->epsilon * atom_j->epsilon / (atom_i->epsilon + atom_j->epsilon);
				double Bij = pair_ptr->epsilon;
				double Aij = pow(pow(Aii * Bii, 1.0 / Bii) * pow(Ajj * Bjj, 1.0 / Bjj), 0.5 / Bij) / Bij;
				pair_ptr->sigma = log(Aij / c) / Bij;
			}

			else if (bohm_ahlrichs_mixing)
			{
				double c = 315.7750382111558307123944638;
				double Aii = c * exp(atom_i->epsilon * atom_i->sigma);
				double Ajj = c * exp(atom_j->epsilon * atom_j->sigma);
				double Bii = atom_i->epsilon;
				double Bjj = atom_j->epsilon;
				pair_ptr->epsilon = 2.0 * atom_i->epsilon * atom_j->epsilon / (atom_i->epsilon + atom_j->epsilon);
				double Bij = pair_ptr->epsilon;
				double Aij = pow(pow(Aii, 1.0 / Bii) * pow(Ajj, 1.0 / Bjj), 0.5 * Bij);
				pair_ptr->sigma = log(Aij / c) / Bij;
			}

			else if (wilson_popelier_mixing)
			{
				double c = 315.7750382111558307123944638;
				double Aii = c * exp(atom_i->epsilon * atom_i->sigma);
				double Ajj = c * exp(atom_j->epsilon * atom_j->sigma);
				double Bii = atom_i->epsilon;
				double Bjj = atom_j->epsilon;
				pair_ptr->epsilon = sqrt(0.5 * (Bii * Bii + Bjj * Bjj));
				double Aij = pow(0.5 * (pow(Aii, 0.4) + pow(Ajj, 0.4)), 1.0 / 0.4);
				pair_ptr->sigma = log(Aij / c) / pair_ptr->epsilon;
			}

			else
			{
				pair_ptr->sigma = 0.5 * (atom_i->sigma + atom_j->sigma);
				pair_ptr->epsilon = 2.0 * atom_i->epsilon * atom_j->epsilon / (atom_i->epsilon + atom_j->epsilon);
			}
			// get mixed dispersion coefficients 
			pair_ptr->c6 = sqrt(atom_i->c6*atom_j->c6)*0.021958709/(3.166811429*0.000001); // Convert H*Bohr^6 to K*Angstrom^6, etc
			pair_ptr->c8 = sqrt(atom_i->c8*atom_j->c8)*0.0061490647/(3.166811429*0.000001); // Dispersion coeffs should be inputed in a.u.

			if(   extrapolate_disp_coeffs   &&   pair_ptr->c6 != 0.0   &&   pair_ptr->c8 != 0.0   )
				pair_ptr->c10  =  49.0 / 40.0 * pair_ptr->c8 * pair_ptr->c8 / pair_ptr->c6;
			else if( extrapolate_disp_coeffs )
				pair_ptr->c10 = 0.0; //either c6 or c8 is zero so lets set c10 to zero too
			else
				pair_ptr->c10 = sqrt(atom_i->c10 * atom_j->c10)  *  0.0017219135/(3.166811429*0.000001);
		}
		else if( c6_mixing ) {
			pair_ptr->sigma = 0.5*(atom_i->sigma + atom_j->sigma);
			if (pair_ptr->sigma!=0.0)
				pair_ptr->epsilon = 64.0  *  sqrt(atom_i->epsilon * atom_j->epsilon)  *  pow(atom_i->sigma,3.0) * pow(atom_j->sigma,3.0) / pow(atom_i->sigma+atom_j->sigma,6.0);
			else
				pair_ptr->epsilon = 0.0;
		}
		else { // lorentz-berthelot
			if(  (atom_i->sigma < 0.0)  ||  (atom_j->sigma < 0.0)  ) {
				pair_ptr->attractive_only = 1;
				pair_ptr->sigma = 0.5 * (fabs(atom_i->sigma) + fabs(atom_j->sigma));
			} else if ((atom_i->sigma == 0 || atom_j->sigma == 0 )) {
				pair_ptr->sigma = 0;
				pair_ptr->epsilon = sqrt(atom_i->epsilon*atom_j->epsilon);
			} else {
				pair_ptr->sigma = 0.5*(atom_i->sigma + atom_j->sigma);
				pair_ptr->epsilon = sqrt(atom_i->epsilon*atom_j->epsilon);
			}
		}
	} // !use_sg

	// ensure that no ES calc for S-S pairs, and ONLY ES for S-everythingelse
	if( spectre ) {

		if( atom_i->spectre && atom_j->spectre ) {  // case for both spectre 

			pair_ptr->rd_excluded = 0;
			pair_ptr->es_excluded = 1;

		} else if(atom_i->spectre || atom_j->spectre) { // case for one spectre 

			pair_ptr->rd_excluded = 1;
			pair_ptr->es_excluded = 0;

		}
	}

	return;
}




void System::minimum_image( Atom *atom_i, Atom *atom_j, Pair *pair_ptr ) {
// perform the modulo minimum image for displacements

	double img[3];
	double d  [3], r,  r2;
	double di [3], ri, ri2;


	// get the real displacement
	pair_ptr->recalculate_energy = 0;	// reset the recalculate flag
	for( int p = 0; p < 3; p++) {
		d[p] = atom_i->pos[p] - atom_j->pos[p];

		// this pair changed and so it will have it's energy recalculated 
		if( d[p] != pair_ptr->d_prev[p]) {
			pair_ptr->recalculate_energy = 1;
			pair_ptr->d_prev[p] = d[p]; // reset w new displacement
		}
	}

	//relative position didn't change. nothing to do here.
	if( pair_ptr->recalculate_energy == 0 )
		return;

	// Matrix multiply:
	// Project the displacement vector into the reciprocal basis and round 
	for( int p = 0; p < 3; p++ ) {
		img[p] = 0;
		for( int q = 0; q < 3; q++ ) {
			img[p] += pbc.reciprocal_basis[q][p] * d[q];
		}
		
		img[p] = rint(img[p]);
	}
	
	// Matrix multiply to project back into our basis 
	for( int p = 0; p < 3; p++) {
		di[p] = 0;
		for( int q = 0; q < 3; q++)
			di[p] += pbc.basis[q][p] * img[q];
	}

	// correct the displacement
	for( int p = 0; p < 3; p++ )
		di[p] = d[p] - di[p];


	// pythagorean terms 
	r2  = 0;
	ri2 = 0;
	for( int p = 0; p < 3; p++) {
		r2 += d[p]*d[p];
		ri2 += di[p]*di[p];
	}
	// finish distance computation for real distance and min image distance
	r  = sqrt(r2 );
	ri = sqrt(ri2);


	// store the results for this pair
	pair_ptr->r = r;

	
	if( isnan(ri) ) {
		// image distance result is bad, use actual distance
		pair_ptr->rimg = r;
		for( int p=0; p<3; p++ )
			pair_ptr->dimg[p] = d[p];
	}
	else {
		// image distance looks good--store it in the node associated with this pair
		pair_ptr->rimg = ri;
		for(int p = 0; p < 3; p++ )
			pair_ptr->dimg[p] = di[p];
	}

	return;
}




void System::flag_all_pairs() {
// flag all pairs to have their energy calculated. 
// flag_all_pairs() needs to be called at simulation start, or can 
// be called periodically to keep the total energy from drifting. 
	
	Molecule *molecule_ptr;
	Atom *atom_ptr;
	Pair *pair_ptr;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
				pair_ptr->recalculate_energy = 1;
}




void System::spectre_wrapall() {
// ensure that the SPECTRE charges are all pulled within the restricted domain 

	Molecule * molecule_ptr;
	Atom * atom_ptr;
	
	double target[3] = {0,0,0};
	double d[3], l;

	// boxlength 
	l = 2.0 * spectre_max_target;

	// get the coordinates of the target particle to wrap around 
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			if( atom_ptr->target )
				for(int p = 0; p < 3; p++)
					target[p] = atom_ptr->pos[p];

		}
	}

	// wrap SPECTRE charges within the box 
	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {

			if( atom_ptr->spectre ) {

				for(int p = 0; p < 3; p++) {
					d[p] = atom_ptr->pos[p] - target[p];
					atom_ptr->pos[p] -= l*rint(d[p]/l);
				}

			}

		} // for atom 
	} // for molecule 

	return;
}




void System::update_com() {
// Computes and updates molecular center of mass for each molecule in the system.
// (spectre and "target" molecules excluded)

	Molecule * molecule_ptr = molecules;
	Atom * atom_ptr;

	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {

		for(int i = 0; i < 3; i++)
			molecule_ptr->com[i] = 0;

		if(   !( molecule_ptr->spectre  ||  molecule_ptr->target )   ) {

			molecule_ptr->mass = 0;
			for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {

				molecule_ptr->mass += atom_ptr->mass;

				for( int i=0; i<3; i++ )
					molecule_ptr->com[i] += atom_ptr->mass * atom_ptr->pos[i];
			}

			for( int i=0; i<3; i++ )
				molecule_ptr->com[i] /= molecule_ptr->mass;
		}
	}
}




int System::wrap_all() {
// Stores PBC-wrapped coords for each (non-Frozen) molecule and
// adjusts the wrappped position of each atom found therein.

	Molecule * molecule_ptr;
	Atom * atom_ptr;
	double d[3], dimg[3];

	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {

		if(  ! molecule_ptr->frozen  ) {
			// get the minimum imaging distance for the com 
			for( int i = 0; i < 3; i++) {
				d[i] = 0;
				for( int j = 0; j < 3; j++) {
					d[i] += pbc.reciprocal_basis[j][i] * molecule_ptr->com[j];
				}
				d[i] = rint(d[i]);
			}

			for( int i = 0; i < 3; i++) {
				dimg[i] = 0;
				for( int j = 0; j < 3; j++)
					dimg[i] += pbc.basis[j][i]*d[j];

				// store the wrapped com coordinate 
				molecule_ptr->wrapped_com[i] = dimg[i];
			}

			// apply the distance to all of the atoms of this molecule 
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
				for( int i = 0; i < 3; i++)
					atom_ptr->wrapped_pos[i] = atom_ptr->pos[i] - dimg[i];
			}

		} else {

			// don't wrap frozen 
			for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
				for( int i = 0; i < 3; i++)
					atom_ptr->wrapped_pos[i] = atom_ptr->pos[i];
			}
		}
	} // molecule 

	return 0;
}




void System::thole_resize_matrices() {
// For uvt runs, resize the A (and B) matrices 

	int i;
	int N_Atoms;
	int dN;
	int oldN;

	// Determine how the number of atoms has changed and realloc matrices 
	oldN = 3*checkpoint->thole_N_atom; //will be set to zero if first time called
	checkpoint->thole_N_atom = countNatoms();
	N_Atoms = 3*checkpoint->thole_N_atom;
	dN = N_Atoms - oldN;

	if( !dN ) return;

	// Grow A matricies by free/malloc (to prevent fragmentation) free the A matrix
	for (i=0; i < oldN; i++) 
		free( A_matrix[i] );
	free( A_matrix );

	// If not iterative, free the B matrix
	if( ! polar_iterative ) {
		for (i=0; i < oldN; i++)
			free( B_matrix[i] );
		free( B_matrix );
	}

	// (RE)allocate the A matrix
	SafeOps::calloc( A_matrix, N_Atoms, sizeof(double*), __LINE__, __FILE__ );
	
	for (i=0; i< N_Atoms; i++ ) 
		SafeOps::malloc( A_matrix[i], N_Atoms * sizeof(double), __LINE__, __FILE__ );

	// (RE)allocate the B matrix if not iterative
	if( ! polar_iterative ) {
		SafeOps::calloc( B_matrix, N_Atoms, sizeof(double*), __LINE__, __FILE__ );
		
		for (i=0; i< N_Atoms; i++ ) 
			SafeOps::malloc(B_matrix[i], N_Atoms * sizeof(double), __LINE__, __FILE__ );
	}

	return;
}




double System::get_rand() {
	
	return dist( mt_rand );

}




int System::calculate_bonds()
{
	int inner_index=0,outer_index=0;
	int bonds=0;
	// double bondlength; (unused variable)
	Molecule * mol;
	Atom * atom;
	Atom * atom2;

	for( mol = molecules; mol; mol=mol->next ) {
		if( mol->frozen ) {
			for( atom=mol->atoms; atom; atom=atom->next, inner_index++ ) {
				if( atom->next ) {
					for( atom2 = atom->next, outer_index=inner_index+1; atom2; atom2=atom2->next, outer_index++ ){
						if( bondlength_check( atom, atom2 ) ) {
							bonds++;
						}
					}
				}
			}
		}
	}
	return bonds;
}




int System::bondlength_check( Atom * atom1, Atom * atom2 )
{
	double distance;
	double gm_mass;      // geometric mean of mass
	int is_bonded;
	double slope=0.0234; // tune to meet bond length expectations
	double yint=0.603;   // tune to meet bond length expectations

	gm_mass=sqrt(atom1->mass * atom2->mass);
	distance=sqrt( pow(atom1->pos[0]-atom2->pos[0],2) + pow(atom1->pos[1] - atom2->pos[1],2) + pow(atom1->pos[2] - atom2->pos[2],2));

	if(  distance  <  (gm_mass * slope + yint) * max_bondlength ) 
		is_bonded=1;
	else
		is_bonded=0;

	return is_bonded;
}




void System::calc_system_mass() {
// calculate and set observable variables for the frozen and total mass of the system

	Molecule * molecule_ptr;
	
	observables->total_mass = 0;
	observables->frozen_mass = 0;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		observables->total_mass += molecule_ptr->mass;
		if( molecule_ptr->frozen || molecule_ptr->adiabatic )
			observables->frozen_mass += molecule_ptr->mass;
	}

	return;
}




void System::count_sorbates() {
	int i;
	Molecule		* molecule_ptr;

	// Zero every count (N) for each sorbate in the averages list
	for ( i=0; i<sorbateCount; i++ ) 
		sorbateInfo[i].currN = 0;

	// Count each sorbate in the system and record the total in the corresponding entry in the sorbate averages list.
	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( i=0; i<sorbateCount; i++ ) {
			if( SafeOps::iequals( sorbateInfo[i].id, molecule_ptr->moleculetype )) {
				sorbateInfo[i].currN++;
				break;
	}}}
}




void System::allocateStatisticsMem() {

	SafeOps::calloc( nodestats,               1, sizeof(nodestats_t      ), __LINE__, __FILE__ );
	SafeOps::calloc( avg_nodestats,           1, sizeof(avg_nodestats_t  ), __LINE__, __FILE__ );
	SafeOps::calloc( observables,             1, sizeof(observables_t    ), __LINE__, __FILE__ );
	SafeOps::calloc( avg_observables,         1, sizeof(avg_observables_t), __LINE__, __FILE__ );
	SafeOps::calloc( checkpoint,              1, sizeof(checkpoint_t     ), __LINE__, __FILE__ );
	SafeOps::calloc( checkpoint->observables, 1, sizeof(observables_t    ), __LINE__, __FILE__ );
	SafeOps::calloc( grids,                   1, sizeof(grid_t           ), __LINE__, __FILE__ );
	SafeOps::calloc( grids->histogram,        1, sizeof(histogram_t      ), __LINE__, __FILE__ );
	SafeOps::calloc( grids->avg_histogram,    1, sizeof(histogram_t      ), __LINE__, __FILE__ );
}



void System::car2basis(double a, double b, double c, double alpha, double beta, double gamma) {
	//converts .car style basis to MPMC style basis if user opts for carbasis.
	//i.e. called when input contains carbasis x x x x x x
	double b0[3] = { 0,0,0 };
	double b1[3] = { 0,0,0 };
	double b2[3] = { 0,0,0 };

	b0[0] = a;
	b0[1] = b * cos(pi / 180.0 * gamma);
	b0[2] = c * cos(pi / 180.0 * beta);

	b1[0] = 0;
	b1[1] = b * sin(pi / 190.0 * gamma);
	b1[2] = ((b * c * cos(pi / 180.0 * alpha)) - (b0[1] * b0[2])) / b1[1];

	b2[0] = 0;
	b2[1] = 0;
	b2[2] = sqrt(c * c - b0[2] * b0[2] - b1[2] * b1[2]);

	//Transposing manually
	pbc.basis[0][0] = b0[0];
	pbc.basis[0][1] = b1[0];
	pbc.basis[0][2] = b2[0];

	pbc.basis[1][0] = b0[1];
	pbc.basis[1][1] = b1[1];
	pbc.basis[1][2] = b2[1];
	
	pbc.basis[2][0] = b0[2];
	pbc.basis[2][1] = b1[2];
	pbc.basis[2][2] = b2[2];

}




int System::printAtoms() {

	char       linebuf[maxLine] = { 0 };
	Molecule * m                = nullptr;
	Atom     * a                = nullptr;
	int        molCount         = 0;
	int        atomNum          = 0;

	for( m=molecules; m; m=m->next ) {
		molCount++;
		atomNum = 0;
		sprintf(linebuf, "\nMOLECULE %d\n", molCount);
		Output::out(linebuf);
		for( a=m->atoms; a; a=a->next) {
			sprintf(linebuf, "M%0dA%0d: (%s)   %lf, %lf, %lf\n", molCount, atomNum, a->atomtype, a->pos[0], a->pos[1], a->pos[2]);
			Output::out(linebuf);
		}
		Output::out("\n");
	}
	return 1;
}