#pragma once
#ifndef __SYSTEM_H 
#define __SYSTEM_H
#ifdef _MPI
	#include <mpi.h>
#endif
#include <random>
#include <stdint.h>
#include <vector>

class Atom;
class Pair;

#include "constants.h"
#include "Molecule.h"
#include "PeriodicBoundary.h"






class System
{
public:

    int printAtoms();

	typedef struct _ptemp {
		int     * index; 
		double  * templist;
	} ptemp_t;

	typedef struct avg_observables {
		// these are uncorrelated averages of the observables 
		double energy, energy_sq, energy_error;
		double N, N_sq, N_error;
		double coulombic_energy, coulombic_energy_sq, coulombic_energy_error;
		double rd_energy, rd_energy_sq, rd_energy_error;
		double polarization_energy, polarization_energy_sq, polarization_energy_error;
		double vdw_energy, vdw_energy_sq, vdw_energy_error;
		double three_body_energy, three_body_energy_sq, three_body_energy_error;
		double dipole_rrms, dipole_rrms_sq, dipole_rrms_error;
		double density, density_sq, density_error;
		double pore_density, pore_density_error;
		double percent_wt, percent_wt_error;
		double percent_wt_me, percent_wt_me_error;
		double excess_ratio, excess_ratio_error;

		// needed for heat capacity error propagation 
		double energy_sq_sq, energy_sq_error;

		// for NVE MC 
		double kinetic_energy, kinetic_energy_sq, kinetic_energy_error;
		double temperature, temperature_sq, temperature_error;

		// for NPT 
		double volume, volume_sq, volume_error;

		// needed for qst 
		double NU;

		// ortho:para spin ratio 
		double spin_ratio, spin_ratio_sq, spin_ratio_error;

		// these quantities are node stats, not observables 
		double boltzmann_factor, boltzmann_factor_sq, boltzmann_factor_error;
		double cavity_bias_probability, cavity_bias_probability_sq, cavity_bias_probability_error;
		double polarization_iterations, polarization_iterations_sq, polarization_iterations_error;

		double acceptance_rate;
		double acceptance_rate_insert, acceptance_rate_remove, acceptance_rate_displace;
		double acceptance_rate_adiabatic, acceptance_rate_spinflip, acceptance_rate_volume, acceptance_rate_ptemp, acceptance_rate_beadPerturb;
		

		// these quantities are based on averages; the error is not easily calculated 
		double qst;
		double heat_capacity;
		double heat_capacity_error;
		double compressibility;
		double compressibility_error;
	} avg_observables_t;

	typedef struct _observables {
		double energy,
		       coulombic_energy,
		       rd_energy,
		       polarization_energy,
		       vdw_energy,
		       three_body_energy,
		       dipole_rrms,
		       kinetic_energy,     // for NVE 
		       temperature,        // for NVE 
		       volume,             // for NPT
		       N,
		       NU,
		       spin_ratio,         // ortho:para spin ratio 
		       frozen_mass,
		       total_mass;         //updated in average.c
	} observables_t;

	typedef struct _checkpoint {
		int              movetype, 
		                 biased_move,
		                 thole_N_atom; //used for keeping track of thole matrix size (allocated)
		Molecule       * molecule_backup, 
		               * molecule_altered,
		               * head,
		               * tail;
		observables_t  * observables;
	} checkpoint_t;

	typedef struct _cavity {
		int occupancy;
		double pos[3];
	} cavity_t;
	
	typedef struct _avg_node_stats {
		int    counter;
		double boltzmann_factor, 
		       boltzmann_factor_sq,
		       acceptance_rate,
		       acceptance_rate_insert,
		       acceptance_rate_remove,
		       acceptance_rate_displace,
		       acceptance_rate_ptemp,
		       acceptance_rate_adiabatic,
		       acceptance_rate_spinflip,
		       acceptance_rate_volume,
		       //acceptance_rate_premp,
			   acceptance_rate_beadPerturb,
		       cavity_bias_probability,
		       cavity_bias_probability_sq,
		       polarization_iterations,
		       polarization_iterations_sq;
	} avg_nodestats_t;
	
	typedef struct _nodestats {
		int    accept, 
		       reject;

		int    accept_insert,
			   accept_remove,
			   accept_displace,
			   accept_adiabatic,
			   accept_spinflip,
			   accept_volume,
			   accept_ptemp,
			   accept_beadPerturb,

			   reject_insert,
			   reject_remove,
			   reject_displace,
			   reject_adiabatic,
			   reject_spinflip,
			   reject_volume,
			   reject_ptemp,
			   reject_beadPerturb;

		double boltzmann_factor,
			   acceptance_rate,
			   acceptance_rate_insert,
			   acceptance_rate_remove,
			   acceptance_rate_displace,
			   acceptance_rate_adiabatic,
			   acceptance_rate_spinflip,
			   acceptance_rate_volume,
			   acceptance_rate_ptemp,
			   acceptance_rate_beadPerturb,
		       cavity_bias_probability,
		       polarization_iterations;
	} nodestats_t;

	typedef struct _histogram {
		int ***grid;
		int x_dim, y_dim, z_dim;
		double origin[3];
		double delta[3][3];
		int count[3];
		int n_data_points;
		int norm_total;
	} histogram_t;
	
	typedef struct _grid {
		histogram_t * histogram;
		histogram_t * avg_histogram;
	} grid_t;

	// used for storing global sorbate averages
	typedef struct _sorbateAverages {
		double avgN, avgN_sq, avgN_err;                            // average sorbate count
		double percent_wt, percent_wt_sq, percent_wt_err;          // weight percent for this sorbate (sorb_mass / total_mass)
		double percent_wt_me, percent_wt_me_sq, percent_wt_me_err; // weight percent for this sorbate (sorb_mass / frozen_mass)
		double excess_ratio, excess_ratio_sq, excess_ratio_err;    // excess adsorption ratio
		double pore_density, pore_density_sq, pore_density_err;    // mass / volume of pore
		double density, density_sq, density_err;                   // mass of sorbate / volume
		double selectivity, selectivity_err;                       // sorbate's selectivity ratio relative to all other sorbates in the insert list.
		struct _sorbateAverages *next;
	} sorbateAverages_t; 

	// local sorbate data array
	typedef struct _sorbateInfo {
		char   id[16];        // identifying tag for the sorbate, e.g. CH4, CO2 or H2
		double mass;          // mass of this sorbate.
		int    currN;         // sorbate count for the current step
		double percent_wt;
		double percent_wt_me;
		double excess_ratio;
		double pore_density;
		double density;
	} sorbateInfo_t;

	//stores vdw energies for each molecule within the coupled dipole model
	typedef struct _vdw {
		char mtype[maxLine];
		double energy;
		struct _vdw * next;
	} vdw_t;

	typedef struct _mtx {
		int dim;
		double * val;
	} mtx_t;

	// For accomodating an arbitrary number of fit_input files
	typedef union {
		int count;
		char *filename;
	} intechar;
	typedef struct _fileNode {
		intechar data;            // data.count is for use in the list head only & indicates the # of list elements
        struct _fileNode *next;   // data.filename is for use in the nodes
	} fileNode_t;

	typedef struct mpiData_t {
		int                msgsize;
		observables_t    * observables;
		avg_nodestats_t  * avg_nodestats;
		sorbateInfo_t    * sinfo;
		double           * temperature;
		char             * snd_strct; 
		char             * rcv_strct;
	} mpiData;

	
	// System.cpp
	~System();
	System();
	System( const System &sd );
	#ifdef _MPI
		MPI_Datatype msgtype;
	#endif		
	bool setup_simulation_box();
	void read_molecules( FILE *fp         );
	void read_molecules( char *input_file );
	void read_pqr_box  ( char *input_file );
	void update_pbc();
	void thole_resize_matrices();
	void rebuild_arrays();
	void countN();
	int  countNatoms();
	void pairs();
	void pair_exclusions( Molecule *molecule_i, Molecule *molecule_j, Atom *atom_i, Atom *atom_j, Pair *pair_ptr);
	void minimum_image( Atom *atom_i, Atom *atom_j, Pair *pair_ptr );
	void flag_all_pairs();
	void spectre_wrapall();
	void update_com();
	int  wrap_all();
	double get_rand();
	int calculate_bonds();
	int bondlength_check( Atom *atom1, Atom *atom2 );
	void calc_system_mass();
	void count_sorbates();
	void allocateStatisticsMem();
	


	// System.Averages.cpp
	void update_root_averages( observables_t *observables );
	void update_sorbate_info();
	void update_root_sorb_averages( sorbateInfo_t * sinfo );
	void clear_avg_nodestats();
	void update_root_nodestats( avg_nodestats_t *avg_nodestats, avg_observables_t *avg_observables );

	
	// System.Cavity.cpp
	void cavity_update_grid();
	void update_cavity_probability();
	void update_cavity_volume();
	bool is_point_empty( double x, double y, double z);
	void setup_cavity_grid();
	double cavity_absolute_check();
	

	// System.Energy.cpp
	double energy();
		
	double * getsqrtKinv( int N );
	double sum_eiso_vdw ( double * sqrtKinv );
	double calc_e_iso( double * sqrtKinv, Molecule * mptr );
	double e2body( Atom * atom, Pair * pair, double r);
	double twobody();

	static mtx_t * build_M( int dim, int offset, double ** Am, double * sqrtKinv );
	static mtx_t * alloc_mtx ( int dim );
	static void free_mtx( mtx_t * M );
	static void printevects( mtx_t * M );
	static double eigen2energy( double * eigvals, int dim ); // , double temperature );
	
	double * lapack_diag( mtx_t * M, int jobtype );

	double vdw();
	double fh_vdw_corr();
	double fh_vdw_corr_2be();
	static void free_vdw_eiso( vdw_t * vdw_eiso_info );
	double lr_vdw_corr();
	
	double anharmonic();
	static double anharmonic_energy(double k, double g, double x);
	static double anharmonic_fk(double temperature, double mass, double k, double g, double x);
	static double anharmonic_fh_second_order(double temperature, double mass, double k, double g, double x);
	static double anharmonic_fh_fourth_order(double temperature, double mass, double k, double g, double x);
	

	// System.Energy.AxilrodTeller.cpp
	double axilrod_teller ();
	

	// System.Energy.Coulombic.cpp
	double coulombic();
	double coulombic_kinetic_gwp();
	static double coulombic_nopbc( Molecule * molecules );
	double coulombic_nopbc_gwp();
	double coulombic_real();
	double coulombic_real_FH( Molecule * molecule_ptr, Pair *pair_ptr, double gaussian_term, double erfc_term );
	double coulombic_reciprocal();
	double coulombic_self();
	double coulombic_wolf();
	
	
	//System.Energy.DispExp.cpp
	double disp_expansion();
	static double tt_damping(int n, double br);
	double disp_expansion_lrc_self( Atom * atom_ptr, const double cutoff );
	double disp_expansion_lrc( Pair * pair_ptr, const double cutoff );
	double exp_fh_corr( Molecule * molecule_ptr, Pair * pair_ptr, int order, double pot );
	double exp_crystal_self( Atom * aptr, double cutoff );
	double exp_lrc_self( Atom * atom_ptr, double cutoff );
	

	//System.Energy.Dreiding.cpp
	double dreiding();
	static double dreiding_nopbc( Molecule *molecules ); 
	

	// System.Energy.ExpRepulsion.cpp
	double exp_repulsion();
	double exp_lrc_corr( Atom * atom_ptr,  Pair * pair_ptr, double cutoff );
	

	// System.Energy.LJ.cpp
	double lj();
	double lj_lrc_corr( Atom * atom_ptr,  Pair * pair_ptr, double cutoff );
	double lj_fh_corr( Molecule * molecule_ptr, Pair * pair_ptr, int order, double term12, double term6 );
	double rd_crystal_self( Atom * aptr, double cutoff );
	double lj_lrc_self( Atom * atom_ptr, double cutoff );
	double lj_buffered_14_7();
	double lj_buffered_14_7_nopbc();
	

	// System.Energy.SG.cpp
	double sg();
	static double sg_nopbc( Molecule *molecules );
		
	
	// System.Energy.Polar.cpp
	double   polar();
	void     thole_amatrix();
	void     zero_out_amatrix ( int N );
	void     ewald_full();
	void     recip_term();
	void     real_term();
	void     init_dipoles_ewald();
	void     clear_ef_induced();
	void     induced_recip_term();
	void     induced_real_term();
	void     induced_corr_term();
	void     calc_dipole_rrms();
	void     new_dipoles(int count);
	int      are_we_done_yet( int iteration_counter );
	void     ewald_palmo_contraction();
	void     thole_field();
	void     thole_field_nopbc();
	void     thole_field_wolf();
	double * polar_wolf_alpha_lookup_init();
	double   polar_wolf_alpha_getval( double r );
	void     ewald_estatic();
	int      thole_iterative();
	void     init_dipoles();
	void     contract_dipoles( int * ranked_array );
	void     palmo_contraction( int * ranked_array );
	void     update_ranking( int * ranked_array );
	double   damp_factor( double t, int i );
	double   get_dipole_rrms();
	void     thole_bmatrix();
	void     thole_bmatrix_dipoles();
	void     thole_polarizability_tensor();
	

	// System.Histogram.cpp
	void setup_histogram();
	void frac2cart(double *answer, double *frac );
	void cart2frac(double *answer, double *cart );
	void setup_dx_variables( histogram_t *hist );
	void allocate_histogram_grid();
	static double magnitude( double *vector );
	void offset_dx_origin( double *real_origin_cartesian, histogram_t *hist );
	void setup_deltas( histogram_t *hist );
	void zero_grid( int ***grid );
	void population_histogram();
	void wrap1coord( double *unwrapped, double *wrapped );
	void compute_bin( double *cart_coords, int *bin_vector );
	void update_root_histogram();
	void write_histogram( FILE *fp_out, int ***grid );
	

	// System.MonteCarlo.cpp
	bool        mc();
	void        do_checkpoint();
	static int  pick_Gibbs_move( std::vector<System*> &sys );
	static void backup_observables(std::vector<System*> &sys );
	void        make_move();
	static void make_move_Gibbs(std::vector<System*> &sys);
	void        enumerate_particles();
	void        displace_1D( Molecule *molecule, double scale );
	void        spectre_displace( Molecule *molecule, double trans_scale, double max_charge, double max_target );
	void        spectre_charge_renormalize();
	void        displace( Molecule *molecule, const PeriodicBoundary &pbc, double trans_scale, double rot_scale );
	void        volume_change();
	static void volume_change_Gibbs( std::vector<System*> &sys );
	void        boltzmann_factor( double initial_energy, double final_energy, double rot_partfunc);
	void        register_accept();
	void        restore();
	void        unupdate_pairs_insert();
	void        unupdate_pairs_remove();
	void        revert_volume_change();
	void        register_reject();
	void        temper_system( double current_energy );
	double      mc_initial_energy();
	mpiData     setup_mpi();
	void        setup_mpi_dataStructs();
	void        setup_mpi_dataStructs( mpiData &md );
	void        do_corrtime_bookkeeping(mpiData &mpi);
	void        output_file_data();


	// System.MPI.cpp
	void mpi_copy_histogram_to_sendbuffer( char *snd, int ***grid );
	void mpi_copy_rcv_histogram_to_data( char *rcv, int ***histogram );


	// System.Output.cpp
	int  open_files();
	void close_files();
	void write_frozen( FILE *fp_frozen );
	int  count_frozen();
	void print_frozen_coords( FILE *fp_frozen );
	void print_frozen_bonds ( FILE *fp_frozen );
	void print_frozen_masses( FILE *fp_frozen );
	void print_frozen_colors( FILE *fp_frozen );
	void write_observables( FILE *fp_energy, observables_t * observables, double core_temp);
	void write_observables_csv( FILE *fp_energy_csv, observables_t * observables, double core_temp);
	int  write_averages();
	int  write_averages(int sysNum);
	int  write_averages(const char *sysID);
	static void track_ar(nodestats_t *ns);
	static void update_nodestats( nodestats_t *nodestats, avg_nodestats_t *avg_nodestats );
	void write_states();
	FILE * open_traj_file();
	int write_molecules_wrapper( char * filename );
	int write_molecules( FILE * fp );
	FILE * open_dipole_file();
	void write_dipole();
	FILE * open_field_file();
	void write_field();
	int write_performance( int i );
	


	// System.Pairs.cpp
	void allocate_pair_lists();
	void update_pairs_insert();
	void update_pairs_remove();

	

	
public:

	// Compilation type flags and data
	int         cuda;
	int         opencl;
	#ifdef OPENCL
		ocl_t * ocl;
	#endif


	int         ensemble;
	char        job_name[maxLine];     // (CRC)
	
	// Monte Carlo Controls
	int         numsteps;              // Total number of MC simulation steps to perform
	int         step;                  // Current MC step
	int         corrtime;              // Number of steps between MC correlation times
	int         ptemp_freq;

	double      move_factor;
	double      rot_factor;
	double      PI_bead_perturb_factor; // option to shrink or expand PI chain perturbations. Default value is 1 (no change)
	double      last_volume;            // NPT option
	double      volume_change_factor;   // NPT option
	double      adiabatic_probability,
	            gwp_probability,
	            insert_probability,
	            spinflip_probability,
	            volume_probability,
		        transfer_probability;

	// Path-Integral settings
	double      bead_perturb_probability;   //  PI option--probability for move that changes bead configuration
	


	// io filenames
	char        dipole_output     [maxLine],
	            energy_output     [maxLine],
	            energy_output_csv [maxLine],
	            field_output      [maxLine],
	            frozen_output     [maxLine],
	            histogram_output  [maxLine],
	            insert_input      [maxLine],
	            pqr_input         [maxLine],
	            pqr_input_B       [maxLine],
	            pqr_output        [maxLine],
	            pqr_restart       [maxLine],
	            surf_output       [maxLine],
	            traj_input        [maxLine],
	            traj_output       [maxLine],
	            virial_output     [maxLine];
	    
	 

	// Observables
	double      temperature,
	            pressure,
	            free_volume,
	            total_energy,
	            N;
	int         fugacitiesCount;
	double      fugacities[maxTokens];
	

	// Inter-node/Inter-system Data relay
	mpiData     mpi_data;


	// Simulation Flags & Data
	int         gwp;
	int         long_output;       // Flag: signals request to print extended (%11.6f) coordinates
	int         parallel_restarts; // Flag: signals that this run is a restart of a parallel job
	double      max_bondlength;    // Bond threshold (re:output files)
	
	

	cavity_t   *** cavity_grid;
	int            cavity_bias, 
	               cavity_grid_size,
	               cavities_open;
	double         cavity_autoreject_repulsion,
	               cavity_autoreject_scale, 
	               cavity_radius, 
	               cavity_volume;
	// Auto-reject Options --- first is in terms of sigma and only
	// applies to LJ; latter is in Angstroms and applies to all pairs
	int            cavity_autoreject;           // Flag: autoreject in terms of sigma--only applies to LJ
	int            cavity_autoreject_absolute;  // Flag: autoreject in Angstroms and applies to all pairs

	// Parallel Tempering Options
	int            parallel_tempering;
	double         max_temperature;
	ptemp_t      * ptemp;
	
	// Info regarding  periodic boundary and unit cell geometry
	int              wrapall;  // Flag: wrap all option requested
	PeriodicBoundary pbc;      // Periodic boundary conditions of THIS system 
	
	
	// (P)RNG
	std::mt19937   mt_rand;
	std::uniform_real_distribution<double> dist{0,1};
	int            preset_seed_on;  //for manual specification of random seeds
	uint32_t       preset_seed;
	
	int            read_pqr_box_on; //read box basis from pqr
	
	// Simulated Annealing
	int            simulated_annealing, 
	               simulated_annealing_linear;
	double         simulated_annealing_schedule,
	               simulated_annealing_target;

	// Spectre
	int            spectre;
	double         spectre_max_charge, 
	               spectre_max_target;


	// Energy-corrections
	int            feynman_hibbs,
	               feynman_kleinert,
	               feynman_hibbs_order;
	int            vdw_fh_2be;          // Flag: 2BE method for polarvdw requested
	int            rd_lrc, 
	               rd_crystal,
	               rd_crystal_order;

	// uVT Fugacity Functions
	int            h2_fugacity, 
	               co2_fugacity, 
	               ch4_fugacity,
	               n2_fugacity,
	               user_fugacities;


	// Force-field Options
	int            rd_only, 
	               rd_anharmonic;
	double         rd_anharmonic_k, 
	               rd_anharmonic_g;
	int            c6_mixing, 
	               damp_dispersion,
	               disp_expansion_mbvdw,
	               use_dreiding, 
	               extrapolate_disp_coeffs, 
	               halgren_mixing,
	         
	               midzuno_kihara_approx,
	               schmidt_ff,
	               use_sg,
	               waldmanhagler;

	bool           using_axilrod_teller, 
	               using_lj_buffered_14_7,
	               using_disp_expansion;

	// ES Options
	int            wolf, 
	               ewald_alpha_set, 
	               ewald_kmax, 
	               polar_ewald_alpha_set;
	double         ewald_alpha,
		           polar_ewald_alpha;
	
	
	// Thole Options
	int            polarization; // Flag signaling that polarization calculation has been requested
	int            polarvdw,
		           polarizability_tensor,
		           cdvdw_exp_repulsion,
		           cdvdw_sig_repulsion,
		           cdvdw_9th_repulsion;
	int            iterator_failed; //flag set when iterative solver fails to converge (when polar_precision is used)
	int            polar_iterative,
		           polar_ewald,
		           polar_ewald_full,
		           polar_zodid,
		           polar_palmo,
		           polar_rrms,
		           polar_gs;

	
	int            polar_gs_ranked;  // Flag indicating if the ranked gauss-seidell algorithm will be used in polar calculations.
	int            polar_sor,
	               polar_esor,
	               polar_max_iter,
	               polar_wolf,
	               polar_wolf_full,
	               polar_wolf_alpha_lookup;
	double         polar_wolf_alpha, 
	               polar_gamma,
	               polar_damp,
	               field_damp,
	               polar_precision,
	             * polar_wolf_alpha_table,
	               polar_wolf_alpha_lookup_cutoff;
	int            polar_wolf_alpha_table_max;     //stores the total size of the array double polar_wolf_alpha_table[]
	int            damp_type;
	double      ** A_matrix;       // A matrix (Thole polarization) 
	double      ** B_matrix;       // B matrix (Thole polarization)
	double         C_matrix[3][3]; // Polarizability tensor 

	vdw_t        * vdw_eiso_info; //keeps track of molecule vdw self energies
	

	//misc
	double         scale_charge;
	int            independent_particle;

	// insertions from a separate linked list
	int            num_insertion_molecules;   // the number of elements found in the Molecule::insertion_molecules and Molecule::insertion_molecules_array
	Molecule     * insertion_molecules;       // linked list of molecules to be randomly inserted
	Molecule    ** insertion_molecules_array; // array providing direct access to the linked list holding molecules from which insertions are chosen


	// quantum rotation stuff
	int            quantum_rotation, 
	               quantum_rotation_hindered,
	               quantum_rotation_l_max,  
	               quantum_rotation_level_max, 
	               quantum_rotation_phi_max,
	               quantum_rotation_sum,
	               quantum_rotation_theta_max,
	               quantum_vibration;
	double         quantum_rotation_B,
	               quantum_rotation_hindered_barrier;
	

	// histogram stuff
	grid_t       * grids;
	int            calc_hist;         // flag to calculate a 3D histogram 
	double         hist_resolution;
	int            n_histogram_bins;

	// atom array	
	int            natoms;
	Atom        ** atom_array;
	Molecule    ** molecule_array;

	//replay option
	int            calc_pressure;
	double         calc_pressure_dv;

	Molecule            * molecules;

	nodestats_t         * nodestats;
	avg_nodestats_t     * avg_nodestats;
	observables_t       * observables;
	avg_observables_t   * avg_observables;

	// Linked list head that will keep track of separate average-observables for each sorbate in the system.
	int                  sorbateCount;   // Number of sorbates in the system.
	sorbateInfo_t      * sorbateInfo;    // stores an array of sorbate Info
	int                  sorbateInsert;  // which sorbate was last inserted
	sorbateAverages_t  * sorbateGlobal;  // where the global average is stored

	checkpoint_t       * checkpoint;

	FILE               * fp_energy;
	FILE               * fp_energy_csv;
	FILE               * fp_field;
	FILE               * fp_histogram;
	FILE               * fp_frozen;
	FILE               * fp_traj_replay;
	FILE               * fp_surf;



	 // surface fitting options
	///////////////////////////////////////////////////////////

	fileNode_t     fit_input_list;
	int            ee_local;                     // Exhaustive Enumeration option
	double         fit_best_square_error,
	               fit_max_energy,
	               fit_schedule,
	               fit_start_temp;
	int            fit_boltzmann_weight;
	int            surf_decomp,
	               surf_fit_arbitrary_configs,
	               surf_qshift_on,
	               surf_weight_constant_on,
	               surf_global_axis_on;
	double         surf_ang,
	               surf_inc,
	               surf_min,
	               surf_max;
	int            surf_descent,
	               surf_preserve,
	               surf_preserve_rotation_on;
	double         surf_preserve_rotation_alpha1,
	               surf_preserve_rotation_alpha2,
	               surf_preserve_rotation_beta1,
	               surf_preserve_rotation_beta2,
	               surf_preserve_rotation_gamma1,
	               surf_preserve_rotation_gamma2;
	double         range_eps, 
	               range_sig, 
	               step_eps,
	               step_sig;
	int            surf_print_level;             // sets the amount of output (1-6) that correspond to the nested loops in surface.c
	int            surf_scale_alpha_on, 
	               surf_scale_epsilon_on,
	               surf_scale_omega_on,
	               surf_scale_sigma_on,
	               surf_scale_pol_on,
	               surf_scale_q_on,
	               surf_scale_r_on,
	               surf_scale_c6_on,
	               surf_scale_c8_on, 
   	               surf_scale_c10_on;
	double         surf_scale_epsilon, 
	               surf_scale_r,
	               surf_scale_omega,
	               surf_scale_sigma,
	               surf_scale_q,
	               surf_scale_alpha,
	               surf_scale_pol,
	               surf_scale_c6,
	               surf_scale_c8,
	               surf_scale_c10;
	int            surf_virial;
	double         surf_quadrupole,
	               surf_weight_constant;

};

#endif // __SYSTEM_H