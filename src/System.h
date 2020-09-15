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

using uint = size_t;


static const double  ewald_alpha_default = 0.5;
static const int     ewald_kmax_default  = 7;
static const int     ptemp_freq_default  = 20; // default frequency for parallel tempering bath swaps
static const double  wolf_alpha_lookup_cutoff_default = 30.0; //angstroms







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
		       kinetic_energy,     // for NVE & PI NVT
		       temperature,        // for NVE 
		       volume,             // for NPT
		       N,
		       NU,
		       spin_ratio,         // ortho:para spin ratio 
		       frozen_mass,
		       total_mass;         //updated in average.c
		inline double potential() { 
			return coulombic_energy + rd_energy + polarization_energy + vdw_energy + three_body_energy;
		}
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
	void recursive_delete_molecules(Molecule * &m);
	System();
	System( const System &sd );
	#ifdef _MPI
		MPI_Datatype msgtype;
	#endif	

	// restore observables from checkpointed backup

	bool setup_simulation_box();
	void read_molecules( FILE *fp         );
	void read_molecules( char *input_file );
	void read_pqr_box  ( char *input_file );
	void update_pbc();
	void thole_resize_matrices();
	void rebuild_arrays();
	unsigned int  countN();
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
	void update_root_nodestats();
	void update_root_nodestats( avg_nodestats_t *avg_nodestats, avg_observables_t *avg_observables );
	void compile_MC_algorithm_stats();

	
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
	
	double axilrod_teller ();
	
	double coulombic();
	double coulombic_kinetic_gwp();
	static double coulombic_nopbc( Molecule * molecules );
	double coulombic_nopbc_gwp();
	double coulombic_real();
	double coulombic_real_FH( Molecule * molecule_ptr, Pair *pair_ptr, double gaussian_term, double erfc_term );
	double coulombic_reciprocal();
	double coulombic_self();
	double coulombic_wolf();
	
	double disp_expansion();
	static double tt_damping(int n, double br);
	double disp_expansion_lrc_self( Atom * atom_ptr, const double cutoff );
	double disp_expansion_lrc( Pair * pair_ptr, const double cutoff );
	double exp_fh_corr( Molecule * molecule_ptr, Pair * pair_ptr, int order, double pot );
	double exp_crystal_self( Atom * aptr, double cutoff );
	double exp_lrc_self( Atom * atom_ptr, double cutoff );
	
	double dreiding();
	static double dreiding_nopbc( Molecule *molecules ); 
	
	double exp_repulsion();
	double exp_lrc_corr( Atom * atom_ptr,  Pair * pair_ptr, double cutoff );
	
	double lj();
	double lj_lrc_corr( Atom * atom_ptr,  Pair * pair_ptr, double cutoff );
	double lj_fh_corr( Molecule * molecule_ptr, Pair * pair_ptr, int order, double term12, double term6 );
	double rd_crystal_self( Atom * aptr, double cutoff );
	double lj_lrc_self( Atom * atom_ptr, double cutoff );
	double lj_buffered_14_7();
	double lj_buffered_14_7_nopbc();
	
	double sg();
	static double sg_nopbc( Molecule *molecules );
	
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
	void        make_move();
	static void make_move_Gibbs(std::vector<System*> &sys);
	void        enumerate_particles();
	void        displace_1D( Molecule *molecule, double scale );
	void        spectre_displace( Molecule *molecule, double trans_scale, double max_charge, double max_target );
	void        spectre_charge_renormalize();
	void        displace( Molecule *molecule, const PeriodicBoundary &pbc, double trans_scale, double rot_scale );
	void        volume_change();
	static void volume_change_Gibbs( std::vector<System*> &sys );
	void        boltzmann_factor( double initial_energy, double final_energy);
	void        restore();
	void        unupdate_pairs_insert();
	void        unupdate_pairs_remove();
	void        revert_volume_change();
	void        register_accept();
	void        register_accept(int movetype);
	void        register_reject();
	void        register_reject(int movetype);
	void        temper_system( double current_energy );
	double      mc_initial_energy();
	mpiData     setup_mpi();
	void        setup_mpi_dataStructs(int qty);
	void        setup_mpi_dataStructs( mpiData &md );
	void        setup_mpi_dataStructs(mpiData& md, size_t qty);
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
	void append_observables_to_output_file(); // "self-documenting" wrapper to write_observables()
	void write_observables( FILE *fp_energy, observables_t * observables, double core_temp);
	void append_observables_to_csv_file(); // "self-documenting" wrapper to write_observables_csv()
	void write_observables_csv( FILE *fp_energy_csv, observables_t * observables, double core_temp);
	int  display_averages();
	int  display_averages(int sysNum);
	int  display_averages(const char *sysID);
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
	int write_performance( unsigned int i );
	


	// System.Pairs.cpp
	void allocate_pair_lists();
	void update_pairs_insert();
	void update_pairs_remove();

	

	
public:
	// a priori system defaults


	// Compilation type flags and data
	int         cuda   = 0;
	int         opencl = 0;
	#ifdef OPENCL
		ocl_t * ocl = nullptr;
	#endif


	int         ensemble = 0;
	char        job_name[maxLine] = { "untitled" }; // (CRC)
	
	// Monte Carlo Controls
	uint32_t     numsteps   = 0;       // Total number of MC simulation steps to perform
	uint32_t     step       = 0;       // Current MC step
	uint32_t     corrtime   = 0;       // Number of steps between MC correlation times
	int          ptemp_freq = 0;

	double       move_factor           = 1.0;
	double       rot_factor            = 1.0;
	double       last_volume           = 0;    // NPT option
	double       volume_change_factor  = 0.25; // set default volume change factor (for NPT) 
	double       adiabatic_probability = 0,
	             gwp_probability       = 0,
	             insert_probability    = 0,
	             spinflip_probability  = 0,
	             volume_probability    = 0,
		         transfer_probability  = 0;

	// Path-Integral settings
	double      bead_perturb_probability = { 0 }; //  PI option--probability for move that changes bead configuration
	


	// io filenames
	char        dipole_output     [maxLine] = { 0 },
	            energy_output     [maxLine] = { 0 },
	            energy_output_csv [maxLine] = { 0 },
	            field_output      [maxLine] = { 0 },
	            frozen_output     [maxLine] = { 0 },
	            histogram_output  [maxLine] = { 0 },
	            insert_input      [maxLine] = { 0 },
	            pqr_input         [maxLine] = { 0 },
	            pqr_input_B       [maxLine] = { 0 },
	            pqr_output        [maxLine] = { 0 },
	            pqr_restart       [maxLine] = { 0 },
	            surf_output       [maxLine] = { 0 },
	            traj_input        [maxLine] = { 0 },
	            traj_output       [maxLine] = { 0 },
	            virial_output     [maxLine] = { 0 };
	    
	 

	// Observables
	double      temperature     = 0,
	            pressure        = 0,
	            free_volume     = 0,
	            total_energy    = 0,
	            N               = 0;
	int         fugacitiesCount = 0;
	double      fugacities[maxTokens] = {0};
	

	// Inter-node/Inter-system Data relay
	mpiData     mpi_data;


	// Simulation Flags & Data
	int         gwp               = 0;
	int         long_output       = 0; // Flag: signals request to print extended (%11.6f) coordinates
	int         parallel_restarts = 0; // Flag: signals that this run is a restart of a parallel job
	double      max_bondlength    = 0; // Bond threshold (re:output files)
	
	

	cavity_t   *** cavity_grid      = nullptr;
	int            cavity_bias      = 0,
	               cavity_grid_size = 0,
	               cavities_open    = 0;
	double         cavity_autoreject_repulsion = 0,
	               cavity_autoreject_scale     = 0,
	               cavity_radius               = 0,
	               cavity_volume               = 0;
	// Auto-reject Options --- first is in terms of sigma and only
	// applies to LJ; latter is in Angstroms and applies to all pairs
	int            cavity_autoreject          = 0;  // Flag: autoreject in terms of sigma--only applies to LJ
	int            cavity_autoreject_absolute = 0;  // Flag: autoreject in Angstroms and applies to all pairs

	// Parallel Tempering Options
	int            parallel_tempering = 0;
	double         max_temperature    = 0;
	ptemp_t      * ptemp              = nullptr;
	
	// Info regarding  periodic boundary and unit cell geometry
	int              wrapall = 1;  // Flag: wrap all option requested
	PeriodicBoundary pbc;          // Periodic boundary conditions of THIS system 
	
	
	// (P)RNG
	std::mt19937   mt_rand;
	std::uniform_real_distribution<double> dist{0,1};
	int            preset_seed_on = 0;  //for manual specification of random seeds
	unsigned int   preset_seed    = 0;  // datatype taken by std::mt19937.seed( unsigned int SEED)
	
	int            read_pqr_box_on = 0; //read box basis from pqr
	
	// Simulated Annealing
	int            simulated_annealing          = 0,
	               simulated_annealing_linear   = 0;
	double         simulated_annealing_schedule = 0,
	               simulated_annealing_target   = 0;

	// Spectre
	int            spectre            = 0;
	double         spectre_max_charge = 0,
	               spectre_max_target = 0;


	// Energy-corrections
	int            feynman_hibbs       = 0,
	               feynman_kleinert    = 0,
	               feynman_hibbs_order = 0;
	int            vdw_fh_2be          = 0; // Flag: 2BE method for polarvdw requested
	int            rd_lrc              = 1, // default rd LRC flag 
	               rd_crystal          = 0,
	               rd_crystal_order    = 0;

	// uVT Fugacity Functions
	int            h2_fugacity     = 0,
	               co2_fugacity    = 0,
	               ch4_fugacity    = 0,
	               n2_fugacity     = 0,
	               user_fugacities = 0;


	// Force-field Options
	int            rd_only         = 0,
	               rd_anharmonic   = 0;
	double         rd_anharmonic_k = 0,
	               rd_anharmonic_g = 0;
	int            c6_mixing       = 0,
	               damp_dispersion = 0,
	               disp_expansion_mbvdw    = 0,
	               use_dreiding            = 0,
	               extrapolate_disp_coeffs = 0,
	               halgren_mixing          = 0,
	         
	               midzuno_kihara_approx   = 0,
	               schmidt_ff              = 0,
	               waldmanhagler           = 0;

	bool           using_axilrod_teller    = false,
	               using_lj_buffered_14_7  = false,
	               using_disp_expansion    = false,
	               use_sg                  = false;

	// ES Options
	int            wolf                  = 0,
	               ewald_alpha_set       = 0,
	               ewald_kmax            = ewald_kmax_default,
	               polar_ewald_alpha_set = 0;
	double         ewald_alpha           = ewald_alpha_default,
		           polar_ewald_alpha     = ewald_alpha_default;
	
	
	// Thole Options
	int            polarization = 0; // Flag signaling that polarization calculation has been requested
	int            polarvdw              = 0,
		           polarizability_tensor = 0,
		           cdvdw_exp_repulsion   = 0,
		           cdvdw_sig_repulsion   = 0,
		           cdvdw_9th_repulsion   = 0;
	int            iterator_failed       = 0; //flag set when iterative solver fails to converge (when polar_precision is used)
	int            polar_iterative       = 0,
		           polar_ewald           = 0,
		           polar_ewald_full      = 0,
		           polar_zodid           = 0,
		           polar_palmo           = 0,
		           polar_rrms            = 0,
		           polar_gs              = 0;

	
	int            polar_gs_ranked         = 0;  // Flag indicating if the ranked gauss-seidell algorithm will be used in polar calculations.
	int            polar_sor               = 0,
	               polar_esor              = 0,
	               polar_max_iter          = 0,
	               polar_wolf              = 0,
	               polar_wolf_full         = 0,
	               polar_wolf_alpha_lookup = 0;
	double         polar_wolf_alpha        = 0,
	               polar_gamma             = 1.0,
	               polar_damp              = 0,
	               field_damp              = 0,
	               polar_precision         = 0,
	             * polar_wolf_alpha_table         = nullptr,
	               polar_wolf_alpha_lookup_cutoff = wolf_alpha_lookup_cutoff_default;
	int            polar_wolf_alpha_table_max     = 0; //stores the total size of the array double polar_wolf_alpha_table[]
	int            damp_type;
	double      ** A_matrix       = nullptr;             // A matrix (Thole polarization) 
	double      ** B_matrix       = nullptr;             // B matrix (Thole polarization)
	double         C_matrix[3][3] = {0}; // Polarizability tensor 

	vdw_t        * vdw_eiso_info = nullptr; //keeps track of molecule vdw self energies
	

	//misc
	double         scale_charge         = 1.0;
	int            independent_particle = 0;

	// insertions from a separate linked list
	int            num_insertion_molecules   = 0;       // the number of elements found in the Molecule::insertion_molecules and Molecule::insertion_molecules_array
	Molecule     * insertion_molecules       = nullptr; // linked list of molecules to be randomly inserted
	Molecule    ** insertion_molecules_array = nullptr; // array providing direct access to the linked list holding molecules from which insertions are chosen


	// quantum rotation stuff
	int            quantum_rotation                  = 0, 
	               quantum_rotation_hindered         = 0,
	               quantum_rotation_l_max            = 0,  
	               quantum_rotation_level_max        = 0, 
	               quantum_rotation_phi_max          = 0,
	               quantum_rotation_sum              = 0,
	               quantum_rotation_theta_max        = 0,
	               quantum_vibration                 = 0;
	double         quantum_rotation_B                = 0,
	               quantum_rotation_hindered_barrier = 0;
	

	// histogram stuff
	grid_t       * grids            = nullptr;
	int            calc_hist        = 0; // flag to calculate a 3D histogram 
	double         hist_resolution  = 0;
	int            n_histogram_bins = 0;

	// atom array	
	int            natoms = 0;
	Atom        ** atom_array     = nullptr;
	Molecule    ** molecule_array = nullptr;

	//replay option
	int            calc_pressure    = 0;
	double         calc_pressure_dv = 0;

	Molecule            * molecules        = nullptr;

	nodestats_t         * nodestats        = nullptr;
	avg_nodestats_t     * avg_nodestats    = nullptr;
	observables_t       * observables      = nullptr;
	avg_observables_t   * avg_observables  = nullptr;

	// Linked list head that will keep track of separate average-observables for each sorbate in the system.
	int                  sorbateCount      = 0;        // Number of sorbates in the system.
	sorbateInfo_t      * sorbateInfo       = nullptr;  // stores an array of sorbate Info
	int                  sorbateInsert     = 0;        // which sorbate was last inserted
	sorbateAverages_t  * sorbateGlobal     = nullptr;  // where the global average is stored

	checkpoint_t       * checkpoint        = nullptr;

	FILE               * fp_energy         = nullptr;
	FILE               * fp_energy_csv     = nullptr;
	FILE               * fp_field          = nullptr;
	FILE               * fp_histogram      = nullptr;
	FILE               * fp_frozen         = nullptr;
	FILE               * fp_traj_replay    = nullptr;
	FILE               * fp_surf           = nullptr;



	 // surface fitting options
	///////////////////////////////////////////////////////////

	fileNode_t     fit_input_list;
	int            ee_local                   = 0; // Exhaustive Enumeration option
	double         fit_best_square_error      = 0,
	               fit_max_energy             = 0,
	               fit_schedule               = 0,
	               fit_start_temp             = 0;
	int            fit_boltzmann_weight       = 0;
	int            surf_decomp                = 0,
	               surf_fit_arbitrary_configs = 0,
	               surf_qshift_on             = 0,
	               surf_weight_constant_on    = 0,
	               surf_global_axis_on        = 0;
	double         surf_ang = 0,
	               surf_inc = 0,
	               surf_min = 0,
	               surf_max = 0;
	int            surf_descent  = 0,
	               surf_preserve = 0,
	               surf_preserve_rotation_on     = 0;
	double         surf_preserve_rotation_alpha1 = 0,
	               surf_preserve_rotation_alpha2 = 0,
	               surf_preserve_rotation_beta1  = 0,
	               surf_preserve_rotation_beta2  = 0,
	               surf_preserve_rotation_gamma1 = 0,
	               surf_preserve_rotation_gamma2 = 0;
	double         range_eps = 0,
	               range_sig = 0,
	               step_eps  = 0,
	               step_sig  = 0;
	int            surf_print_level      = 0; // sets the amount of output (1-6) that correspond to the nested loops in surface.c
	int            surf_scale_alpha_on   = 0,
	               surf_scale_epsilon_on = 0,
	               surf_scale_omega_on   = 0,
	               surf_scale_sigma_on   = 0,
	               surf_scale_pol_on     = 0,
	               surf_scale_q_on       = 0,
	               surf_scale_r_on       = 0,
	               surf_scale_c6_on      = 0,
	               surf_scale_c8_on      = 0,
   	               surf_scale_c10_on     = 0;
	double         surf_scale_epsilon    = 0,
	               surf_scale_r          = 0,
	               surf_scale_omega      = 0,
	               surf_scale_sigma      = 0,
	               surf_scale_q          = 0,
	               surf_scale_alpha      = 0,
	               surf_scale_pol        = 0,
	               surf_scale_c6         = 0,
	               surf_scale_c8         = 0,
	               surf_scale_c10        = 0;
	int            surf_virial           = 0;
	double         surf_quadrupole       = 0,
	               surf_weight_constant  = 0;

};

#endif // __SYSTEM_H
