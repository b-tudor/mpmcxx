#pragma once

#include <stdint.h>
#include <stdio.h>
#include <map>
#include <string>
#include <vector>

#include "constants.h"
#include "System.h"
#include "Vector3D.h"






class SimulationControl
{

public:

	System sys; // Template system for creating system images. Can be used to retrieve certain info common to all systems (e.g. temperature)
	int nSys; // Trotter number for Path Integral runs, i.e. number of beads/images being used to represent each quantum object
	int PI_trial_chain_length; // PI option--when perturbing COM configuration, the number of beads to move at a time

	SimulationControl( char * inFilename, int P, bool write_PI_frames, char *fname );
	~SimulationControl();
	bool runSimulation();
	void initializeSimulationObjects();


	inline bool PI_xyz_corrtime_frames_requested() { return write_PI_frames; }
	
	
private:

	typedef struct  _Boltzmann_Factor_contributor {
		double init;
		double trial;
		double current;
		inline double change() { return trial - init; }
	} BoltzFactor_contributor;
	typedef struct _PI_NVT_BFContributors {
		BoltzFactor_contributor potential;       // Potential energy in Kelvin
		BoltzFactor_contributor chain_mass_len2; // kg * m^2
		BoltzFactor_contributor orient_mu_len2;  // kg * m^2
	} PI_NVT_BFContributors;
	PI_NVT_BFContributors BFC;
	

	bool write_PI_frames; // do we want to write the PI frames for visualization?
	char *PI_frames_filename; // filename for PI xyz images (image, as in 'bead') file

	
	std::vector<System *> systems;                  // collection of systems for use with multi-system coordinated ensembles (Gibbs, Path Integral)
	std::vector<Vector3D> orientations;             // orientations to be applied to a bead-representation of a molecule during bead perturbation

	std::vector<std::string> pqr_input_filenames;   // Input geometry filenames for each PI system
	std::vector<std::string> pqr_restart_filenames; // Output geometry filenames for each PI system
	std::vector<std::string> pqr_final_filenames;   // Output geometry filenames for each PI system--final output upon Sim normal exit
	std::vector<std::string> traj_filenames;        // Output trajectory files for each PI system
	

	// Structures used to hold data about linear molecules (H2, CO2) that will be used to orient them when they are 
	// subject to a perturbation in their Path Integral representation. These values will also be needed to compute
	// the PI energy associated with a molecule's distribution of orientations. 
	typedef struct _molecular_metadata {
		int    orientation_site;
		double bond_length;
		double reduced_mass;
	} molecular_metadata;
	static std::vector<molecular_metadata> sorbate_data;
	static std::map<std::string, uint32_t> sorbate_data_index;
	
	

	// arrays to accumulate energies across different systems, when using MPI (or for easy access on single threads)
	double * rd_energies;           // Repulsion/dispersion energies of all individual systems, indexed by system #
	double * coulombic_energies;    //            Coulombic energies of all individual systems, indexed by system #
	double * polarization_energies; //         Polarization energies of all individual systems, indexed by system #
	double * vdw_energies;          //        Van der Waals energies of all individual systems, indexed by system #
	
		
	// read_config() parses a simulation input file and populates the Simulation Controller
	// with the data found therein. Data read in is largely unvalidated.
	void read_config(char *inFilename );
	

	// Read each input line from the input file and set individual system flags accordingly. Minimal
	// error checking performed, mostly in the form of checking for malformed command sequences. 
	bool process_command( char token[maxTokens][maxLine] );
	

	// The following functions validate the system, checking for missing,
	// invalid or incompatible options & settings.
	bool check_system();
	bool check_mc_options();
	bool check_spectre_options();
	bool check_io_files_options();
	bool check_feynman_hibbs_options();
	bool check_simulated_annealing_options();
	bool check_hist_options();
	bool check_polarization_options();
	bool check_qrot_options();


	
	void initialize_PI_NVT_Systems();        // Allocate data structures and initialize system data for multi-system ensembles
	void backup_observables_ALL_SYSTEMS();   // backup observables for all systems in the system vector as well as the sys variable
	void backup_observables_SYS_VECTOR();    // backup observables only for the systems in the system vector
	void restore_observables_ALL_SYSTEMS();  // restore observables for all systems in the system vector as well as the sys variable
	void restore_observables_SYS_VECTOR();   // restore observables only for the systems in the system vector

	// Standard 
	inline bool mc()                    { return sys.mc(); } // SimulationControl wrapper function for System.mc()
	inline bool surface()               { return false; }    // sys.surface(); } // Simulation Control wrapper for System.surface()
	inline bool surface_fit()           { return false; }    // sys.surface_fit(); }
	inline bool surface_fit_arbitrary() { return false; }    // sys.surface_fit_arbitrary(); }
	inline bool calculate_te()          { return false; }    // sys.calculate_te (); }
	inline bool replay_trajectory()     { return false; }    // sys.replay_trajectory(); }
	
	

	// SimulationControl.PathIntegral.cpp
	bool   PI_nvt_mc();
	void   do_PI_corrtime_bookkeeping();
	bool   check_PI_options();
	
	void   assert_all_perturb_targets_exist(); // Check all systems for an existing target molecule for the current MC move
	double PI_calculate_energy(); // Calculate the total energy for the aggregate PI system for the current point in phase space
	double PI_calculate_potential(); // Calculate the potential energy for the aggregate PI system
	double PI_calculate_kinetic(); // Calculate the kinetic energy for the aggregate PI system
	double PI_chain_mass_length2_ENTIRE_SYSTEM(); // Return sum total of all mass-weighted chain length measures for all the molecules in aggregate system
	double PI_chain_mass_length2(); // Return COM PI "polymer" mass-weighted chain length measure for molecule targeted by MC move
	double PI_chain_mass_length2(std::vector<Molecule*> &m); // Return COM PI "polymer" mass-weighted chain length for molecule represented by m
	double PI_orientational_mu_length2_ENTIRE_SYSTEM(); // Return "orientational chain" length measure for all molecules in the system
	double PI_orientational_mu_length2(); // Return "orientational chain" length measure for the MC move's target molecule
	double PI_orientational_mu_length2(std::vector<Vector3D*> &o); // Return "orientational chain" length measure for given vector of orientations
	double PI_NVT_boltzmann_factor( PI_NVT_BFContributors BF );
	void   PI_calc_system_mass(); // Calculate the system mass from a PI image and copy the resulting values into the PI observables struct (sys.observables)
	int    PI_pick_NVT_move();
	void   PI_make_move( int move );
	void   PI_flip_spin();
	void   PI_displace();
	void   PI_perturb_beads();
	void   PI_perturb_bead_COMs_ENTIRE_SYSTEM();
	void   PI_perturb_bead_COMs();      // perturb the user-specified number of beads
	void   PI_perturb_bead_COMs(int n); // specify number of beads to perturb
	void   PI_perturb_beads_orientations();
	void   generate_orientation_configs();
	void   generate_orientation_configs(unsigned int start, unsigned int end, unsigned int P, unsigned int numBeads, double b2, double uMkT );
	void   apply_orientation_configs();
	void   restore_PI_systems();
	void   average_current_observables_into_PI_avgObservables();
	
	void   write_PI_frame();

	static void   add_orientation_site_entry(const char *id, const int site);
	static int    get_orientation_site( std::string id );
	static void   add_bond_length_entry(const char *id, const double bond_length);
	static double get_bond_length(std::string molecule_id);
	static void   add_reduced_mass_entry(const char *id, const double reduced_mass);
	static double get_reduced_mass(std::string molecule_id);


	// SimulationControl.Gibbs.cpp
	bool Gibbs_mc();
	bool check_Gibbs_options();
	void initialize_Gibbs_systems();
	static void boltzmann_factor_NVT_Gibbs(System &sys1, double initEnergy_1, double finalEnergy_1, System &sys2, double initEnergy_2, double finalEnergy_2);

};