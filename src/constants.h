#pragma once

constexpr int maxLine   = 512;
constexpr int maxID     =  25;
constexpr int maxTokens =  10;

// MPI on/off switch
// #define MPI



// Physical Constants
const double pi        = 3.141592653589793238462643383279502884L;
const double h         = 6.626068e-34;    // Planck's constant in J s 
const double hBar      = 1.054571e-34;    // above divided by 2pi in J s 
const double c_hBar    = 7.63822291e-12;  // Ks //HBAR is already taken to be in Js
const double hBar2     = 1.11211999e-68;  // hBar^2 in (Js)^2
const double hBar4     = 1.23681087e-136; // hBar^4 in (Js)^4
const double half_hBar = 3.81911146e-12;  // hBar/2 in Ks
const double kB        = 1.3806503e-23;   // Boltzmann's constant in J/K 
const double kB2       = 1.90619525e-46;  // kB^2
const double NA        = 6.0221415e23;    // Avogadro's number
const double c         = 2.99792458e8;    // speed of light in vacuum in m/s



// conversion factors 
const double au2invseconds     = 4.13412763705666648752113572754445220741745180640e16; //s^-1 a.u.^-1
const double AU2ANGSTROM       = 0.529177249;    // convert from Bohr radii to angstroms 
const double METER2ANGSTROM    = 1.0e10;         // convert from meters to angstroms 
const double ANGSTROM2METER    = 1.0e-10;        // convert from meters to angstroms 
const double M2A2              = 1.0e20;
const double M2A4              = 1.0e40;
const double HARTREE2KELVIN    = 3.15774655e5;   // convert from Hartrees to Kelvin 
const double E2REDUCED         = 408.7816;       // convert from e to sqrt(K*A) 
const double ATM2REDUCED       = 0.0073389366;   // convert from atm to K/A^3 
const double ATM2PASCALS       = 101325.0;       // convert from atm to Pascals 
const double ATM2PSI           = 14.6959488;     // convert from atm to psi 
const double A32CM3            = 1.0e-24;        // convert from A^3 to cm^3 
const double AMU2KG            = 1.66053873e-27; // convert amu's to kg 
const double DEBYE2SKA         = 85.10597636;    // convert from Debye to sqrt(KA)*A 
const double EV2K              = 1.160444e4;     // convert eV to K 
const double K2WN              = 0.695039;       // convert K to cm^-1 
const double KoverANGcubed2ATM = 136.259;        // convert K/A^3 to ATM 
const double LITER2A3          = 1.0e27;         // convert Liters to Ang^3 
const double GASCONSTANT       = 0.8205746;

const double OneOverSqrtPi     = 0.5641895835477562869480794515607725858440506293289988;
const double SqrtPi            = 1.77245385091;
const double twoPi             = 2.0L * pi;

const double MAX_ITERATION_COUNT = 128;
const double MAXVALUE            = 1.0e40;
const double SMALL_dR            = 1.0e-12;

const double FEYNMAN_KLEINERT_TOLERANCE = 1.0e-12;  // tolerance in A^2



// Enumerated sets //////////////////////////////////////////////////////////

enum { 
	REAL, 
	IMAGINARY 
};
enum { 
	DAMPING_OFF, 
	DAMPING_LINEAR, 
	DAMPING_EXPONENTIAL 
};
enum { 
	NUCLEAR_SPIN_PARA, 
	NUCLEAR_SPIN_ORTHO
};
enum {
	ENSEMBLE_UVT,
	ENSEMBLE_NVT,
	ENSEMBLE_SURF,
	ENSEMBLE_SURF_FIT,
	ENSEMBLE_NVE,
	ENSEMBLE_TE,
	ENSEMBLE_NPT,
	ENSEMBLE_REPLAY,
	ENSEMBLE_PATH_INTEGRAL_NVT,
	ENSEMBLE_NVT_GIBBS
};
enum {
	MOVETYPE_INSERT,
	MOVETYPE_REMOVE,
	MOVETYPE_DISPLACE,
	MOVETYPE_ADIABATIC,
	MOVETYPE_SPINFLIP,
	MOVETYPE_VOLUME,
	MOVETYPE_PERTURB_BEADS
};





// Function return codes //////////////////////////////////////////////////////
const bool ok   = true;
const bool fail = false;



// Error Codes for exceptions and program exit //////////////////////////////
const int exception_ok                  =   100; // Program finished normally, albeit via the exception mechanism.
const int internal_error                =   101; // Error in code. That is, this is in no way a result of user input.
const int invalid_monte_carlo_move      =   102; // Error in code: invalid Monte Carlo move code generated/requested.
const int invalid_quaternion_mode       =   103; // Error in code: invalid Quaternion construction mode code specified.
const int interrupt_signal_received     =   104; // Program received an external request to end execution.

const int fopen_fail_read               =  1000; // fopen failed (read mode)
const int fopen_fail_write              =  1001; // fopen failed (write mode)
const int fopen_fail_append             =  1002; // fopen failed (append mode)
const int fopen_fail_unknown            =  1003; // fopen failed (unknown mode)
const int unknown_file_error            =  1004; // unspecified failed disk operation 
const int null_file_ptr_error           =  1005; // null file pointer sent to a write routine

const int memory_request_fail           =  2000; // memory was requested of system, but request failed
const int memory_request_invalid        =  2001; // memory request is nonsensical

const int invalid_input                 =  3000; // nonsensical/unreadable input, settings (input file)
const int no_molecules_in_system        =  3001; // the pqr file for the physical system has no molecules in a system that requires them

const int invalid_setting               =  4000; // input correctly parsed, but requested value/option not valid for setting
const int invalid_ensemble              =  4001; // specified ensemble code is not recognized (or was never set?)
const int incompatible_settings         =  4002; // two or more options specified cannot be used with each other
const int missing_setting               =  4003; // a required input setting was not provided
const int unsupported_setting           =  4004; // a valid feature has been requested, but it is not yet supported

const int missing_required_datum        =  6000; // e.g. z coordinate of an atom
const int invalid_datum                 =  6001; // nonsensical/unreadale input, data (geometry file)
const int molecule_wo_atoms             =  6002; // a molecule in the geometry file contains no atoms
const int invalid_box_dimensions        =  6004; // simulation box was an invalid size (vol or cutoff probably < 0 )
const int incongruent_bead_states       =  6005; // the beads representing a molecule have become dissimilar in an invalid way (missing beads/corresponding beads are different elements/etc)

const int surface_fit_module_fail       = 10001; // surface fitting function failed to complete

const int infinite_energy_calc          = 11000; // energy calculation was not finite
const int lapack_error                  = 11001; // lapack function call returned error 
const int attempted_singular_mtx_inv    = 11002; // function attempted to invert a matrix whose determinant is 0

const int invalid_MPI_size_for_PI       = 12000; // Path Integral "bead chains" require at least three systems

const int internal_err_invalid_mc_move  = 20000; // Moves in different systems within a coordinated ensemble are incompatible