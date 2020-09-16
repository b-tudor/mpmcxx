#include "System.h"

#include <cstring>
// (OS Dependent) Timing Includes
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32)
#	include <time.h>
#	include <windows.h>
#else
#	include <sys/time.h>
#endif

#include "Atom.h"
#include "Molecule.h"
#include "Output.h"
#include "SafeOps.h"
#include "UsefulMath.h"


// MPI Includes
#ifdef _MPI
#	include <mpi.h>
#endif
extern int rank, size;
extern bool mpi;




int System::open_files()
{

	if( energy_output[0] ) {
		fp_energy = SafeOps::openFile( energy_output, "w", __LINE__, __FILE__ );
		fprintf( fp_energy, "#step #energy #coulombic #rd #polar #vdw #kinetic #kin_temp #N #spin_ratio #volume #core_temp\n" );
	}

	if( energy_output_csv[0] ) {
		fp_energy_csv = SafeOps::openFile( energy_output_csv, "w", __LINE__, __FILE__ );
		//filecheck(system->file_pointers.fp_energy_csv,system->energy_output_csv,WRITE);
		fprintf( fp_energy_csv,	"#step,#energy,#coulombic,#rd,#polar,#vdw,#kinetic,#kin_temp,#N,#spin_ratio,#volume,#core_temp\n" );
	}

	// if we're just calculating energy or replaying a trajectory, we need no other output files
	if( ensemble == ENSEMBLE_REPLAY || ensemble == ENSEMBLE_TE ) 
		return 0;

	if( histogram_output[0] ) {
		fp_histogram = SafeOps::openFile( histogram_output, "w", __LINE__, __FILE__ );
	}

	if( frozen_output[0] ) {
		fp_frozen = SafeOps::openFile( frozen_output, "w", __LINE__, __FILE__ );
		//go ahead and write the frozen lattice configuration now
		if (fp_frozen) {
			write_frozen(fp_frozen);
			fclose(fp_frozen);
			fp_frozen = nullptr;
		}
	}

	return 0;
}




void System::close_files()
{
	if( fp_energy      ) fclose( fp_energy      );
	if( fp_energy_csv  ) fclose( fp_energy_csv  );
	if( fp_histogram   ) fclose( fp_histogram   );
	if( fp_traj_replay ) fclose( fp_traj_replay );
	if( fp_surf        ) fclose( fp_surf        );

	fp_energy      = nullptr;
	fp_energy_csv  = nullptr;
	fp_histogram   = nullptr;
	fp_traj_replay = nullptr;
	fp_surf        = nullptr;
}




void System::write_frozen( FILE *fp )
{
	if( !fp )
		throw null_file_ptr_error;

	int numatoms;
	int numbonds;

	numatoms = count_frozen();
	numbonds = calculate_bonds();
	
	rewind( fp );
	fprintf( fp, "# OpenDX format coordinate file for frozen atoms\n" );
	fprintf( fp, "object 1 class array type float rank 1 shape 3 items %d data follows\n", numatoms );
	print_frozen_coords( fp );
	fprintf( fp, "object 2 class array type int rank 1 shape 2 items %d data follows\n", numbonds );
	print_frozen_bonds( fp );
	fprintf( fp, "attribute \"element type\" string \"lines\"\n" );
	fprintf( fp, "attribute \"ref\" string \"positions\"\n" );
	fprintf( fp, "object 3 class array type float rank 0 items %d data follows\n", numatoms );
	print_frozen_masses( fp );
	fprintf( fp, "attribute \"dep\" string \"positions\"\n" );
	fprintf( fp, "object 4 class array type float rank 1 shape 3 items %d data follows\n", numatoms );
	print_frozen_colors( fp );
	fprintf( fp, "object \"irregular positions irregular connections\" class field\n" );
	fprintf( fp, "component \"positions\" value 1\n" );
	fprintf( fp, "component \"connections\" value 2\n" );
	fprintf( fp, "component \"data\" value 3\n" );
	fprintf( fp, "component \"colors\" value 4\n" );
	fprintf( fp, "end\n" );

}




int System::count_frozen()
{
	int         count = 0;
	Molecule  * mol   = nullptr;
	Atom      * atom  = nullptr;

	for( mol = molecules; mol; mol=mol->next ){
		if(mol->frozen){
			for(atom = mol->atoms; atom; atom=atom->next)
				count++;
		}
	}
	return count;
}




void System::print_frozen_coords( FILE *fp )
{
	if( !fp )
		throw null_file_ptr_error;

	Molecule * mol  = nullptr;
	Atom     * atom = nullptr;

	for( mol = molecules; mol; mol=mol->next ) {
		if( mol->frozen ) {
			for( atom = mol->atoms; atom; atom=atom->next ) {
				fprintf( fp, "%f %f %f\n", atom->pos[0], atom->pos[1], atom->pos[2] );
			}
		}
	}
}




void System::print_frozen_bonds( FILE *fp )
{
	if( !fp )
		throw null_file_ptr_error;

	int inner_index=0,outer_index=0;
	int bonds=0;
	
	Molecule  * mol   = nullptr;
	Atom      * atom  = nullptr,
	          * atom2 = nullptr;

	for( mol = molecules; mol; mol = mol->next ) {
		if( mol->frozen ) {
			for(atom=mol->atoms; atom; atom=atom->next, inner_index++){
				if( atom->next ) {
					for( atom2 = atom->next, outer_index = inner_index+1; atom2; atom2 = atom2->next, outer_index++ ) {
						if( bondlength_check(atom, atom2) ) {
							bonds++;
							fprintf( fp, "%d %d\n", inner_index, outer_index );
						}
					}
				}
			}
		}
	}
}




void System::print_frozen_masses( FILE *fp )
{
	if( !fp )
		throw null_file_ptr_error;

	Molecule  * mol  = nullptr;
	Atom      * atom = nullptr;

	for( mol = molecules; mol; mol = mol->next ) {
		if( mol->frozen ) {
			for( atom=mol->atoms; atom; atom = atom->next ){
				fprintf( fp, "%f\n", atom->mass );
			}
		}
	}
}




void System::print_frozen_colors( FILE *fp )
{
	if( !fp )
		throw null_file_ptr_error;

	Molecule  * mol  = nullptr;
	Atom      * atom = nullptr;
	double      mass = 0;
	
	const std::string COLOR_H  ( "0.2 0.2 0.2" );
	const std::string COLOR_C  ( "0.1 0.5 0.1" );
	const std::string COLOR_N  ( "0.2 0.2 1.0" );
	const std::string COLOR_O  ( "1.0 0.0 0.0" );
	const std::string COLOR_XXX( "0.1 0.1 0.1" );

	for( mol = molecules; mol; mol=mol->next ) {
		if( mol->frozen ) {
			for( atom=mol->atoms; atom; atom=atom->next ) {
				mass=atom->mass;
				if( mass < 1.1 )
					fprintf( fp, "%s\n", COLOR_H.c_str() );
				else if( mass < 12.2 )
					fprintf( fp, "%s\n", COLOR_C.c_str() );
				else if( mass < 14.1 )
					fprintf( fp, "%s\n", COLOR_N.c_str() );
				else if(mass < 16.1)
					fprintf( fp, "%s\n", COLOR_O.c_str() );
				else
					fprintf( fp, "%s\n", COLOR_XXX.c_str() );
			}
		}
	}

}




void System::append_observables_to_output_file() {
	write_observables(fp_energy, observables, temperature);
}
void System::write_observables(FILE* fp, observables_t* obs, double core_temp) {
	
	if( !fp )
		throw null_file_ptr_error;

	fprintf( fp, "%d %f %f %f %f %f %f %f %f %f %f %f", 
		step,
		obs->energy, 
		obs->coulombic_energy, 
		obs->rd_energy, 
		obs->polarization_energy, 
		obs->vdw_energy, 
		obs->kinetic_energy, 
		obs->temperature, 
		obs->N, 
		obs->spin_ratio, 
		obs->volume,
		core_temp );
	fprintf( fp, "\n" );
	fflush( fp );
}




void System::append_observables_to_csv_file() {
	write_observables_csv(fp_energy_csv, observables, temperature);
}
void System::write_observables_csv( FILE *fp, observables_t * obs, double core_temp) {
	
	if( !fp )
		throw null_file_ptr_error;

	fprintf( fp, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", 
		step,
		obs->energy, 
		obs->coulombic_energy, 
		obs->rd_energy, 
		obs->polarization_energy, 
		obs->vdw_energy, 
		obs->kinetic_energy, 
		obs->temperature, 
		obs->N, 
		obs->spin_ratio, 
		obs->volume,
		core_temp );
	fprintf( fp, "\n" );
	fflush( fp );
}




int System::display_averages() {
	return display_averages("");
}
int System::display_averages(int sysNum) {
	char sysID[16] = { 0 };
	sprintf(sysID, "_%d", sysNum);
	return display_averages(sysID);
}
int System::display_averages( const char *sysID ) {

	int i;
	char linebuf[maxLine];
	avg_observables_t &avg =*avg_observables;

	if(avg.boltzmann_factor > 0.0) {
		sprintf( linebuf, "OUTPUT%s: BF = %.5lg +- %.5lg\n", sysID, avg.boltzmann_factor, avg.boltzmann_factor_error );
		Output::out1(linebuf);
	}

	if(avg.acceptance_rate > 0.0) {
		sprintf(
			linebuf, 
			"OUTPUT%s: AR = %.5lf (%.5lf I/ %.5lf R/ %.5lf D",
			sysID,
			avg.acceptance_rate,
			avg.acceptance_rate_insert, 
			avg.acceptance_rate_remove,
			avg.acceptance_rate_displace
		);
		
		if(avg.acceptance_rate_adiabatic > 0.0) 
			sprintf(&linebuf[strlen(linebuf)], "/ %.5lf A", avg.acceptance_rate_adiabatic);

		if(avg.acceptance_rate_spinflip > 0.0) 
			sprintf(&linebuf[strlen(linebuf)], "/ %.5lf S", avg.acceptance_rate_spinflip);

		if(avg.acceptance_rate_volume > 0.0) 
			sprintf(&linebuf[strlen(linebuf)], "/ %.5lf V", avg.acceptance_rate_volume);

		if(avg.acceptance_rate_ptemp > 0.0) 
			sprintf(&linebuf[strlen(linebuf)], "/ %.5lf PT", avg.acceptance_rate_ptemp);
		
		if (avg.acceptance_rate_beadPerturb > 0.0) 
			sprintf(&linebuf[strlen(linebuf)], "/ %.5lf BEAD", avg.acceptance_rate_beadPerturb);
		
		sprintf( &linebuf[strlen(linebuf)], "\n");
		Output::out1(linebuf);
	}

	//print node's current temperature if doing SA or PT
	if( simulated_annealing ) {
		sprintf( linebuf, "OUTPUT%s: Simulated Annealing Temperature = %.5f K\n", sysID, temperature );
		Output::out1(linebuf);
	}

	if( avg.cavity_bias_probability > 0.0 ) {
		sprintf( linebuf, "OUTPUT%s: Cavity bias probability = %.5f +- %.5f\n", 
			sysID,
			avg.cavity_bias_probability, avg.cavity_bias_probability_error);
		Output::out1(linebuf);
	}


	if (gwp)
		sprintf( linebuf, "OUTPUT%s: total energy = %.5lf +- %.5lf eV\n", sysID, avg.energy / EV2K, avg.energy_error / EV2K );
	else if (ensemble == ENSEMBLE_PATH_INTEGRAL_NVT) {
		sprintf( linebuf, "OUTPUT%s: total energy          = %.5lf +- %.5lf K\n", sysID, avg.energy, avg.energy_error );
		sprintf( linebuf+strlen(linebuf), "OUTPUT%s: total energy (virial) = %.5lf +- %.5lf K\n", sysID, avg.energy, avg.energy_error );
	} 
	else
		sprintf( linebuf, "OUTPUT%s: potential energy = %.5lf +- %.5lf K\n", sysID, avg.energy, avg.energy_error );
	Output::out1(linebuf);


	if( avg.coulombic_energy != 0.0 ) {
		if( gwp )
			sprintf( linebuf, "OUTPUT%s: electrostatic energy = %.5lf +- %.5lf eV\n",
				sysID, 
				avg.coulombic_energy/EV2K, avg.coulombic_energy_error/EV2K);
		else
			sprintf( linebuf, "OUTPUT%s: electrostatic energy = %.5lf +- %.5lf K\n",
				sysID, 
				avg.coulombic_energy, avg.coulombic_energy_error);
		Output::out1(linebuf);
	}


	if(avg.rd_energy != 0.0) {
		sprintf( linebuf, "OUTPUT%s: repulsion/dispersion energy = %.5lf +- %.5lf K\n", 
			sysID,
			avg.rd_energy, avg.rd_energy_error);
		Output::out1(linebuf);
	}

	if( avg.polarization_energy != 0.0 ) {
		sprintf( linebuf, "OUTPUT%s: polarization energy = %.5f +- %.5f K", 
			sysID,
			avg.polarization_energy, avg.polarization_energy_error );
		Output::out1(linebuf);

		if( avg.dipole_rrms_error != 0.0    &&    polar_rrms ) {
			sprintf( linebuf, " (iterations = %.1f +- %.1f rrms = %e +- %e)", 
				avg.polarization_iterations, avg.polarization_iterations_error, 
				avg.dipole_rrms, avg.dipole_rrms_error);
			Output::out1(linebuf);
		}
		else if( avg.polarization_iterations != 0.0 ) {
			sprintf( linebuf, " (iterations = %.1f +- %.1f)", 
				avg.polarization_iterations, avg.polarization_iterations_error);
			Output::out1(linebuf);
		}

		Output::out("\n");
	}

	#ifdef VDW
		sprintf( linebuf, "OUTPUT%s: (coupled-dipole) vdw energy = %.5f +- %.5f K\n", 
			sysID,
			avg.vdw_energy, avg.vdw_energy_error);
		Output::out1(linebuf);
	#endif

	if(avg.kinetic_energy > 0.0) {
		if( gwp ) {

			sprintf( linebuf, "OUTPUT%s: kinetic energy = %.5lf +- %.5lf eV\n", 
				sysID,
				avg.kinetic_energy/EV2K, 
				avg.kinetic_energy_error/EV2K
			);
			Output::out1(linebuf);

		} else {
			sprintf( linebuf, 
				"OUTPUT%s: kinetic energy = %.5lf +- %.5lf K\n", 
				sysID,
				avg.kinetic_energy, 
				avg.kinetic_energy_error
			);
			Output::out1(linebuf);
		}

		sprintf( linebuf, 
			"OUTPUT%s: kinetic temperature = %.5lf +- %.5lf K\n", 
			sysID,
			avg.temperature, 
			avg.temperature_error
		);
		Output::out1(linebuf);
	}

	sprintf( linebuf, "OUTPUT%s: N = %.5lf +- %.5lf molecules\n", sysID, avg.N, avg.N_error);
	Output::out1(linebuf);

	if( sorbateCount == 1 ) { //all based on calculations which assume only one type of sorbate

		sprintf( linebuf, "OUTPUT%s: density = %.5f +- %.5f g/cm^3\n", sysID, avg.density, avg.density_error );
		Output::out1(linebuf);

		if( avg.pore_density != 0.0   &&   ensemble != ENSEMBLE_NPT ){
			sprintf( linebuf, "OUTPUT%s: pore density = %.5f +- %.5f g/cm^3\n", sysID, avg.pore_density, avg.pore_density_error);
			Output::out1(linebuf);
		}
		if(avg.percent_wt > 0.0 ) {
			sprintf( linebuf, "OUTPUT%s: wt %% = %.5f +- %.5f %%\n", sysID, avg.percent_wt, avg.percent_wt_error);
			Output::out1(linebuf);
			sprintf( linebuf, "OUTPUT%s: wt %% (ME) = %.5f +- %.5f %%\n", sysID, avg.percent_wt_me, avg.percent_wt_me_error);
			Output::out1(linebuf);
		}
		if(avg.excess_ratio > 0.0) {
			sprintf( linebuf, "OUTPUT%s: excess adsorption ratio = %.5f +- %.5f mg/g\n", sysID, avg.excess_ratio, avg.excess_ratio_error);
			Output::out1(linebuf);
		}
		if(  (avg.qst > 0.0)   &&   std::isfinite(avg.qst)  ) {
			sprintf( linebuf, "OUTPUT%s: qst = %.5lf kJ/mol\n", sysID, avg.qst);
			Output::out1(linebuf);
		}
		if( (avg.compressibility > 0.0) && std::isfinite(avg.compressibility) ) {
			sprintf( linebuf, "OUTPUT%s: compressibility = %.6g +- %.6g atm^-1\n", sysID, avg.compressibility, avg.compressibility_error);
			Output::out1(linebuf);
			sprintf( linebuf, "OUTPUT%s: bulk modulus = %.6g +- %.6g GPa\n", sysID, ATM2PASCALS*1.0e-9/avg.compressibility,
				ATM2PASCALS*1.0e-9*avg.compressibility_error/avg.compressibility/avg.compressibility);
			Output::out1(linebuf);
		}
	}

	if( (avg.heat_capacity > 0.0)   &&   (std::isfinite(avg.heat_capacity)) ) {
		sprintf( linebuf, "OUTPUT%s: heat capacity = %.5g +- %.5g kJ/mol K\n", sysID, avg.heat_capacity, avg.heat_capacity_error );
		Output::out1(linebuf);
	}

	if( ensemble == ENSEMBLE_NPT  ||  ensemble == ENSEMBLE_REPLAY ) {
		sprintf( linebuf, "OUTPUT%s: volume = %.5f +- %.5f A^3\n", sysID, avg.volume, avg.volume_error );
		Output::out1(linebuf);
	}

	if(avg.spin_ratio > 0.0) {
		sprintf( linebuf, "OUTPUT%s: ortho spin ratio = %.5lf +- %.5lf %%\n", sysID, avg.spin_ratio*100.0, avg.spin_ratio_error*100.0 );
		Output::out1(linebuf);
	}
	
	if( sorbateCount > 1 ){

		for( i=0; i < sorbateCount; i++ ) {

			sprintf( linebuf, "OUTPUT%s: Stats for %s\n", sysID, sorbateInfo[i].id );
			Output::out1(linebuf);
			sprintf( linebuf,
			         "             Average_N(%s)= %.5lf +- %.5lf\n", 
			         sorbateInfo[i].id, sorbateGlobal[i].avgN, 
			         sorbateGlobal[i].avgN_err);
			Output::out1(linebuf);
			sprintf( linebuf,
			         "             Sorbed_Mass(%s)= %.5lf +- %.5lf g/mol\n",
			         sorbateInfo[i].id, sorbateGlobal[i].avgN*sorbateInfo[i].mass, 
			         sorbateGlobal[i].avgN_err*sorbateInfo[i].mass
			);
			Output::out1(linebuf);
			sprintf( linebuf,
			         "             density(%s)= %.5le +- %.5le g/cm^3\n",
			         sorbateInfo[i].id, sorbateGlobal[i].density, 
			         sorbateGlobal[i].density_err
			);
			Output::out1(linebuf);

			if( observables->frozen_mass > 0 ) {
				sprintf( linebuf,
				         "             pore_density(%s)= %.5le +- %.5le g/cm^3\n",
				         sorbateInfo[i].id, sorbateGlobal[i].pore_density, 
				         sorbateGlobal[i].pore_density_err
				);
				sprintf( &linebuf[strlen(linebuf)],
				         "             excess_ratio(%s)= %.5le +- %.5le g/cm^3\n",
				         sorbateInfo[i].id,
				         sorbateGlobal[i].excess_ratio, 
				         sorbateGlobal[i].excess_ratio_err
				);
				sprintf( &linebuf[strlen(linebuf)],
				         "             wt_%%(%s)= %.5lf +- %.5le %%\n",
				         sorbateInfo[i].id,
				         sorbateGlobal[i].percent_wt, 
				         sorbateGlobal[i].percent_wt_err
				);
				sprintf( &linebuf[strlen(linebuf)],
				         "             wt_%%(%s)(ME)= %.5lf +- %.5le %%\n",
				         sorbateInfo[i].id, 
				         sorbateGlobal[i].percent_wt_me, 
				         sorbateGlobal[i].percent_wt_me_err
				);
			}
			sprintf( &linebuf[strlen(linebuf)],
			         "             Selectivity(%s)= %.4lf +- %.4lf\n", 
			         sorbateInfo[i].id,
			         sorbateGlobal[i].selectivity, 
			         sorbateGlobal[i].selectivity_err
			);

			Output::out1(linebuf);
		}
	}
	
	Output::out1("\n");
	return 0;
}



// determine the acceptance rate
void System::track_ar(nodestats_t *ns) {

	if(ns->accept + ns->reject)
		ns->acceptance_rate = ns->accept / ((double)(ns->accept) + ns->reject);
	else
		ns->acceptance_rate = 0;

	if(ns->accept_insert + ns->reject_insert)
		ns->acceptance_rate_insert = ns->accept_insert / ((double)(ns->accept_insert) + ns->reject_insert);
	else
		ns->acceptance_rate_insert = 0;

	if(ns->accept_remove + ns->reject_remove)
		ns->acceptance_rate_remove = ns->accept_remove / ((double)(ns->accept_remove) + ns->reject_remove);
	else
		ns->acceptance_rate_remove = 0;

	if(ns->accept_displace + ns->reject_displace)
		ns->acceptance_rate_displace = ns->accept_displace / ((double)(ns->accept_displace) + ns->reject_displace);
	else
		ns->acceptance_rate_displace = 0;

	if(ns->accept_adiabatic + ns->reject_adiabatic)
		ns->acceptance_rate_adiabatic = ns->accept_adiabatic / ((double)(ns->accept_adiabatic) + ns->reject_adiabatic);
	else
		ns->acceptance_rate_adiabatic = 0;

	if(ns->accept_spinflip + ns->reject_spinflip)
		ns->acceptance_rate_spinflip = ns->accept_spinflip / ((double)(ns->accept_spinflip) + ns->reject_spinflip);
	else
		ns->acceptance_rate_spinflip = 0;

	if(ns->accept_volume + ns->reject_volume)
		ns->acceptance_rate_volume = ns->accept_volume / ((double)(ns->accept_volume) + ns->reject_volume);
	else
		ns->acceptance_rate_volume = 0;

	if(ns->accept_ptemp + ns->reject_ptemp)
		ns->acceptance_rate_ptemp = ns->accept_ptemp / ((double)(ns->accept_ptemp) + ns->reject_ptemp);
	else
		ns->acceptance_rate_ptemp = 0;

	if (ns->accept_beadPerturb + ns->reject_beadPerturb)
		ns->acceptance_rate_beadPerturb = ns->accept_beadPerturb / ((double)(ns->accept_beadPerturb) + ns->reject_beadPerturb);
	else
		ns->acceptance_rate_beadPerturb = 0;
}



// update node statistics related to the processing
void System::update_nodestats( nodestats_t *nstats, avg_nodestats_t *avg_ns ) {

	static int counter = 0;
	double    quantity = 0;
	
	counter++;
		
	double factor   = (counter - 1.0) / counter;   // Weight current average carries in the new/updated average
	double new_fctr =            1.0  / counter;   // Weight new data will carries in the average

	quantity = nstats->boltzmann_factor;
	avg_ns->boltzmann_factor    = factor*avg_ns->boltzmann_factor    + new_fctr * quantity;
	avg_ns->boltzmann_factor_sq = factor*avg_ns->boltzmann_factor_sq + new_fctr * quantity*quantity;

	quantity = nstats->cavity_bias_probability;
	avg_ns->cavity_bias_probability    = factor*avg_ns->cavity_bias_probability    + new_fctr * quantity;
	avg_ns->cavity_bias_probability_sq = factor*avg_ns->cavity_bias_probability_sq + new_fctr * quantity*quantity;

	quantity = nstats->polarization_iterations;
	avg_ns->polarization_iterations    = factor*avg_ns->polarization_iterations    + new_fctr * quantity;
	avg_ns->polarization_iterations_sq = factor*avg_ns->polarization_iterations_sq + new_fctr * quantity*quantity;

	// the remaining items aren't really averages, but cumulative values
	avg_ns->acceptance_rate             = nstats->acceptance_rate;
	avg_ns->acceptance_rate_insert      = nstats->acceptance_rate_insert;
	avg_ns->acceptance_rate_remove      = nstats->acceptance_rate_remove;
	avg_ns->acceptance_rate_displace    = nstats->acceptance_rate_displace;
	avg_ns->acceptance_rate_adiabatic   = nstats->acceptance_rate_adiabatic;
	avg_ns->acceptance_rate_spinflip    = nstats->acceptance_rate_spinflip;
	avg_ns->acceptance_rate_volume      = nstats->acceptance_rate_volume;
	avg_ns->acceptance_rate_ptemp       = nstats->acceptance_rate_ptemp;
	avg_ns->acceptance_rate_beadPerturb = nstats->acceptance_rate_beadPerturb;

}




void System::write_states() {

	Molecule  * molecule_ptr           = nullptr;
	Atom      * atom_ptr               = nullptr;
	FILE      * fp                     = nullptr;
	int         num_frozen_molecules   = 0,
	            num_moveable_molecules = 0,
	            num_frozen_atoms       = 0,
	            num_moveable_atoms     = 0,
	            ext_output             = 0; // By default, PDB compliant coordinates are printed (%8.3f), else extended output is used (%11.6f)
	
	//don't bother if we'd be writing to /dev/null
	if( ! strncmp("/dev/null", traj_output, 9 ))
		return;
 	else
		fp = open_traj_file();

	// count the number of molecules, atoms, etc.
	for( molecule_ptr = molecules;   molecule_ptr;   molecule_ptr = molecule_ptr->next ) {

		if(molecule_ptr->frozen) {
			++num_frozen_molecules;
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				++num_frozen_atoms;
		} else {
			++num_moveable_molecules;
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				++num_moveable_atoms;
		}

	}

	fprintf(fp, "REMARK step=%d\n", step);
	#ifdef _MPI
		if (mpi) 
			fprintf(fp, "REMARK node=%d\n", (int) rank);
		if (parallel_tempering)
			fprintf(fp, "REMARK temperature=%.6lf\n", (double) ptemp->templist[ptemp->index[rank]]);
	#endif

	fprintf(fp, "REMARK total_molecules=%d, total_atoms=%d\n", 
		(num_frozen_molecules + num_moveable_molecules), (num_frozen_atoms + num_moveable_atoms));
	fprintf(fp, "REMARK frozen_molecules=%d, moveable_molecules=%d\n", 
		num_frozen_molecules, num_moveable_molecules);
	fprintf(fp, "REMARK frozen_atoms=%d, moveable_atoms=%d\n", 
		num_frozen_atoms, num_moveable_atoms);

	// Check if extended coordinate output is needed (CRC)
	if( long_output )
		ext_output = 1;
	else if( (pbc.basis[0][0] >= 200.0) || (pbc.basis[0][1] >= 200.0) || (pbc.basis[0][2] >= 200.0) || 
	         (pbc.basis[1][0] >= 200.0) || (pbc.basis[1][1] >= 200.0) || (pbc.basis[1][2] >= 200.0) || 
	         (pbc.basis[2][0] >= 200.0) || (pbc.basis[2][1] >= 200.0) || (pbc.basis[2][2] >= 200.0) )
		ext_output = 1;

	// write PBC data 
	fprintf(fp,"CRYST1");
	fprintf(fp,"%9.3f",sqrt(UsefulMath::dddotprod(pbc.basis[0], pbc.basis[0])));
	fprintf(fp,"%9.3f",sqrt(UsefulMath::dddotprod(pbc.basis[1], pbc.basis[1])));
	fprintf(fp,"%9.3f",sqrt(UsefulMath::dddotprod(pbc.basis[2], pbc.basis[2])));
	fprintf(fp,"%7.2f", 180.0/pi*acos( UsefulMath::dddotprod(pbc.basis[1],pbc.basis[2]) / sqrt( UsefulMath::dddotprod(pbc.basis[1], pbc.basis[1]) * UsefulMath::dddotprod(pbc.basis[2], pbc.basis[2]) ) ));
	fprintf(fp,"%7.2f", 180.0/pi*acos( UsefulMath::dddotprod(pbc.basis[2],pbc.basis[0]) / sqrt( UsefulMath::dddotprod(pbc.basis[0], pbc.basis[0]) * UsefulMath::dddotprod(pbc.basis[2], pbc.basis[2]) ) ));
	fprintf(fp,"%7.2f", 180.0/pi*acos( UsefulMath::dddotprod(pbc.basis[0],pbc.basis[1]) / sqrt( UsefulMath::dddotprod(pbc.basis[1], pbc.basis[1]) * UsefulMath::dddotprod(pbc.basis[0], pbc.basis[0]) ) ));
	fprintf(fp,"\n");

	// write pqr formatted states
	int i=1, j=1;
	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next, j++) {
		for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			fprintf(fp, "ATOM  ");
			fprintf(fp, "%5d", i);		// give each one a unique id
			fprintf(fp, " %-4.45s", atom_ptr->atomtype);
			fprintf(fp, " %-3.3s ", molecule_ptr->moleculetype);
			if(atom_ptr->adiabatic)
				fprintf(fp, "%-1.1s", "A");
			else if(atom_ptr->frozen)
				fprintf(fp, "%-1.1s", "F");
			else if(atom_ptr->spectre)
				fprintf(fp, "%-1.1s", "S");
			else if(atom_ptr->target)
				fprintf(fp, "%-1.1s", "T");
			else
				fprintf(fp, "%-1.1s", "M");
			fprintf(fp, "%4d    ", j);		// give each molecule a unique id 

			if(ext_output == 0) {
				// Regular (PDB compliant) Coordinate Output 
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[0]);
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[1]);
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[2]);
			} else {
				// Extended (PQR) Coordinate Output 
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[0]);
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[1]);
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[2]);
			}

			fprintf(fp, " %8.4f", atom_ptr->mass);
			fprintf(fp, " %8.4f", atom_ptr->charge/E2REDUCED);	// convert charge back to real units
			fprintf(fp, " %8.5f", atom_ptr->polarizability);
			fprintf(fp, " %8.5f", atom_ptr->epsilon);
			fprintf(fp, " %8.5f", atom_ptr->sigma);
			fprintf(fp, " %8.5f", atom_ptr->omega);
			fprintf(fp, " %8.5f", atom_ptr->gwp_alpha);
			fprintf(fp, " %8.5f", atom_ptr->c6);
			fprintf(fp, " %8.5f", atom_ptr->c8);
			fprintf(fp, " %8.5f", atom_ptr->c10);
			fprintf(fp, " %8.5f", atom_ptr->c9);
			fprintf(fp, "\n");

		}
	}

	//write basis to the output file. needed for restarting NPT jobs, or whatever.
	fprintf(fp, "REMARK BOX BASIS[0] = %20.14lf %20.14lf %20.14lf\n", 
		pbc.basis[0][0], pbc.basis[0][1], pbc.basis[0][2]);
	fprintf(fp, "REMARK BOX BASIS[1] = %20.14lf %20.14lf %20.14lf\n",
		pbc.basis[1][0], pbc.basis[1][1], pbc.basis[1][2]);
	fprintf(fp, "REMARK BOX BASIS[2] = %20.14lf %20.14lf %20.14lf\n", 
		pbc.basis[2][0], pbc.basis[2][1], pbc.basis[2][2]);

	fprintf(fp, "ENDMDL\n");
	fflush(fp);
	fclose(fp);

}




FILE * System::open_traj_file() {
	FILE * fp;
	char * filename;
	static int clobber = 1; //if clobber is set, we will overwrite old files

		//open files for append
	if( traj_output ) {

		if (mpi) {
			#ifdef _MPI // each node will write it's own file
				if (parallel_tempering)
					filename = Output::make_filename(traj_output, ptemp->index[rank]); //append bath index to filename
				else
					filename = Output::make_filename(traj_output, (int) rank); //append core index to filename
			#endif

		} else {

			filename = traj_output;
		}



		if ( clobber == 1 ) {
			fp = SafeOps::openFile(filename, "w", __LINE__, __FILE__ );
			clobber = 0; //don't clobber again
		}
		else {
			fp = SafeOps::openFile(filename, "a", __LINE__, __FILE__ );
		}

		
		#ifdef _MPI
		if (mpi) { SafeOps::free(filename); }
		#endif

		return fp;
	}

	return nullptr;
}




int System::write_molecules_wrapper( char * filename ) {
	int rval = -1;
	char  filenameold[maxLine]; 
	FILE * fp;

	if (mpi) {

		#ifdef _MPI

			int j = 0;
			char* filenameno = nullptr;

			//make a new filename with the core/or bath number appended
			if (parallel_tempering)
				filenameno = Output::make_filename(filename, ptemp->index[rank]); //append bath index to filename
			else
				filenameno = Output::make_filename(filename, (int) rank); //append core index to filename

			//move most recent state file to file.last
			if (SafeOps::file_exists(filename)) {
				sprintf(filenameold, "%s.last", filenameno);
				int code = rename(filenameno, filenameold); // fxn marked _Check_return_ in VS
				if (code) Output::err("WARNING: Unable to rename .last file.");
			}

			// open the file and free the filename string
			fp = SafeOps::openFile(filenameno, "w", __LINE__, __FILE__);
			SafeOps::free(filenameno);

			// we write files one at a time to avoid disk congestion
			for (j = 0; j < size; j++) {
				MPI_Barrier(MPI_COMM_WORLD);
				if (j == rank)
					rval = write_molecules(fp);
			}

			//free the file pointer
			fclose(fp);
			
		#endif // _MPI

	} else {

		//move most recent state "filename" to "filename.last"
		if (SafeOps::file_exists(filename)) {
			sprintf(filenameold, "%s.last", filename);
			int code = rename(filename, filenameold); // fxn marked _Check_return_ in VS
			if (code) Output::err("WARNING: Unable to rename .last file.");
		}

		fp = SafeOps::openFile(filename, "w", __LINE__, __FILE__);

		//write the file
		rval = write_molecules(fp);

		fclose(fp);
	}

	return rval;
}



// write out the final system state as a PQR file 
int System::write_molecules(FILE * fp) {

	int    atom_box            =  0,
	       molecule_box        =  0,
	       p                   =  0,
	       q                   =  0;
	double box_pos[3]          = {0},
	       box_occupancy[3]    = {0};
	int    l                   =  0,
	       m                   =  0,
	       n                   =  0,
	       box_labels[2][2][2] = {0},
	       diff                =  0;
	
	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	int        i            = 0,
	           j            = 0,
	           k            = 0;
	int        ext_output   = 0;
	
	// Check if extended coordinate output is needed (CRC) 
	// By default, PDB compliant coordinates are printed (%8.3f), else extended output is used (%11.6f)
	if(long_output)
		ext_output = 1;
	else if( (pbc.basis[0][0] >= 100.0) || (pbc.basis[0][1] >= 100.0) || (pbc.basis[0][2] >= 100.0) || 
	         (pbc.basis[1][0] >= 100.0) || (pbc.basis[1][1] >= 100.0) || (pbc.basis[1][2] >= 100.0) || 
	         (pbc.basis[2][0] >= 100.0) || (pbc.basis[2][1] >= 100.0) || (pbc.basis[2][2] >= 100.0) )
		ext_output = 1;
	else
		ext_output = 0;

	// write PBC data 
	// VMD uses a weird convention which essentially reverses alpha <-> beta
	fprintf(fp,"CRYST1");
	fprintf(fp,"%9.3f",sqrt(UsefulMath::dddotprod(pbc.basis[0], pbc.basis[0])));
	fprintf(fp,"%9.3f",sqrt(UsefulMath::dddotprod(pbc.basis[1], pbc.basis[1])));
	fprintf(fp,"%9.3f",sqrt(UsefulMath::dddotprod(pbc.basis[2], pbc.basis[2])));
	fprintf(fp,"%7.2f", 180.0/pi*acos( UsefulMath::dddotprod(pbc.basis[2],pbc.basis[0]) / sqrt( UsefulMath::dddotprod(pbc.basis[0], pbc.basis[0]) * UsefulMath::dddotprod(pbc.basis[2], pbc.basis[2]) ) ));
	fprintf(fp,"%7.2f", 180.0/pi*acos( UsefulMath::dddotprod(pbc.basis[1],pbc.basis[2]) / sqrt( UsefulMath::dddotprod(pbc.basis[1], pbc.basis[1]) * UsefulMath::dddotprod(pbc.basis[2], pbc.basis[2]) ) ));
	fprintf(fp,"%7.2f", 180.0/pi*acos( UsefulMath::dddotprod(pbc.basis[0],pbc.basis[1]) / sqrt( UsefulMath::dddotprod(pbc.basis[1], pbc.basis[1]) * UsefulMath::dddotprod(pbc.basis[0], pbc.basis[0]) ) ));
	fprintf(fp,"\n");
	

	// write pqr
	for(molecule_ptr = molecules, i = 1, j = 1; molecule_ptr; molecule_ptr = molecule_ptr->next, j++) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			fprintf(fp, "ATOM  ");
			fprintf(fp, "%5d", i);		/* give each one a unique id */
			fprintf(fp, " %-4.45s", atom_ptr->atomtype);
			fprintf(fp, " %-3.3s ", molecule_ptr->moleculetype);
			if(atom_ptr->adiabatic)
				fprintf(fp, "%-1.1s", "A");
			else if(atom_ptr->frozen)
				fprintf(fp, "%-1.1s", "F");
			else if(atom_ptr->spectre)
				fprintf(fp, "%-1.1s", "S");
			else if(atom_ptr->target)
				fprintf(fp, "%-1.1s", "T");
			else
				fprintf(fp, "%-1.1s", "M");
			if(independent_particle)
				fprintf(fp, " %4d   ", i);
			else
				fprintf(fp, " %4d   ", j);		// give each molecule a unique id 

			// Regular (PDB compliant) Coordinate Output 
			if( (wrapall) && (ext_output == 0) ) {
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[0]);
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[1]);
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[2]);
			} else if(ext_output == 0){
				fprintf(fp, "%8.3f", atom_ptr->pos[0]);
				fprintf(fp, "%8.3f", atom_ptr->pos[1]);
				fprintf(fp, "%8.3f", atom_ptr->pos[2]);
			}

			// Extended (PQR) Coordinate Output
			if( (wrapall) && (ext_output == 1) ) {
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[0]);
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[1]);
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[2]);
			} else if (ext_output == 1) {
				fprintf(fp, "%11.6f ", atom_ptr->pos[0]);
				fprintf(fp, "%11.6f ", atom_ptr->pos[1]);
				fprintf(fp, "%11.6f ", atom_ptr->pos[2]);
			}
			fprintf(fp, " %8.5f", atom_ptr->mass);
			fprintf(fp, " %8.5f", atom_ptr->charge/E2REDUCED);	// convert charge back to real units
			fprintf(fp, " %8.5f", atom_ptr->polarizability);
			fprintf(fp, " %8.5f", atom_ptr->epsilon);
			fprintf(fp, " %8.5f", atom_ptr->sigma);
			fprintf(fp, " %8.5f", atom_ptr->omega);
			fprintf(fp, " %8.5f", atom_ptr->gwp_alpha);
			fprintf(fp, " %8.5f", atom_ptr->c6);
			fprintf(fp, " %8.5f", atom_ptr->c8);
			fprintf(fp, " %8.5f", atom_ptr->c10);
			fprintf(fp, " %8.5f", atom_ptr->c9);
			fprintf(fp, "\n");

		}
	}

	if(wrapall) {

		// output the box coords as virtual particles for visualization
		atom_box = i;
		molecule_box = j;
		for(i = 0; i < 2; i++) {
			for(j = 0; j < 2; j++) {
				for(k = 0; k < 2; k++) {

					// make this frozen
					fprintf(fp, "ATOM  ");
					fprintf(fp, "%5d", atom_box);
					fprintf(fp, " %-4.45s", "X");
					fprintf(fp, " %-3.3s ", "BOX");
					fprintf(fp, "%-1.1s", "F");
					fprintf(fp, " %4d   ", molecule_box);

					// box coords
					box_occupancy[0] = ((double)i) - 0.5;
					box_occupancy[1] = ((double)j) - 0.5;
					box_occupancy[2] = ((double)k) - 0.5;

					for(p = 0; p < 3; p++)
						for(q = 0, box_pos[p] = 0; q < 3; q++)
							box_pos[p] += pbc.basis[q][p]*box_occupancy[q];

					for(p = 0; p < 3; p++)
						if(ext_output == 0)
							fprintf(fp, "%8.3f", box_pos[p]);
						else
							fprintf(fp, "%11.6f ", box_pos[p]);

					// null interactions
					fprintf(fp, " %8.4f", 0.0);
					fprintf(fp, " %8.4f", 0.0);
					fprintf(fp, " %8.5f", 0.0);
					fprintf(fp, " %8.5f", 0.0);
					fprintf(fp, " %8.5f", 0.0);
					fprintf(fp, "\n");

					box_labels[i][j][k] = atom_box;
					++atom_box;

				}
			}
		}

		for(i = 0; i < 2; i++) 
			for(j = 0; j < 2; j++) 
				for(k = 0; k < 2; k++) 
					for(l = 0; l < 2; l++) 
						for (m = 0; m < 2; m++) 
							for (n = 0; n < 2; n++) {

								diff = (int)fabs(i - l) + (int)fabs(j - m) + (int)fabs(k - n);
								if (diff == 1)
									fprintf(fp, "CONECT %4d %4d\n", box_labels[i][j][k], box_labels[l][m][n]);
							}

		

	} // if wrapall 

	// write basis to the output file. needed for restarting NPT jobs, or whatever.
	fprintf(fp, "REMARK BOX BASIS[0] = %20.14lf %20.14lf %20.14lf\n", 
		pbc.basis[0][0], pbc.basis[0][1], pbc.basis[0][2]);
	fprintf(fp, "REMARK BOX BASIS[1] = %20.14lf %20.14lf %20.14lf\n",
		pbc.basis[1][0], pbc.basis[1][1], pbc.basis[1][2]);
	fprintf(fp, "REMARK BOX BASIS[2] = %20.14lf %20.14lf %20.14lf\n", 
		pbc.basis[2][0], pbc.basis[2][1], pbc.basis[2][2]);

	// if surface fitting, write some surface fit info as a remark
	if( ensemble == ENSEMBLE_SURF_FIT ) {
		if( fit_boltzmann_weight )
			fprintf(fp,"REMARK SURFACE FITTING BOLTZMANN TEMPERATURE = %lf\n", fit_max_energy);
		else
			fprintf(fp,"REMARK SURFACE FITTING MAX_ENERGY = %lf\n", fit_max_energy);
		fprintf(fp,"REMARK SURFACE FITTING SQUARE_ERROR = %lf\n", fit_best_square_error);
	}    

	// output the connectivity information 
	fprintf(fp, "END\n");
	fflush(fp);

	return(0);

}




FILE * System::open_dipole_file() {
	FILE     * fp       = nullptr;
	char     * filename = nullptr;
	static int clobber  = 1; //if clobber is set, we will overwrite old files

	//open files for append
	if( dipole_output ) {
		#ifdef _MPI // each  will write it's own file
			if( parallel_tempering )
				filename = Output::make_filename( dipole_output, ptemp->index[rank]); //append bath index to filename
			else 
				filename = Output::make_filename( dipole_output, (int) rank); //append core index to filename
		#else
			filename = dipole_output;
		#endif // MPI

		if( clobber == 1 ) {
			fp = SafeOps::openFile( filename, "w", __LINE__, __FILE__ );
			clobber = 0; //don't clobber again
		}
		else {
			fp = SafeOps::openFile( filename, "a", __LINE__, __FILE__ );
		}

		#ifdef _MPI
			SafeOps::free(filename);
		#endif
		return fp;
	}

	return nullptr;
}



// output each molecular dipole (in debye) per line 
void System::write_dipole() {

	FILE     * fp           = nullptr;
	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	double     dipole[3]    = {0};

	//don't bother if we'd be writing to /dev/null
	if( ! strncmp( "/dev/null", dipole_output, 9 ))
		return;
 	else
		fp = open_dipole_file();

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for( int p = 0; p < 3; p++ )
			dipole[p] = 0;
		for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for( int p = 0; p < 3; p++ )
				dipole[p] += atom_ptr->mu[p];
		}
		if(!molecule_ptr->frozen) fprintf(fp, "%f %f %f\n", dipole[0]/DEBYE2SKA, dipole[1]/DEBYE2SKA, dipole[2]/DEBYE2SKA);
	}

	fflush(fp);
	fclose(fp);

	return;
}




FILE * System::open_field_file() {
	FILE     * fp       = nullptr;
	char     * filename = nullptr;
	static int clobber  = 1; //if clobber is set, we will overwrite old files

	//open files for append
	if( field_output ) {
		#ifdef _MPI // each node will write it's own file
			if ( parallel_tempering )
				filename = Output::make_filename( field_output, ptemp->index[rank] ); //append bath index to filename
			else 
				filename = Output::make_filename( field_output, (int) rank ); //append core index to filename
		#else
				filename = field_output;
		#endif // MPI

		if ( clobber == 1 ) {
			fp = SafeOps::openFile( filename, "w", __LINE__, __FILE__ );
			clobber = 0; //don't clobber again
		}
		else {
			fp = SafeOps::openFile( filename, "a", __LINE__, __FILE__ );
		}

		#ifdef _MPI
			SafeOps::free( filename );
		#endif
		return fp;
	}

	return nullptr;
}



// output the total molecular electrostatic field (in e/A) per line)
void System::write_field() {

	FILE * fp;
	Molecule *molecule_ptr;
	
	Atom *atom_ptr;
	double field[3] = {0};

	//don't bother if we'd be writing to /dev/null
	if ( ! strncmp("/dev/null", field_output, 9 ))
		return;
 	else
		fp = open_field_file();

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for( int p = 0; p < 3; p++ )
			field[p] = 0;
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for( int p = 0; p < 3; p++ )
				field[p] += atom_ptr->ef_static[p] + atom_ptr->ef_induced[p];
		}
		if( !molecule_ptr->frozen )
			fprintf(fp, "%f %f %f\n", field[0]/E2REDUCED, field[1]/E2REDUCED, field[2]/E2REDUCED);
	}

	fflush(fp);
	fclose(fp);

	return;
}




int System::write_performance(  unsigned int i ) {

	static struct timeval current_time, last_time;

	char       linebuf[maxLine] = {'\0'};
	double     sec_step         = 0;
	static int last_step        = 0;

	Output::GetTimeOfDay( &current_time );
	

	if( i > corrtime ) {
		double di = (double) i;
		sec_step = Output::calctimediff(current_time, last_time) / (di - last_step);

		if( ensemble == ENSEMBLE_UVT ) {
			sprintf(linebuf, "OUTPUT: Grand Canonical Monte Carlo simulation running on %d core(s)\n", (int)(mpi?size:1));
			Output::out1( linebuf );
		} else {
			sprintf( linebuf, "OUTPUT: Canonical Monte Carlo simulation running on %d core(s)\n", (int)(mpi?size:1) );
			Output::out1( linebuf );
		}


		#if defined(WIN32) || defined(_WIN32) || defined(__WIN32)
			time_t time = current_time.tv_sec;
			sprintf(linebuf, "OUTPUT: Root collecting statistics at %s", ctime(&time));
		#else
			sprintf(linebuf, "OUTPUT: Root collecting statistics at %s", ctime( &(current_time.tv_sec)) );
		#endif
		
		Output::out1( linebuf );
		sprintf(linebuf, "OUTPUT: Completed step %d/%d  (%.3f %%)\n", i, numsteps, (di/numsteps)*100);
		Output::out1( linebuf );
		sprintf(linebuf, "OUTPUT: %.3lf sec/step, ETA = %.3lf hrs\n", sec_step, sec_step*(numsteps - di)/3600.0);
		Output::out1( linebuf );

	}	

	last_step = i;
	last_time.tv_sec  = current_time.tv_sec;
	last_time.tv_usec = current_time.tv_usec;

	return 0;

}