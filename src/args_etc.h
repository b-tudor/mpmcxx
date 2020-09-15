#pragma once
#ifndef ARGS_ETC_H
#define	ARGS_ETC_H

#include <cstring>
#include <iostream>
#include <sstream>
#include <stdlib.h>

#ifdef _OMP
#	include <omp.h>
#endif




#include "Output.h"
#include "SafeOps.h"
#include "SimulationControl.h"

using uidx = size_t;
extern uidx rank;
extern uidx size;
extern bool mpi;






typedef struct _parameters {
	char* prog_name;
	char* in_filename;
	uidx  Ptrotter_number; // P AKA trotter number (this variable is prounounced "Trotter" -- think pterodactyl).
	bool  write_PI_Frames_at_corrtime;
	char* PI_frame_file; // file name of PI frames file
} params;







//\//// Clean up & exit /////////////////////////////////////////////////
void die(int success_code) {

	#ifdef _MPI
		// Kill MPI before quitting, when neccessary
		if (mpi) { MPI_Finalize(); }
	#endif

	if (success_code)
		exit(EXIT_SUCCESS);
	else
		exit(EXIT_FAILURE);
} //\////////////////////////////////////////////////////////////////////




//\/// Install the signal handler for clean exits (on Unix systems) /////
#ifdef __unix__
#	include <unistd.h>
#	include <signal.h>
#endif 
SimulationControl* sc;
void install_signal_handler(SimulationControl* simControl) {
	sc = simControl;
#if defined( _POSIX_VERSION )  
	signal(SIGTERM, signal_handler);
	signal(SIGUSR1, signal_handler);
	signal(SIGUSR2, signal_handler);
#endif
}
//\//////////////////////////////////////////////////////////////////////




// Print brief application  user instructions and exit
void displayUsageAndDie( params &p ) {

	char linebuf[maxLine];

	if (!rank) {
		sprintf( linebuf, "MPMC++\nMassively Parallel Monte Carlo: Multi-System Edition, v%s -- 2012-2019 GNU Public License\n", VERSION );
		Output::out(linebuf);
		std::cout << "\nUsage:" << std::endl;
		std::cout << "\t" << p.prog_name << " INPUT_FILE [options]" << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "\t-P X        Where X is the trotter number for a non-MPI path integral job." << std::endl;
		std::cout << "\t-xyz  FILE  Write .xyz PI visualization frames at corrtime steps to file FILE." << std::endl << std::endl;
		std::cout << "\nWhen using MPI with path integral ensembles, the number of MPI processes will determine the" << std::endl;
		std::cout << "number of path integral \"beads\" (i.e. the Trotter number). It is therefore the uneccessary to " << std::endl;
		std::cout << "use the -P option for multi-threaded MPI runs. For multi-site/multi-atom sorbate molecule, the" << std::endl;
		std::cout << "Trotter number must be a power of 2 greater than or equal to 4." << std::endl;
		
		std::cout << "\nExample:\n\t" << p.prog_name << " -P 8 <my_path_integral_sim_file>" << std::endl;
		std::cout << "\n        \tmpirun -np 8 " << p.prog_name << " <my_path_integral_sim_file>" << std::endl;
		std::cout << "\nSee https://github.com/mpmccode/mpmc for input file specification." << std::endl;
		#ifndef _MPI
			std::cout << "\nNOTE: MPI Functionality is not compiled into this executable." << std::endl;
		#endif
	}
	
	if (p.in_filename) { SafeOps::free(p.in_filename); }

	die(EXIT_FAILURE);
}




//\/// Display program name and welcome message. ////////////////////////
void introduce_self() {

	char linebuf[maxLine];
	time_t t = time(nullptr);
	struct tm tm = *localtime(&t);

	char mpi_msg[20] = { 0 };
	#ifdef _MPI
		sprintf(mpi_msg, " (MPI enabled)");
	#endif

	sprintf(
		linebuf,
		"MPMC++\nMassively Parallel Monte Carlo: Multi-System Edition%s, v%s -- 2012-2019 GNU Public License\n",
		mpi_msg, 
		VERSION
	);
	Output::out1(linebuf);

	if (mpi) {
		sprintf(
			linebuf,
			"MAIN: process%s started on %d MPI thread%s @ %d-%d-%d %d:%d:%d\n",
			((size>1) ? "es" : ""),
			(int)size,
			((size>1) ? "s" : ""),
			tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec
		);
	}
	else {

		uidx thread_count = 1;
				
		#pragma omp parallel 
		{
			int np = 1;
			int id = 0;
			#ifdef _OMP
				np = omp_get_num_threads();
				id = omp_get_thread_num();
			#endif
			if (!id)
				thread_count = (uidx) np;
		}
		sprintf(
			linebuf,
			"MAIN: process%s started on %d thread%s @ %d-%d-%d %d:%d:%d\n",
			((thread_count > 1) ? "es" : ""),
			(int)thread_count,
			((thread_count > 1) ? "s" : ""),
			tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec
		);
	}
	Output::out1(linebuf);
}



//\/// Check if MPI or openMP service is available at runtime and configure job accordingly
void parallel_introspection_and_initialization(int& argc, char* argv[], uidx P) {

	#ifdef _MPI	 // Start up the MPI chain
		if (MPI_Init(&argc, &argv) == MPI_SUCCESS)
		{
			int r, s;
			MPI_Comm_rank(MPI_COMM_WORLD, &r);
			MPI_Comm_size(MPI_COMM_WORLD, &s);
			rank = (uidx) r;
			size = (uidx) s;
		}
	#endif

	mpi = size > 1; // signals a multi-threaded run
	if( mpi ) {
		if (P) {
			// If we get here, we are using MPI, but a Trotter number has been manually specified 
			char linebuf[maxLine];
			sprintf(linebuf, "MPMC++\nMassively Parallel Monte Carlo: Multi-System Edition, v%s -- 2012-2019 GNU Public License\n", VERSION);
			Output::out(linebuf);
			Output::err("Do not explicitly set a Trotter number (-P) for MPI jobs. The number of MPI threads\n");
			Output::err("(i.e. MPI 'size') implicitly sets the Trotter number. E.g. for P=8, try:\n");
			Output::err("\tmpirun -np 8 mpmc++ my-input-file\n");
			die(fail);
		}
	} 

	#ifdef _MPI
		// finalize MPI services if MPI is available, but we are only running on a single MPI thread
		else { MPI_Finalize(); }
	#endif

	#ifdef _OMP
		if (!mpi) 
			omp_set_num_threads((int)P);
	#endif
}



//\/// Removes the path from a filename, leaving on the filename ////////
char* stripPath(char* full_path) {

	// Select either MS Windows or Unix file separator
	#ifndef _WIN32
		const char file_separator = '/';
	#else
		const char file_separator = '\\';
	#endif

	char* simple_filename = full_path;
	char* path = full_path;

	// Skip to the last separator in the filename
	while (*path != '\0') {
		path++;
		if ((*path) == file_separator)
			simple_filename = path + 1;
	}

	return simple_filename;
}



//\/// Parse command line arguments /////////////////////////////////////
void processArgs(int argc, char* argv[], params &p) {
	
	std::string fname;
	char linebuf[maxLine];
	int n = 1;         // index into argv for argument we are currently parsing

	// DEFAULT VALUES
	p.prog_name                   = stripPath(argv[0]);  // filename of the executable on host system
	p.in_filename                 = nullptr;
	p.Ptrotter_number             = 0;
	p.write_PI_Frames_at_corrtime = false;
	p.PI_frame_file               = nullptr;



	if (argc >= 2) {

		while (n < argc) {

			std::istringstream issOptionToken(argv[n]);

			// report acceptance/rejection 
			if (!strncmp(issOptionToken.str().c_str(), "-P", 3)) {
				n++;
				int trotter;
				if (SafeOps::atoi(argv[n], trotter )) {
					p.Ptrotter_number = (uidx) trotter;
					n++;
					continue;
				} 
			}

			// write xyz frames for visualization (at corr time)
			else if (!strncmp(issOptionToken.str().c_str(), "-xyz", 5)) {
				n++;
				p.write_PI_Frames_at_corrtime = true;
				std::istringstream issArgToken(argv[n]);
				issArgToken >> fname;
				if (!issArgToken)
					displayUsageAndDie(p);
				p.PI_frame_file = (char*)calloc(fname.size() + 1, sizeof(char));
				if (p.PI_frame_file)
					std::memcpy(p.PI_frame_file, fname.c_str(), fname.size());
				else {
					Output::err("Unable to allocate xyz PI filename buffer.\n");
					die(fail);
				}
				n++;
				continue;
			}

			else // if token is not a flag (-x) then it must be the filename
			{
				// INPUT FILENAME
				if (p.in_filename) {
					// if filename is already set, then the command is malformed. There should
					// only be one argument that is not preceeded by an option flag (e.g. "-P")
					sprintf(linebuf, "ERROR: Multiple input files specified(or unrecognized option \"%s\")\n", argv[n]);
					Output::err(linebuf);
					displayUsageAndDie(p);
				}

				issOptionToken >> fname;
				if (!issOptionToken)
					displayUsageAndDie(p);
				p.in_filename = (char*)calloc(fname.size() + 1, sizeof(char));
				if (p.in_filename) 
					std::memcpy(p.in_filename, fname.c_str(), fname.size());
				else {
					Output::err("Unable to allocate input filename buffer.\n");
					die(fail);
				}
				n++;
			}

		} // end arg processing
	} // end arg count

	else displayUsageAndDie(p);
}


  /////////////////////////////////////////////////////////////
 ////    ADD SUPPORT FOR CLEAN UP IN signal_handler()     ////
/////////////////////////////////////////////////////////////



// on SIGTERM, cleanup and exit 
void signal_handler(int sigtype) {
	char linebuf[maxLine];




	#ifdef __unix__
		if (sigtype == SIGTERM)
			Output::out(" ************ SIGTERM received, exiting *************\n");
		else if (sigtype == SIGUSR1)
			Output::out(" ************ SIGUSR1 received, exiting *************\n");
		else if (sigtype == SIGUSR2)
			Output::out(" ************ SIGUSR2 received, exiting *************\n");
		else
			Output::out(" ************ Unknown interrupt signal received, exiting *************\n");

		//sc->close_files();
		//sc->cleanup();
	#ifdef _MPI
	//	if(!rank)
	//		sc->close_files();
	#else
		die(EXIT_FAILURE);
	#endif // MPI
	#else
		sprintf(linebuf, " ************ Interrupt signal ( %d ) received, exiting *************\n", sigtype);
		Output::out(linebuf);
	#endif	
	return;
}



#endif	/* ARGS_H */