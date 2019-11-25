#pragma once
#ifndef ARGS_ETC_H
#define	ARGS_ETC_H

#include <cstring>
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include "Output.h"
#include "SafeOps.h"
#include "SimulationControl.h"

#ifdef __unix__
#include <unistd.h>
#include <signal.h>
#endif



extern int rank;
extern int size;
extern bool mpi;



// Includes for interrupt signal handler 
SimulationControl* sc;

typedef struct _parameters {
	char* in_filename;
	int   Ptrotter_number; // P AKA trotter number (this variable is prounounced "Trotter" -- think pterodactyl).
	bool  write_PI_Frames_at_corrtime;
} params;



void  die(int code);                                                 // Kill any MPI Processes and stop execution
void  displayUsageAndDie(char* progname, params* p);                 // Print brief application  user instructions and exit
void  install_signal_handler(SimulationControl *sc);                 // Initialize signal handlers for clean exits (on Posix systems)
void  introduce_self();                                              // Display program name and welcome message. 
void  mpi_introspection_and_initialization(int& argc, char* argv[]); // Check if MPI service is available at runtime and configure job accordingly
void  processArgs(int argc, char* argv[], params* p);                // Parse command line arguments
void  signal_handler(int sigtype);                                   // Function to handle interrupt signals (on Posix systems)
char* stripPath(char* full_path);                                    // Removes the path from a filename, leaving on the filename





// Kill MPI before quitting, when neccessary
void die(int success_code) {

	#ifdef _MPI
		if (mpi) { MPI_Finalize(); }
	#endif

	if (success_code)
		exit(EXIT_SUCCESS);
	else
		exit(EXIT_FAILURE);

}



// Print brief application  user instructions and exit
void displayUsageAndDie(char* progname_w_path, params &p) {

	if (!rank) {
		char linebuf[maxLine];
		char* progName = stripPath(progname_w_path);
		sprintf( linebuf, "MPMC++\nMassively Parallel Monte Carlo: Multi-System Edition, v%s -- 2012-2019 GNU Public License\n", VERSION );
		Output::out(linebuf);
		std::cout << "\nUsage:" << std::endl;
		std::cout << "\t" << progName << " INPUT_FILE [options]" << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "\t-P X        Where X is the trotter number for a single threaded path integral job." << std::endl;
		std::cout << "\t-xyz        Write .xyz PI visualization frames at corrtime steps." << std::endl << std::endl;
		std::cout << "\nWhen using MPI with path integral ensembles, the number of MPI processes will determine the" << std::endl;
		std::cout << "number of path integral \"beads\" (i.e. the Trotter number). It is therefore the uneccessary to " << std::endl;
		std::cout << "use the -P option for multi-threaded/MPI runs. For multi-site/multi-atom sorbate molecule, the" << std::endl;
		std::cout << "Trotter number must be a power of 2 greater than or equal to 4." << std::endl;
		
		std::cout << "\nExample:\n\t" << progName << " -n 8 <my_path_integral_sim_file>" << std::endl;
		std::cout << "\n        \tmpirun -n 8 " << progName << " <my_path_integral_sim_file>" << std::endl;
		std::cout << "\nSee https://github.com/mpmccode/mpmc for input file specification." << std::endl;
		#ifndef _MPI
			std::cout << "\nNOTE: MPI Functionality is not compiled into this executable." << std::endl;
		#endif
	}
	
	if (p.in_filename) { SafeOps::free(p.in_filename); }

	die(EXIT_FAILURE);
}



// Install the signal handler for clean exits (on Unix systems)
void install_signal_handler(SimulationControl* simControl) {
	sc = simControl;
	#if defined( _POSIX_VERSION )  
		signal(SIGTERM, signal_handler);
		signal(SIGUSR1, signal_handler);
		signal(SIGUSR2, signal_handler);
	#endif
}



// Display program name and welcome message. 
void introduce_self() {

	char linebuf[maxLine];
	time_t t = time(nullptr);
	struct tm tm = *localtime(&t);

	char mpi_msg[20] = { 0 };
	#ifdef _MPI
		sprintf(mpi_msg, " (MPI eneabled)");
	#endif

	sprintf(
		linebuf,
		"MPMC++\nMassively Parallel Monte Carlo: Multi-System Edition%s, v%s -- 2012-2019 GNU Public License\n",
		mpi_msg, 
		VERSION
	);
	Output::out1(linebuf);

	sprintf(
		linebuf,
		"MAIN: process%s started on %d thread%s @ %d-%d-%d %d:%d:%d\n",
		(mpi ? "es" : ""),
		(mpi ? size : 1 ),
		(mpi ?  "s" : ""),
		tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec
	);
	Output::out1(linebuf);
}



// Check if MPI service is available at runtime and configure job accordingly
void mpi_introspection_and_initialization(int& argc, char* argv[], int P) {

	#ifdef _MPI	 // Start up the MPI chain
		if (MPI_Init(&argc, &argv) == MPI_SUCCESS)
		{
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			MPI_Comm_size(MPI_COMM_WORLD, &size);
		}
	#endif

	mpi = size > 1; // signals a multi-threaded run
	if( mpi ) {
		if (P) {
			// If we get here, we are using MPI, but a Trotter number has been manually specified 
			Output::err("Do not explicitly set a Trotter number (-P) for MPI jobs. The number of MPI threads\n");
			Output::err("(i.e. MPI 'size') implicitly sets the Trotter number. E.g. for P=8, try:\n");
			Output::err("\tmpirun -P 8 mpmc++ my-input-file\n");
			die(fail);
		}
	} else {
		MPI_Finalize();
	}
}



// Parse command line arguments
void processArgs(int argc, char* argv[], params &p) {

	std::string fname;
	char linebuf[maxLine];
	int n = 1;         // index into argv for argument we are currently parsing

	// DEFAULT VALUES
	p.in_filename                 = nullptr;
	p.Ptrotter_number             = 0;
	p.write_PI_Frames_at_corrtime = false;



	if (argc >= 2) {

		while (n < argc) {

			std::istringstream issOptionToken(argv[n]);

			// report acceptance/rejection 
			if (!strncmp(issOptionToken.str().c_str(), "-P", 3)) {
				n++;
				if (SafeOps::atoi(argv[n], p.Ptrotter_number)) {
					n++;
					continue;
				} 
			}

			// write xyz frames for visualization (at corr time)
			else if (!strncmp(issOptionToken.str().c_str(), "-xyz", 5)) {
				n++;
				p.write_PI_Frames_at_corrtime = true;
			}

			else // if token is not a flag (-x) then it must be the filename
			{
				// INPUT FILENAME
				if (p.in_filename) {
					// if filename is already set, then the command is malformed. There should
					// only be one argument that is not preceeded by an option flag (e.g. "-P")
					sprintf(linebuf, "ERROR: Multiple input files specified(or unrecognized option \"%s\")\n", argv[n]);
					Output::err(linebuf);
					displayUsageAndDie(argv[0], p);
				}

				issOptionToken >> fname;
				if (!issOptionToken)
					displayUsageAndDie(argv[0], p);
				p.in_filename = (char*)calloc(fname.size() + 1, sizeof(char));
				if (p.in_filename) 
					std::memcpy(p.in_filename, fname.c_str(), fname.size());
				else {
					Output::err("Unable to allocate filename buffer.\n");
					die(fail);
				}
				n++;
			}

		} // end arg processing
	} // end arg count

	else displayUsageAndDie(argv[0], p);
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



// Removes the path from a filename, leaving on the filename
char* stripPath(char* full_path) {

	// Select either MS Windows or Unix file separator
	#ifndef _WIN32
		const char file_separator = '/';
	#else
		const char file_separator = '\\';
	#endif

	char* simple_filename = full_path;
	char* path            = full_path;

	// Skip to the last separator in the filename
	while (*path != '\0') {
		path++;
		if ((*path) == file_separator)
			simple_filename = path + 1;
	}

	return simple_filename;
}


#endif	/* ARGS_H */