#pragma once
#ifndef ARGS_H
#define	ARGS_H

#include <cstring>
#include <iostream>
#include <sstream>
#include "Output.h"
#include "SafeOps.h"
#include "SimulationControl.h"

#ifdef __unix__
#include <unistd.h>
#include <signal.h>
#endif


#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 1
#endif
#ifndef EXIT_FAILURE 
#define EXIT_FAILURE 0
#endif


extern int rank;
extern int size;
extern bool mpi;

void die(int code);                                   // Kill any MPI Processes and stop execution
void install_signal_handler(SimulationControl *sc);   // Initialize signal handlers for clean exits (on Posix systems)
void signal_handler(int sigtype);                     // Function to handle interrupt signals (on Posix systems)



typedef struct _parameters {
	char *in_filename;
	bool  reportAR;
	bool  write_PI_Frames_at_corrtime;
} params;



// Includes for interrupt signal handler 
SimulationControl *sc;



char * stripPath(char * full_path) {

	// Select either MS Windows or Unix file separator
#ifndef _WIN32
	const char file_separator = '/';
#else
	const char file_separator = '\\';
#endif

	char *simple_name = full_path;
	char *i = full_path;

	// Skip to the last separator in the filename
	while (*i != '\0') {
		i++;
		if ((*i) == file_separator)
			simple_name = i + 1;
	}

	return simple_name;
}





void displayUsageAndDie(char* progname, params *p) {

	if (!rank) {
		char *progName_wo_path = stripPath(progname);
		std::cout << "\nUsage:" << std::endl;
		std::cout << "\t" << progName_wo_path << " INPUT_FILE [options]" << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "\t-xyz                 Write .xyz PI visualization frames at corrtime steps." << std::endl;
		std::cout << "\nWhen using MPI with path integral ensembles, the number of MPI processes will determine the" << std::endl;
		std::cout << "number of path integral \"beads\" (i.e. the Trotter number). For multi-site/multi-atom sorbate" << std::endl;
		std::cout << "molecules, the Trotter number must be a power of 2 greater than or equal to 4." << std::endl;
		std::cout << "\nExample:\n\tmpirun -n 8 " << progName_wo_path << " my_path_integral_sim_file" << std::endl;
		std::cout << "\nSee https://github.com/mpmccode/mpmc for input file specification." << std::endl;
		#ifndef _MPI
			std::cout << "\nNOTE: MPI Functionality is NOT compiled into this executable." << std::endl;
		#endif
	}
	if (p->in_filename)
		SafeOps::free(p->in_filename);

	die(EXIT_FAILURE);
}





//////////////////////////////////////////////////////////////
////    ADD SUPPORT FOR CLEAN UP IN signal_handler()     ////
////////////////////////////////////////////////////////////

// Install the signal handler for clean exits (on Unix systems)
void install_signal_handler(SimulationControl *simControl) {
	sc = simControl;
#if defined( _POSIX_VERSION )  
	signal(SIGTERM, signal_handler);
	signal(SIGUSR1, signal_handler);
	signal(SIGUSR2, signal_handler);
#endif
	return;
}




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




// Kill MPI before quitting, when neccessary
void die(int code) {

	#ifdef _MPI
	if (mpi) { MPI_Finalize(); }
	#endif

	if (code)
		exit(EXIT_SUCCESS);
	else
		exit(EXIT_FAILURE);
}





void processArgs(int argc, char * argv[], params *p) {
	std::string fname;
	int n = 1;

	// DEFAULT VALUES
	p->in_filename = nullptr;
	p->reportAR = false;
	p->write_PI_Frames_at_corrtime = false;



	if (argc >= 2) {

		while (n < argc) {

			std::istringstream issOptionToken(argv[n]);

			// report acceptance/rejection 
			if (!strncmp(issOptionToken.str().c_str(), "-r", 3)) {
				n++;
				p->reportAR = true;
			}

			// write xyz frames for visualization (at corr time)
			else if (!strncmp(issOptionToken.str().c_str(), "-xyz", 5)) {
				n++;
				p->write_PI_Frames_at_corrtime = true;
			}

			else // if token is not a flag (-x) then it must be the filename
			{
				// INPUT FILENAME
				if (p->in_filename) {
					// if filename is already set, then the command is malformed. There should
					// only be one argument that is not preceeded by an option flag (e.g. "-P")
					if (!rank) std::cout << "ERROR: Multiple input files specified (or unrecognized option \"" << argv[n] << "\")" << std::endl;
					displayUsageAndDie(argv[0], p);
				}

				issOptionToken >> fname;
				if (!issOptionToken)
					displayUsageAndDie(argv[0], p);
				p->in_filename = (char *)calloc(fname.size() + 1, sizeof(char));
				std::memcpy(p->in_filename, fname.c_str(), fname.size());

				n++;
			}

		} // end arg processing
	} // end arg count
	else displayUsageAndDie(argv[0], p);

}
#endif	/* ARGS_H */

