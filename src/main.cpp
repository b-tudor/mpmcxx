// Space Research Group
// Department of Chemistry
// University of South Florida

#ifndef VERSION
#define VERSION "0.81"
#endif

int rank = 0;
int size = 0;
bool mpi = false;

#include <stdio.h>
#include <time.h>
#include "args_etc.h"
#include "constants.h"
#include "Output.h"
#include "SafeOps.h"
#include "SimulationControl.h"



int main(int argc, char * argv[])
{
	//  Parse command line
	params args;
	processArgs(argc, argv, args);

	//  Detect/start MPI services 
	mpi_introspection_and_initialization(argc, argv, args.Ptrotter_number);

	//  Say hello
	introduce_self();

	
	
	
	try {

		// Read the input file
		SimulationControl simController(args.in_filename, args.Ptrotter_number, args.write_PI_Frames_at_corrtime ); 
		Output::out1("MAIN: Simulation parameters established.\n");

		install_signal_handler(&simController);

		// Check simulation parameters, read geometry and allocate/populate required data structures
		simController.initializeSimulationObjects();
		Output::out1("MAIN: Input parameters checked. System data structures allocated and initialized.\n");

		// Execute the simulation
		simController.runSimulation();

	}
	catch (int e) {

		char linebuf[maxLine] = { 0 };
		sprintf(linebuf, "MPMC exiting with error code: %d.\n", e);
		Output::err(linebuf);

		if (args.in_filename) { SafeOps::free(args.in_filename); }

		die(fail);
	}

	die(ok);
}
