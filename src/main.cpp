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
	params args;
	char linebuf[maxLine] = { 0 };

	time_t t = time(nullptr);
	struct tm tm = *localtime(&t);

	mpi_introspection_and_initialization(argc, argv); // detect/start MPI services
	processArgs(argc, argv, &args); // Parse command line

	sprintf(  linebuf, "MPMC++\nMassively Parallel Monte Carlo: Multi-System Edition%s, v%s -- 2012-2019 GNU Public License\n", (size>0) ? " (MPI enabled)" : "", VERSION);
	sprintf( &linebuf[strlen(linebuf)], "MAIN: processes started on %d threads(s) @ %d-%d-%d %d:%d:%d\n", mpi ? size : 1, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
	Output::out1(linebuf);

	
	try {
		// Read the input file
		SimulationControl simController(args.in_filename, args.reportAR, args.write_PI_Frames_at_corrtime); 
		Output::out1("MAIN: Simulation parameters established.\n");

		install_signal_handler(&simController);

		// Check simulation parameters, read geometry and allocate/populate required data structures
		simController.initializeSimulationObjects();
		Output::out1("MAIN: Input parameters checked. System data structures allocated and initialized.\n");

		// Execute the simulation
		simController.runSimulation();

	}
	catch (int e) {

		sprintf(linebuf, "MPMC exiting with error code: %d.\n", e);
		Output::err(linebuf);
		if (args.in_filename) { SafeOps::free(args.in_filename); }
		die(fail);
	}

	die(ok);
}
