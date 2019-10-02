// Space Research Group
// Department of Chemistry
// University of South Florida


#include "System.h"


void System::mpi_copy_histogram_to_sendbuffer( char *snd, int ***grid ) {
	
	int x_dim = grids->histogram->x_dim;
	int y_dim = grids->histogram->y_dim;
	int z_dim = grids->histogram->z_dim;

	int * sndcast = (int *)snd; // cast the send buffer address to an (int *)

    // histogram
    for( int k=0; k<z_dim; k++ ) 
        for( int j=0; j<y_dim; j++ ) 
		    for( int i=0; i<x_dim; i++ ) 
				sndcast[i + j * x_dim + k * x_dim*y_dim] = grid[i][j][k];
				
}




void System::mpi_copy_rcv_histogram_to_data( char *rcv, int ***histogram )
{
	int x_dim = grids->histogram->x_dim;
	int y_dim = grids->histogram->y_dim;
	int z_dim = grids->histogram->z_dim;
	int *rcvcast=(int *)rcv;

        // histogram
        for( int k=0; k<z_dim; k++ ){
        	for( int j=0; j<y_dim; j++ ){
		        for( int i=0; i<x_dim; i++ ){
				histogram[i][j][k]=rcvcast[i + j*x_dim + k*x_dim*y_dim];
			}
		}
	}
}

