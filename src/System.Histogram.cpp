#include "SafeOps.h"
#include "System.h"
#include "Output.h"

extern int rank;


void System::setup_histogram() 
// Setup the various quantities that define the histogram grid 
{

	char   linebuf[256]  = {0};
	double trial_vec1[3] = {0};
	double trial_vec2[3] = {0};
	double magA=0,
	       magB=0,
	       magC=0;
	int    Nbins=0,
	       x_dim=0,
	       y_dim=0,
	       z_dim=0;

	// get the magnitudes of all the basis vectors and test the frac2cart
	// routine. Define a fractional vector (1,0,0) and transform it with 
	// our basis. Then calculate its magnitude. Do this in all 3 dimensions.
	trial_vec1[0]=1.0;	trial_vec1[1]=0.0;	trial_vec1[2]=0.0;
	frac2cart( trial_vec2, trial_vec1 );
	magA = magnitude( trial_vec2 );

	trial_vec1[0]=0.0;	trial_vec1[1]=1.0;	trial_vec1[2]=0.0;
	frac2cart( trial_vec2, trial_vec1 );
	magB = magnitude( trial_vec2 );

	trial_vec1[0]=0.0;	trial_vec1[1]=0.0;	trial_vec1[2]=1.0;
	frac2cart( trial_vec2, trial_vec1 );
	magC = magnitude( trial_vec2 );

	// calculate the number of bins in each fractional coordinate 
	x_dim = (int) rint( magA / hist_resolution);
	y_dim = (int) rint( magB / hist_resolution);
	z_dim = (int) rint( magC / hist_resolution);
	sprintf( linebuf, "HISTOGRAM: %f resolution -> %d bins(A) * %d bins(B) * %d bins(C)\n", hist_resolution, x_dim, y_dim, z_dim);
	Output::out( linebuf ); 
	
	Nbins=x_dim * y_dim * z_dim;
	sprintf(linebuf,"HISTOGRAM: Total Bins = %d\n",Nbins);
	Output::out( linebuf );


	grids->histogram->x_dim = x_dim;
	grids->histogram->y_dim = y_dim;
	grids->histogram->z_dim = z_dim;
	grids->avg_histogram->x_dim = x_dim;
	grids->avg_histogram->y_dim = y_dim;
	grids->avg_histogram->z_dim = z_dim;
	grids->histogram->n_data_points = Nbins;
	n_histogram_bins = Nbins;
	grids->avg_histogram->norm_total=0;
	grids->histogram->norm_total=0;

	setup_dx_variables( grids->histogram );
	allocate_histogram_grid();
}



// to use this routine, compile with -g defined.
 // call this routine at any point in the execution
 // that you would like to attach.  when execution
 // begins, the process should write out the PIDS
 // of all processes and hang. open a new terminal
 // and start gdb. Attach to the PID of choice. go
 // up the function stack until you are in the scope 
 // of the 
 // gdb> attach 21456
 
/*
void attach(void)
{
	// debug routine to attach to gdb 
	static int i=0;
	printf("RANK: %d PID: %d\n",rank,getpid());
	fflush(NULL);
	while(!i) sleep(1);
}

*/



void System::update_root_histogram() 
// update the avg_histogram stored on root 	
{
	int xdim = grids->histogram->x_dim;
	int ydim = grids->histogram->y_dim;
	int zdim = grids->histogram->z_dim;

	for( int k=0; k < zdim; k++ ){
		for( int j=0; j < ydim; j++ ){
			for( int i=0; i < xdim; i++ ){
				grids->avg_histogram->grid[i][j][k] += grids->histogram->grid[i][j][k];
				// norm_total is updated here to normalize upon write out 
				grids->avg_histogram->norm_total    += grids->histogram->grid[i][j][k];
			}
		}
	}
}




void System::zero_grid( int ***grid ) 
// zero out the histogram grid 
{
	int xdim = grids->histogram->x_dim;
	int ydim = grids->histogram->y_dim;
	int zdim = grids->histogram->z_dim;

	for( int k=0; k < zdim; k++){
		for( int j=0; j < ydim; j++){
			for( int i=0; i < xdim; i++){
				grid[i][j][k]=0;
			}
		}
	}
}




void System::compute_bin( double *cart_coords, int *bin_vector ) 
// compute the histogram bin number to place this molecule in 
{
	double frac_coords[3] = {0};
	int    Abin           =  0,
	       Bbin           =  0,
	       Cbin           =  0;


	// we need the fractional coords to simplify the bin calculation
	cart2frac( frac_coords, cart_coords );

	// the coordinate system in the simulation is from -0.5 to 0.5 
	// so we need to correct for this: add 0.5 to each dimension 
	frac_coords[0] += 0.5;
	frac_coords[1] += 0.5;
	frac_coords[2] += 0.5;
	
	// compute bin in each dimension 
	Abin = (int)floor( frac_coords[0] * grids->histogram->x_dim );
	Bbin = (int)floor( frac_coords[1] * grids->histogram->y_dim );
	Cbin = (int)floor( frac_coords[2] * grids->histogram->z_dim );

	// return result to the bin_vector passed in 
	bin_vector[0] = Abin;
	bin_vector[1] = Bbin;
	bin_vector[2] = Cbin;
}



void System::wrap1coord( double *unwrapped, double *wrapped ) 
{
	double unit  [3] = {0},
	       offset[3] = {0},
	       frac  [3] = {0};

	// put coords in fractional representation 
	for( int i=0; i<3; i++ )
		for( int j=0; j<3; j++ ) 
			//we use transpose(recip_basis), because transpose(recip_basis).basis_vector = <1,0,0> , <0,1,0> or <0,0,1>
			frac[i] += pbc.reciprocal_basis[j][i]*unwrapped[j];

	// any fractional coord > .5 or < -.5 round to 1,-1 etc. 
	for( int i=0; i<3; i++ )
		unit[i] = rint(frac[i]);

	// multiply this rounded fractional unit vector by basis 
	for( int i=0; i<3; i++ )
		for( int j=0; j<3; j++ )
			offset[i] += pbc.basis[j][i] * unit[j];

	// subtract this distance from the incoming vector 
	for( int i=0; i<3; i++ )
		wrapped[i] = unwrapped[i] - offset[i];

}


void System::population_histogram() 
// Population histogram should be performed only every corr time. Individual nodes store the data on grids->histogram
// Only root will compile the sum into grids->avg_histogram. The node only grid->histogram will be rezeroed at every corr 
// time to prevent overflow. NOTE: The histogram is normalized to 1 (e.g. if there are a total of 328 bin counts -> every
// bin is divided by 328).
{
	Molecule  * mol_p            = nullptr;
	int         bin[3]           = {0};
	double      wrappedcoords[3] = {0};

	for( mol_p = molecules; mol_p; mol_p = mol_p->next ){
		if( ! mol_p->frozen ) {
			// wrap the coordinates of mol_p 
			wrap1coord( mol_p->com,wrappedcoords );
			// compute what bin to increment. store answer in bin[] 
			compute_bin( wrappedcoords, bin );
			// increment the bin returned in bin[] 
			(grids->histogram-> grid [(bin[0])] [(bin[1])] [(bin[2])] )++;
		} 
	}
}


void System::write_histogram( FILE *fp_out, int ***grid ) 
// This writes out the grid with the Cbin (last index) varying the fastest.  Line break between dimensions, double line break 
// between ZY sheets, ####'s between complete sets. Remember, for this to work, we had to offset the origin by 1/2 a bin width 
{
	int xdim  = 0,
	    ydim  = 0,
	    zdim  = 0,
	    count = 0;

	xdim = grids->histogram->x_dim;
	ydim = grids->histogram->y_dim;
	zdim = grids->histogram->z_dim;

	rewind( fp_out );
	fprintf( fp_out, "# OpenDX format population histogram\n" );
	fprintf( fp_out, "object 1 class gridpositions counts %d %d %d\n",xdim,ydim,zdim );
	fprintf( fp_out, "origin\t%f\t%f\t%f\n", grids->histogram->origin[0],   grids->histogram->origin[1],   grids->histogram->origin[2]   );
	fprintf( fp_out, "delta \t%f\t%f\t%f\n", grids->histogram->delta[0][0], grids->histogram->delta[0][1], grids->histogram->delta[0][2] );
	fprintf( fp_out, "delta \t%f\t%f\t%f\n", grids->histogram->delta[1][0], grids->histogram->delta[1][1], grids->histogram->delta[1][2] );
	fprintf( fp_out, "delta \t%f\t%f\t%f\n", grids->histogram->delta[2][0], grids->histogram->delta[2][1], grids->histogram->delta[2][2] );
	fprintf( fp_out, "\n" );
	fprintf( fp_out, "object 2 class gridconnections counts %d %d %d\n", xdim, ydim, zdim );
	fprintf( fp_out, "\n" );
	//This line is deprecated for viewing .dx files in current VMD software.
	//fprintf( fp_out, "object 3 class array type float rank 0 items %d data follows\n", grids->histogram->n_data_points );

	for( int i=0; i < xdim; i++ ) {
		for( int j=0; j < ydim; j++ ) {
			for( int k=0; k < zdim; k++ ) {
				fprintf( fp_out, "%f ", (float)(grid[i][j][k]) / (float)(grids->avg_histogram->norm_total) );
				count += grid[i][j][k];
			}
			fprintf(fp_out,"\n");
		}
		fprintf( fp_out, "\n" );
	}

	fprintf( fp_out, "# count=%d\n", count );
	fprintf( fp_out, "attribute \"dep\" string \"positions\"\n" );
	fprintf( fp_out, "object \"regular positions regular connections\" class field\n" );
	fprintf( fp_out, "component \"positions\" value 1\n" );
	fprintf( fp_out, "component \"connections\" value 2\n" );
	fprintf( fp_out, "component \"data\" value 3\n" );
	fprintf( fp_out, "\nend\n" );
	fflush(  fp_out );
}




void System::frac2cart(double *answer, double *frac ) 
// Take a vector in fractional coordinates (frac[]) and convert it to cartesian coordinates. Store the answer in answer[]     
{
	for( int i=0; i<3; i++ ) 
		answer[i]=0.0;	

	for(int i=0; i<3; i++ )
		for( int j=0; j<3; j++ )
			answer[i] += pbc.basis[j][i] * frac[j];
}




void System::cart2frac(double *answer, double *cart ) 
// Take a vector in cartesian coordinates (cart[]) and convert it to fractional coordinates. Store the answer in anwer[]
{
	for( int i=0; i<3; i++ )
		answer[i] = 0.0;

	for( int i=0; i<3; i++ )
		for(int j=0; j<3; j++ )
			//we use transpose(recip_basis), because transpose(recip_basis).basis_vector = <1,0,0> , <0,1,0> or <0,0,1>
			answer[i] += pbc.reciprocal_basis[j][i] * cart[j];

	return;
}



double System::magnitude( double *vector ) 
// Returns the magnitude of a 3-vector 
{
	double mag=0;
	for( int i=0; i<3; i++ )
		mag += vector[i] * vector[i];

	return sqrt(mag);
}



void System::allocate_histogram_grid()
{
	int x_dim = grids->histogram->x_dim;
	int y_dim = grids->histogram->y_dim;
	int z_dim = grids->histogram->z_dim;

	// allocate a 3D grid for the histogram 
	SafeOps::calloc( grids->histogram->grid, x_dim, sizeof(int **), __LINE__, __FILE__ );
	
	for( int i=0; i<x_dim; i++ ) 
		SafeOps::calloc( grids->histogram->grid[i], y_dim, sizeof(int *), __LINE__, __FILE__ );
	
	for( int i=0; i<x_dim; i++ ) 
		for( int j=0; j<y_dim; j++ ) 
			SafeOps::calloc( grids->histogram->grid[i][j], z_dim, sizeof(int), __LINE__, __FILE__ );

	// if root, allocate an avg_histogram grid 
	if( !rank ) {
		SafeOps::calloc( grids->avg_histogram->grid, x_dim, sizeof(int **), __LINE__, __FILE__ );
		
		for( int i=0; i<x_dim; i++) 
			SafeOps::calloc( grids->avg_histogram->grid[i], y_dim, sizeof(int *), __LINE__, __FILE__ );

		for( int i=0; i<x_dim; i++ ) 
			for( int j=0; j<y_dim; j++ ) 
				SafeOps::calloc( grids->avg_histogram->grid[i][j], z_dim, sizeof(int), __LINE__, __FILE__ );
	}


}




void System::offset_dx_origin( double *real_origin_cartesian, histogram_t *hist )
// Because OpenDX puts our data values at the points of the grid, we need to offset them by half a 
// bin width to reflect the fact that we have a histogram. Our data really lies between each point. 
{
	double fractional_binwidth[3];
	double cart_halfbin       [3];

	// figure out how wide each bin is in each dimension 
	fractional_binwidth[0] = 1.0 / hist->x_dim;
	fractional_binwidth[1] = 1.0 / hist->y_dim;
	fractional_binwidth[2] = 1.0 / hist->z_dim;
	
	// figure out how wide half a binwidth is in cartesians
	fractional_binwidth[0]/=2.0;
	fractional_binwidth[1]/=2.0;
	fractional_binwidth[2]/=2.0;
	frac2cart(cart_halfbin, fractional_binwidth );

	// add this value to the origin
	real_origin_cartesian[0] += cart_halfbin[0];
	real_origin_cartesian[1] += cart_halfbin[1];
	real_origin_cartesian[2] += cart_halfbin[2];
	
}	




void System::setup_dx_variables( histogram_t *hist )
// Variables needed upon printing in openDX native format. see DX userguide.pdf appendix B for details. 
{
	double vec   [3],
	       origin[3];

	// setup counts 
	hist->count[0] = hist->x_dim;
	hist->count[1] = hist->y_dim;
	hist->count[2] = hist->z_dim;

	// setup origin
	vec[0] = -0.5;
	vec[1] = -0.5;
	vec[2] = -0.5;
	frac2cart( origin, vec );

	// IMPORTANT!!! I am offsetting the origin by 1/2 a bin in each dimension 
	// the result of origin is not the true origin!!! 
	offset_dx_origin( origin, hist );

	hist->origin[0]=origin[0];
	hist->origin[1]=origin[1];
	hist->origin[2]=origin[2];

	// setup deltas
	setup_deltas( grids->histogram );

	// setup N data points
	hist->n_data_points = hist->x_dim * hist->y_dim * hist->z_dim;
}




void System::setup_deltas( histogram_t *hist )
// These are needed by dxwrite routines. Delta is a transformation matrix that defines the step
// size for setting up a grid in OpenDX (see the DX userguide.pdf appendix B for more details).
{
	for( int i=0; i<3; i++ )
		for( int j=0; j<3; j++ )
			hist->delta[i][j] = pbc.basis[j][i] / hist->count[i]; 
	// we divide by the count to get our actual step size in each dimension
	//Correction for non-orthorhombic systems
	hist->delta[1][0] = hist->delta[0][1];
	hist->delta[0][1] = 0.0;
	hist->delta[2][0] = hist->delta[0][2];
	hist->delta[0][2] = 0.0;
}