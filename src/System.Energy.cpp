#include <cstring>
#include <math.h>

#include "Atom.h"
#include "Molecule.h"
#include "Output.h"
#include "Pair.h"
#include "SafeOps.h"
#include "System.h"
#include "UsefulMath.h"
#include "Vector3D.h"






// returns the total potential energy for the system and updates our observables
double System::energy() {

	double potential_energy  = 0, 
	       rd_energy         = 0, 
	       coulombic_energy  = 0,
	       polar_energy      = 0,
	       vdw_energy        = 0,
	       three_body_energy = 0,
	       kinetic_energy    = 0;

	#ifdef POLARTIMING
		static double timing = 0;
		static int count = 0;
	#endif

	natoms = countNatoms();

	// get the pairwise terms necessary for the energy calculation 
	pairs();

	// Only on the first simulation step, make sure that all recalculate flags are set if we made a 
	// volume change (or just reverted from one) set recalculate flags OR if replaying a trajectory 
	// we set last_volume at the end of this function
	if(   last_volume != pbc.volume   ||   ensemble == ENSEMBLE_REPLAY   ||   observables->energy == 0.0  )
		flag_all_pairs();

	// get the electrostatic potential
	if(  ! ( use_sg || rd_only )  ) {

		if( spectre )
			coulombic_energy = coulombic_nopbc( molecules );
		else if( gwp ) {
			coulombic_energy = coulombic_nopbc_gwp();
			kinetic_energy   = coulombic_kinetic_gwp();
			observables->kinetic_energy = kinetic_energy;
		} else
			coulombic_energy = coulombic();

		observables->coulombic_energy = coulombic_energy;

		// get the polarization potential
		if( polarization ) {

			#ifdef POLARTIMING
				// get timing of polarization energy function for cuda comparison 
				Output::GetTimeOfDay( &old_time );
			#endif

			#ifdef CUDA
				if(system->cuda)
					polar_energy = (double)polar_cuda(system);
				else
					polar_energy = polar(system);
			#else
				polar_energy = polar();
			#endif // CUDA 

			#ifdef POLARTIMING
				Output::GetTimeOfDay( &new_time );
				timing = timing * (double)count/((double)count+1.0) 
					+ (double)((new_time.tv_sec-old_time.tv_sec)*1e6+(new_time.tv_usec-old_time.tv_usec)) * 1.0/((double)count +1.0);
				count++;
				if ( system->corrtime ) {
					if ( count % system->corrtime == 0 ) sprintf(linebuf, "OUTPUT: Polarization energy function took %lf us\n", timing);
					output(linebuf);
				}
				else	{
					sprintf(linebuf, "OUTPUT: Polarization energy function took %lf us\n", timing);
					output(linebuf);
				}
			#endif

			observables->polarization_energy = polar_energy;

		}
		if( polarvdw ) {
			#ifdef CUDA
				if (system->cuda) {
					error("error: cuda polarvdw not yet implemented!\n");
					die(-1);
				}
				else
				vdw_energy = vdw(system);
			#else
				vdw_energy = vdw();
			#endif
				observables->vdw_energy = vdw_energy;
		}

	}



	// get the repulsion/dispersion potential
	if( rd_anharmonic )
		rd_energy = anharmonic();
	else if( use_sg )
		rd_energy = sg();
	else if( use_dreiding )
		rd_energy = dreiding();
	else if( using_lj_buffered_14_7 )
		rd_energy = lj_buffered_14_7();
	else if( using_disp_expansion )
		rd_energy = disp_expansion();
	else if( cdvdw_exp_repulsion )
		rd_energy = exp_repulsion();
	else if( !gwp )
		rd_energy = lj();
	observables->rd_energy = rd_energy;
    
	if( using_axilrod_teller )
	{
		three_body_energy = axilrod_teller();
		observables->three_body_energy = three_body_energy;
	}

	// sum the total potential energy 
	potential_energy = rd_energy + coulombic_energy + polar_energy + vdw_energy + three_body_energy;
	
	// not truly potential, but stick it there for convenience of MC 
	//
	//    POSSIBLE BUG: kinetic_energy was uninitialized, and previously only given a value inside the conditional: 
	//
	//        if(!(system->use_sg || system->rd_only)) { HERE }
	//
	//    If this conditional fails, but (system->gwp) is true (idk if this is possible), an un-initialized value would have been
	//    added to potential_energy. Now, 0 is being added, but am not sure if this is the desired behavior. -bt
	 
	
	if( gwp )
		potential_energy += kinetic_energy;
	observables->energy = potential_energy;

	countN();
	observables->spin_ratio /= observables->N;

	// for NVE
	if(ensemble == ENSEMBLE_NVE) {
		observables->kinetic_energy = total_energy - potential_energy;
		observables->temperature = (2.0/3.0) * observables->kinetic_energy/observables->N;
	}

	// need this for the isosteric heat 
	observables->NU = observables->N * observables->energy;

	// set last known volume
	last_volume = pbc.volume;

	if( cavity_autoreject_absolute )
		potential_energy += cavity_absolute_check();

	return potential_energy;
}



double System::vdw() 
//returns interaction VDW energy
{

	int       NAtoms = natoms; // number of atoms
	double    e_total,         // total energy
	          e_iso;           // isolation energy (atoms @ infinity)
	double  * sqrtKinv;        // matrix K^(-1/2); cholesky decomposition of K
	double ** Am = A_matrix;   // A_matrix
	mtx_t   * Cm;              // C_matrix (we use single pointer due to LAPACK requirements)
	double  * eigvals;         // eigenvales
	double    fh_corr, lr_corr;


	//allocate arrays. sqrtKinv is a diagonal matrix. d,e are used for matrix diag.
	sqrtKinv = getsqrtKinv( NAtoms );

	//calculate energy vdw of isolated molecules
	e_iso = sum_eiso_vdw ( sqrtKinv );

	//Build the C_Matrix
	Cm = build_M (3*NAtoms, 0, Am, sqrtKinv);

	//setup and use lapack diagonalization routine dsyev_()
	eigvals = lapack_diag (Cm, polarvdw ); //eigenvectors if system->polarvdw == 2
	if( polarvdw == 2 )
		printevects(Cm);

	//return energy in inverse time (a.u.) units
	e_total = eigen2energy(eigvals, Cm->dim); // , temperature );
	e_total *= au2invseconds * half_hBar; //convert a.u. -> s^-1 -> K

	//vdw energy comparison
	if ( polarvdw == 3 )
		printf("VDW Two-Body | Many Body = %lf | %lf\n", twobody(), e_total-e_iso );

	if( feynman_hibbs ) {
		if( vdw_fh_2be ) fh_corr = fh_vdw_corr_2be(); //2be method
		else fh_corr = fh_vdw_corr(); //mpfd
	}
	else fh_corr=0;

	if ( rd_lrc ) lr_corr = lr_vdw_corr();
	else lr_corr=0;

	//cleanup and return
	free(sqrtKinv);
	free(eigvals);
	free_mtx(Cm);

	return e_total - e_iso + fh_corr + lr_corr;

}


//build the matrix K^(-1/2) -- see the PDF
double * System::getsqrtKinv( int NAtoms ) {
	double   * sqrtKinv;
	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	int i = 0;

	//malloc 3*N wastes an insignificant amount of memory, but saves us a lot of index management
	 SafeOps::malloc( sqrtKinv, 3 * NAtoms * sizeof(double), __LINE__, __FILE__ );
	
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr=molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			//Seek through atoms, calculate sqrtKinv*
			sqrtKinv[i] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			sqrtKinv[i+1] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			sqrtKinv[i+2] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			i+=3;
		}
	}
	
	return sqrtKinv;
}


//go through each molecule and determine the VDW energy associated with each isolated molecule
double System::sum_eiso_vdw ( double * sqrtKinv ) {

	char       linebuf[maxLine];
	double     e_iso = 0;
	Molecule * mp;
	vdw_t    * vp;
	vdw_t    * vpscan;

	//loop through molecules. if not known, calculate, store and count. otherwise just count.
	for ( mp = molecules; mp; mp=mp->next ) {
		for ( vp = vdw_eiso_info; vp != NULL; vp=vp->next ) { //loop through all vp's
			if ( strncmp(vp->mtype,mp->moleculetype, maxLine ) == 0 ) {
					e_iso += vp->energy; //count energy
					break; //break out of vp loop. the current molecule is accounted for now. go to the next molecule
			}
			else continue; //not a match, check the next vp
		} //vp loop

		if ( vp == NULL ) { //if the molecule was unmatched, we need to grow the list
			// end of vp list and we haven't matched yet -> grow vdw_eiso_info
			// scan to the last non-NULL element
			if ( vdw_eiso_info == NULL ) {
				SafeOps::calloc( vdw_eiso_info, 1, sizeof(vdw_t), __LINE__, __FILE__ );
				vpscan = vdw_eiso_info; //set scan pointer
			} else {
				for ( vpscan = vdw_eiso_info; vpscan->next != NULL; vpscan=vpscan->next );
				SafeOps::calloc( vpscan->next, 1, sizeof(vdw_t), __LINE__, __FILE__ ); //allocate space
				vpscan = vpscan->next;
			} //done scanning and malloc'ing
		
			//set values
			strncpy( vpscan->mtype, mp->moleculetype, maxLine ); //assign moleculetype
			vpscan->energy = calc_e_iso( sqrtKinv, mp ); //assign energy
			if ( std::isfinite(vpscan->energy) == 0 ) { //if nan, then calc_e_iso failed
				sprintf(linebuf,"VDW: Problem in calc_e_iso.\n");
				Output::out( linebuf );
				throw infinite_energy_calc;
			}
			//otherwise count the energy and move to the next molecule
			e_iso += vpscan->energy;

		} //vp==NULL
	} //mp loop	

	////all of this logic is actually really bad if we're doing surface fitting, since omega will change... :-(
	//free everything so we can recalc next step
	if( ensemble == ENSEMBLE_SURF_FIT )  {
		free_vdw_eiso( vdw_eiso_info );
		vdw_eiso_info = NULL;
	}
	
	return e_iso;
}

//calculate energies for isolated molecules
//if we don't know it, calculate it and save the value
double System::calc_e_iso ( double * sqrtKinv, Molecule * mptr ) {

	int        nstart, nsize; 
	double     e_iso;         //total vdw energy of isolated molecules
	mtx_t    * Cm_iso;        //matrix Cm_isolated
	double   * eigvals;       //eigenvalues of Cm_cm
	Molecule * molecule_ptr;
	Atom     * atom_ptr;

	nstart=nsize=0; //loop through each individual molecule
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr=molecule_ptr->next ) {
		if ( molecule_ptr != mptr ) {  //count atoms then skip to next molecule
			for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) nstart++;
			continue;
		}

		//now that we've found the molecule of interest, count natoms, and calc energy
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) nsize++;

		//build matrix for calculation of vdw energy of isolated molecule
		Cm_iso  = build_M(3*(nsize), 3*nstart, A_matrix, sqrtKinv);
		//diagonalize M and extract eigenvales -> calculate energy
		eigvals = lapack_diag( Cm_iso, 1 ); //no eigenvectors
		e_iso = eigen2energy(eigvals, Cm_iso->dim); // , temperature );

		//free memory
		free(eigvals);
		free_mtx(Cm_iso);

		//convert a.u. -> s^-1 -> K
		return e_iso * au2invseconds * half_hBar;
	}	

	//unmatched molecule
	return NAN; //we should never get here
}


//free vdw pointer which keeps track of e_iso energies
void System::free_vdw_eiso(vdw_t * vdw_eiso_info) {

	vdw_t * vp;
	vdw_t ** varray = NULL;
	int i=0;

	for ( vp = vdw_eiso_info; vp; vp=vp->next ) {
		SafeOps::realloc( varray, sizeof(vdw_t *)*(i+1), __LINE__, __FILE__ );
		varray[i]=vp;
		i++;
	}

	while ( i>0 ) {
		i--;
		free(varray[i]);
	}

	free(varray);

	return;
}


//build C matrix for a given molecule/system, with atom indicies (offset)/3..(offset+dim)/3
System::mtx_t * System::build_M ( int dim, int offset, double ** Am, double * sqrtKinv ) {
	int i; //dummy
	int iA, jA; //Am indicies
	int iC, jC; //Cm indicies
	int nonzero; //non-zero col/rows in Am
	mtx_t * Cm; //return matrix Cm

	//count non-zero elements
	nonzero=0;
	for ( i=offset; i<dim+offset; i++ )
		if ( sqrtKinv[i] != 0 ) nonzero++;

	//allocate
	Cm = alloc_mtx(nonzero);

	// build lapack compatible matrix from Am[offset..dim, offset..dim]
	iC=jC=-1; //C index
	for ( iA=offset; iA<dim+offset; iA++ ) {
		if ( sqrtKinv[iA] == 0 ) continue; //suppress rows/cols full of zeros
		iC++; jC=-1;
		for ( jA=offset; jA<=iA; jA++ ) {
			if ( sqrtKinv[jA] == 0 ) continue; //suppress
			jC++;
			(Cm->val)[iC+jC*(Cm->dim)]=
				Am[iA][jA]*sqrtKinv[iA]*sqrtKinv[jA];
		}
	}

	return Cm;
}



System::mtx_t * System::alloc_mtx( int dim ) {
	//alloc matrix variable and set dim
	mtx_t * M = nullptr;
	SafeOps::malloc( M, sizeof(mtx_t), __LINE__, __FILE__ );
	M->dim=dim;
	//alloc matrix storage space
	SafeOps::calloc( M->val, dim*dim, sizeof(double), __LINE__, __FILE__ );
	
	return M;
}



void System::free_mtx( mtx_t * M ) {
	free(M->val);
	free(M);
	return;
}



void System::printevects( mtx_t * M ) {
	int R, C;

	printf("%%vdw === Begin Eigenvectors ===\n");
	for( R=0; R < (M->dim); R++ ) {
		for ( C=0; C < (M->dim); C++ ) {
			printf("%.2le ", (M->val)[R+C*M->dim]);
		}
		printf("\n");
	}
	printf("%%vdw=== End Eigenvectors ===\n");
}


//double wtanh ( double w, double T ) {
// not needed unless T >> 300 
//	if ( w < 1.0E-10 )
//		TWOoverHBAR*T/au2invsec; //from Taylor expansion
//	if ( T == 0 ) return w;
//	return w/tanh(halfHBAR*w*au2invsec/T);
//}

double System::eigen2energy( double * eigvals, int dim ) { //, double temperature ) {
	int i;
	double rval=0;
	
	if ( eigvals == NULL ) return 0;

	for ( i=0; i<dim; i++ ) {
		if ( eigvals[i] < 0 ) eigvals[i]=0;
		//rval += wtanh(sqrt(eigvals[i]), temperature);
		rval += sqrt(eigvals[i]);
	}
	return rval;
}


double System::twobody() 
//with damping
{
	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;
	double     energy = 0;

	//for each pair
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
				//skip if frozen
				if ( pair_ptr->frozen ) continue;
				//skip if they belong to the same molecule
				if ( molecule_ptr == pair_ptr->molecule ) continue;
				//skip if distance is greater than cutoff
				if ( pair_ptr->rimg > pbc.cutoff ) continue;
				//check if fh is non-zero
				if ( atom_ptr->polarizability == 0 || pair_ptr->atom->polarizability == 0 || 
					atom_ptr->omega == 0 || pair_ptr->atom->omega == 0 ) continue; //no vdw energy

				//calculate two-body energies
				energy += e2body(atom_ptr,pair_ptr,pair_ptr->rimg);
			}
		}
	}

	return energy;
}


//calculate T matrix element for a particular separation
double System::e2body( Atom * atom, Pair * pair, double r ) {
	double energy;
	double lr  = polar_damp * r;
	double lr2 = lr*lr;
	double lr3 = lr*lr2;
	double Txx = pow(r,-3)*(-2.0+(0.5*lr3+lr2+2*lr+2)*exp(-lr));
	double Tyy = pow(r,-3)*(1-(0.5*lr2+lr+1)*exp(-lr));
	double * eigvals;
	mtx_t * M = alloc_mtx(6);
	
	//only the sub-diagonals are non-zero
	M->val[1]=M->val[2]=M->val[4]=M->val[5]=M->val[6]=M->val[8]=M->val[9]=M->val[11]=0;
	M->val[12]=M->val[13]=M->val[15]=M->val[16]=M->val[19]=M->val[20]=M->val[22]=M->val[23]=0;
	M->val[24]=M->val[26]=M->val[27]=M->val[29]=M->val[30]=M->val[31]=M->val[33]=M->val[34]=0;

	//true diagonals
	M->val[0]=M->val[7]=M->val[14]=(atom->omega)*(atom->omega);
	M->val[21]=M->val[28]=M->val[35]=(pair->atom->omega)*(pair->atom->omega);

	//sub-diagonals
	M->val[3]=M->val[18]=
		(atom->omega)*(pair->atom->omega)*sqrt(atom->polarizability*pair->atom->polarizability)*Txx;
	M->val[10]=M->val[17]=M->val[25]=M->val[32]=
		(atom->omega)*(pair->atom->omega)*sqrt(atom->polarizability*pair->atom->polarizability)*Tyy;

	eigvals=lapack_diag(M,1);
	energy = eigen2energy(eigvals, 6); // , temperature);

	//subtract energy of atoms at infinity
	//energy -= 3*wtanh(atom->omega, system->temperature);
	energy -= 3*atom->omega;
	//energy -= 3*wtanh(pair->atom->omega, system->temperature);
	energy -= 3*pair->atom->omega;

	free(eigvals);
	free_mtx(M);

  return energy * au2invseconds * half_hBar;
}



double * System::lapack_diag ( mtx_t * M, int jobtype )
//    LAPACK using 1D arrays for storing matricies.
//    / 0  3  6 \
//    | 1  4  7 |    =    [ 0 1 2 3 4 5 6 7 8 ]
//    \ 2  5  8 /									
{
	
	char     job;      //job type
	char     uplo;
	double * work;     //working space for dsyev
	int      lwork;    //size of work array
	int      rval=0;   //returned from dsyev_
	double * eigvals;
	char     linebuf[maxLine];

	uplo = 'L';                      //operate on lower triangle
	job = (jobtype==2) ? 'V' : 'N';  //eigenvectors or no?
	
	if ( M->dim == 0 ) return NULL;

	//allocate eigenvalues array
	SafeOps::malloc( eigvals, M->dim*sizeof(double), __LINE__, __FILE__ );
	
	//optimize the size of work array
	lwork = -1;
	SafeOps::malloc( work, sizeof(double), __LINE__, __FILE__ );
	//dsyev_(&job, &uplo, &(M->dim), M->val, &(M->dim), eigvals, work, &lwork, &rval); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//now optimize work array size is stored as work[0]
	lwork = (int)work[0];
	SafeOps::realloc( work, lwork*sizeof(double), __LINE__, __FILE__ );
	//diagonalize
	//dsyev_(&job, &uplo, &(M->dim), M->val, &(M->dim), eigvals, work, &lwork, &rval); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if ( rval != 0 ) {
		sprintf(linebuf,"error: LAPACK: dsyev returned error: %d\n", rval);
		Output::err(linebuf);
		throw lapack_error;
	}

	free(work);

	return eigvals;
}


// long-range correction
double System::lr_vdw_corr() {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	double w1, w2;   //omegas
	double a1, a2;   //alphas
	double cC;       //leading coefficient to r^-6
	double corr = 0; //correction to the energy

	//skip if PBC isn't set-up
	if ( pbc.volume == 0 ) {
		Output::err( "VDW: PBC not set-up. Did you define your basis? Skipping LRC.\n" );
		return 0;
	}

	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
					//skip if frozen
					if ( pair_ptr->frozen ) continue;
					// skip if same molecule  // don't do this... this DOES contribute to LRC
					// if ( molecule_ptr == pair_ptr->molecule ) continue;
					// fetch alphas and omegas
					a1 = atom_ptr->polarizability;
					a2 = pair_ptr->atom->polarizability;
					w1 = atom_ptr->omega;
					w2 = pair_ptr->atom->omega;
					if ( w1 == 0 || w2 == 0 || a1 == 0 || a2 == 0 ) continue; //no vdw energy
					// 3/4 hbar/k_B(Ks) omega(s^-1)  Ang^6
					cC=1.5 * c_hBar * w1*w2/(w1+w2) * au2invseconds * a1 * a2;

					// long-range correction
					corr += -4.0/3.0 * pi * cC * pow(pbc.cutoff,-3) / pbc.volume;
			}
		}
	}

	return corr;
}


// feynman-hibbs correction - molecular pair finite differencing method
double System::fh_vdw_corr() {

	const double FINITE_DIFF = 0.01; //too small -> vdw calc noises becomes a problem

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;
	double     rm;                //reduced mass
	double     E[5];              //energy at five points, used for finite differencing
	double     dv, d2v, d3v, d4v; //derivatives
	double     corr = 0;          //correction to the energy
	double     corr_single;       //single vdw interaction energy
	double     H = FINITE_DIFF ;  //small dr used for finite differencing //too small -> vdw calculation noise becomes a problem

	//for each pair
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
				//skip if frozen
				if ( pair_ptr->frozen ) continue;
				//skip if they belong to the same molecule
				if ( molecule_ptr == pair_ptr->molecule ) continue;
				//skip if distance is greater than cutoff
				if ( pair_ptr->rimg > pbc.cutoff ) continue;
				//check if fh is non-zero
				if ( atom_ptr->polarizability == 0 || pair_ptr->atom->polarizability == 0 || 
					atom_ptr->omega == 0 || pair_ptr->atom->omega == 0 ) continue; //no vdw energy

				//calculate two-body energies
				E[0]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg-H-H); //smaller r
				E[1]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg-H); 
				E[2]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg); //current r
				E[3]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg+H); //larger r
				E[4]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg+H+H);

				//derivatives (Numerical Methods Using Matlab 4E 2004 Mathews/Fink 6.2)
				dv = (E[3]-E[1])/(2.0*H);
				d2v = (E[3]-2.0*E[2]+E[1])/(H*H);
				d3v = (E[4]-2*E[3]+2*E[1]-E[0])/(2*pow(H,3));
				d4v = (E[4]-4*E[3]+6*E[2]-4*E[1]+E[0])/pow(H,4);
				
				// reduced mass
				rm=AMU2KG*(molecule_ptr->mass)*(pair_ptr->molecule->mass)/
					((molecule_ptr->mass)+(pair_ptr->molecule->mass));

				//2nd order correction
				corr_single = pow(METER2ANGSTROM, 2) * (hBar*hBar/(24.0*kB*temperature*rm)) * (d2v + 2.0*dv/pair_ptr->rimg);
				//4th order correction
				if( feynman_hibbs_order >= 4 )
					corr_single += pow(METER2ANGSTROM, 4)*(pow(hBar, 4) /
						(1152.0*pow(kB*temperature*rm, 2))) *
						(15.0*dv/pow(pair_ptr->rimg, 3) + 4.0*d3v/pair_ptr->rimg + d4v);

				corr += corr_single;
			}
		}
	}

	return corr;
}


// feynman-hibbs using 2BE (shitty)
double System::fh_vdw_corr_2be() {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	double rm; //reduced mass
	double w1, w2; //omegas
	double a1, a2; //alphas
	double cC; //leading coefficient to r^-6
	double dv   = 0;
	double d2v  = 0;
	double d3v  = 0;
	double d4v  = 0; //derivatives
	double corr = 0; //correction to the energy
	double corr_single; //single vdw interaction energy

	//for each pair
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) { 
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) { 
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) { 
				//skip if frozen
				if ( pair_ptr->frozen ) continue;
				//skip if they belong to the same molecule
				if ( molecule_ptr == pair_ptr->molecule ) continue;
				//skip if distance is greater than cutoff
				if ( pair_ptr->rimg > pbc.cutoff ) continue;
				//fetch alphas and omegas
				a1=atom_ptr->polarizability;
				a2=pair_ptr->atom->polarizability;
				w1=atom_ptr->omega;
				w2=pair_ptr->atom->omega;

				if ( w1 == 0 || w2 == 0 || a1 == 0 || a2 == 0 ) continue; //no vdw energy
				// 3/4 hbar/k_B(Ks) omega(s^-1)  Ang^6
				cC=1.5 * c_hBar * w1*w2/(w1+w2) * au2invseconds * a1 * a2;
				// reduced mass
				rm = AMU2KG * (molecule_ptr->mass) * (pair_ptr->molecule->mass)/
					((molecule_ptr->mass)+(pair_ptr->molecule->mass));

				//derivatives 
				dv = 6.0*cC*pow(pair_ptr->rimg,-7);
				d2v= dv * (-7.0)/pair_ptr->rimg;
				if( feynman_hibbs_order >= 4 ) {
					d3v= d2v* (-8.0)/pair_ptr->rimg;
					d4v= d3v* (-9.0)/pair_ptr->rimg;
				}

				//2nd order correction
				corr_single = pow(METER2ANGSTROM, 2)*(hBar*hBar/(24.0*kB*temperature*rm))*(d2v + 2.0*dv/pair_ptr->rimg);
				//4th order correction
				if ( feynman_hibbs_order >= 4 )
					corr_single += pow(METER2ANGSTROM, 4)   *   (  pow(hBar, 4) / (1152.0 * pow(kB*temperature*rm, 2))  ) 
					               * ( 15.0 * dv/pow(pair_ptr->rimg, 3)   +   4.0*d3v/pair_ptr->rimg   +   d4v );
				corr += corr_single;
			}
		}
	}

	return corr;
}


// energy of molecule in a 1D anharmonic well
double System::anharmonic() {

	double     k, g, x;
	double     energy;
	Molecule * molecule_ptr;
	Atom     * atom_ptr;

	k = rd_anharmonic_k;
	g = rd_anharmonic_g;

	energy = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			x = atom_ptr->pos[0];
			energy += anharmonic_energy(k, g, x);

			if( feynman_hibbs ) {

				if( feynman_kleinert ) {

					energy = anharmonic_fk( temperature, atom_ptr->mass, k, g, x );

				} else {

					if( feynman_hibbs_order == 2 )
						energy += anharmonic_fh_second_order( temperature, atom_ptr->mass, k, g, x );
					else if( feynman_hibbs_order == 4 )
						energy += anharmonic_fh_fourth_order( temperature, atom_ptr->mass, k, g, x );
				}
			}
		}
	}

	return(energy);
}


// anharmonic potential
double System::anharmonic_energy(double k, double g, double x) {

	double potential;

	potential =   0.5 * k * pow(x, 2)    +    0.25 * g * pow(x, 4);

	return potential;
}


// Feynman-Kleinert iterative method of effective potential */
double System::anharmonic_fk(double temperature, double mass, double k, double g, double x) {

	int    keep_iterating;
	double a_sq;               // width a^2 (A^2) 
	double omega, omega_sq;    // spacing Omega^2 (K/A^2) 
	double prev_a_sq;          // last a_sq 
	double tolerance;          // iterative tolerance 
	double V_a;                // V_a^2 (K) 
	double potential;          // W_1 (K) 
	double conversion_factor;  // hbar^2*(m2A)^2/(k*m)

	// convert the mass to kg 
	mass *= AMU2KG;

	conversion_factor = pow(METER2ANGSTROM, 2)*pow(hBar, 2)/(kB*mass);

	// initial guess a^2 = beta/12 
	a_sq = pow(METER2ANGSTROM, 2)*pow(hBar, 2)/(12.0*kB*temperature*mass);

	// solve self-consistently 
	keep_iterating = 1;
	while(keep_iterating) {

		// save the last a_sq for tolerance 
		prev_a_sq = a_sq;

		omega_sq = conversion_factor*(k + 3.0*g*a_sq + 3.0*g*pow(x, 2)); omega = sqrt(omega_sq);
		a_sq = conversion_factor*(temperature/omega_sq)*((omega/(2.0*temperature))*(1.0/tanh(omega/(2.0*temperature))) - 1.0);

		tolerance = fabs(prev_a_sq - a_sq);
		if(tolerance < FEYNMAN_KLEINERT_TOLERANCE)
			keep_iterating = 0;

	}

	V_a = 0.5*a_sq*k + 0.75*g*pow(a_sq, 2) + 0.5*(k + 3.0*g*a_sq)*pow(x, 2) + 0.25*g*pow(x, 4);

	potential = temperature*log(sinh(omega/(2.0*temperature))/(omega/(2.0*temperature))) - 0.5*omega_sq*a_sq/conversion_factor + V_a;

	return potential;
}


// up to FH h^2 effective potential term 
double System::anharmonic_fh_second_order(double temperature, double mass, double k, double g, double x) {

	double first_derivative;
	double second_derivative;
	double potential;

	mass *= AMU2KG;

	first_derivative = k*x + g*pow(x, 3);
	second_derivative = k + 3.0*g*pow(x, 2);

	potential = pow(METER2ANGSTROM, 2)*pow(hBar, 2)/(24.0*kB*temperature*mass)*(second_derivative + 2.0*first_derivative/x);

	return potential;
}


// up to FH h^4 effective potential term 
double System::anharmonic_fh_fourth_order(double temperature, double mass, double k, double g, double x) {

	double first_derivative, second_derivative;
	double other_derivatives;
	double potential;

	mass *= AMU2KG;

	first_derivative = k*x + g*pow(x, 3);
	second_derivative = k + 3.0*g*pow(x, 2);
	other_derivatives = 15.0*k/pow(x, 2) + 45.0*g;

	potential = pow(METER2ANGSTROM, 2)*pow(hBar, 2)/(24.0*kB*temperature*mass)*(second_derivative + 2.0*first_derivative/x);
	potential += pow(METER2ANGSTROM, 4)*pow(hBar, 4)/(1152.0*pow(kB*temperature*mass, 2))*other_derivatives;

	return potential;
}






//   Lennard-Jones repulsion/dispersion    //////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



double System::lj()
{

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	Pair     * pair_ptr     = nullptr;
	double     sigma_over_r        = 0,
	           term12              = 0,
	           term6               = 0,
	           sigma_over_r6       = 0,
	           sigma_over_r12      = 0,
	           r                   = 0,
	           potential           = 0,
	           potential_classical = 0,
	           cutoff              = 0;
	int        i[3] = { 0, 0, 0 };
	double     a[3] = { 0, 0, 0 };

	//set the cutoff
	if (rd_crystal)
		cutoff = 2.0 * pbc.cutoff * ((double)rd_crystal_order - 0.5);
	else
		cutoff = pbc.cutoff;

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {

					pair_ptr->rd_energy = 0;

					// pair LRC 
					if (rd_lrc)
						pair_ptr->lrc = lj_lrc_corr(atom_ptr, pair_ptr, cutoff);

					// to include a contribution, we require
					if ((pair_ptr->rimg - SMALL_dR < cutoff) &&   // inside cutoff?
						(!pair_ptr->rd_excluded || rd_crystal) &&   // either not excluded OR rd_crystal is ON
						(!pair_ptr->frozen)
						) { //not frozen

						//loop over unit cells
						if (rd_crystal) {
							sigma_over_r6 = 0;
							sigma_over_r12 = 0;
							for (i[0] = -(rd_crystal_order - 1); i[0] <= rd_crystal_order - 1; i[0]++)
								for (i[1] = -(rd_crystal_order - 1); i[1] <= rd_crystal_order - 1; i[1]++)
									for (i[2] = -(rd_crystal_order - 1); i[2] <= rd_crystal_order - 1; i[2]++) {
										if (!i[0] && !i[1] && !i[2] && pair_ptr->rd_excluded)
											continue; //no i=j=k=0 for excluded pairs (intra-molecular)
										//calculate pair separation (atom with it's image)
										for (int p = 0; p < 3; p++) {
											a[p] = 0;
											for (int q = 0; q < 3; q++)
												a[p] += pbc.basis[q][p] * i[q];
											a[p] += atom_ptr->pos[p] - pair_ptr->atom->pos[p];
										}
										r = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

										if (r > cutoff)
											continue;
										sigma_over_r = fabs(pair_ptr->sigma) / r;
										sigma_over_r6 += pow(sigma_over_r, 6);
										sigma_over_r12 += pow(sigma_over_r, 12);
									}
						}
						else { //otherwise, calculate as normal
							sigma_over_r = fabs(pair_ptr->sigma) / pair_ptr->rimg;
							sigma_over_r6 = sigma_over_r * sigma_over_r*sigma_over_r;
							sigma_over_r6 *= sigma_over_r6;
							sigma_over_r12 = sigma_over_r6 * sigma_over_r6;
						}

						// the LJ potential 
						if (spectre) {
							term6 = 0;
							term12 = sigma_over_r12;
							potential_classical = term12;
						}
						else {

							if (polarvdw)
								term6 = 0; //vdw calc'd by vdw.c	
							else term6 = sigma_over_r6;


							if (pair_ptr->attractive_only)
								term12 = 0;
							else
								term12 = sigma_over_r12;


							if (cdvdw_sig_repulsion)
								potential_classical = pair_ptr->sigrep*term12; //C6*sig^6/r^12
							else
								potential_classical = 4.0*pair_ptr->epsilon*(term12 - term6);
						}

						pair_ptr->rd_energy += potential_classical;

						if (feynman_hibbs)
							pair_ptr->rd_energy += lj_fh_corr(molecule_ptr, pair_ptr, feynman_hibbs_order, term12, term6);

						// if cavity_autoreject is on (cavity_autoreject_absolute is performed in energy.c)
						if (cavity_autoreject)
							if (pair_ptr->rimg < cavity_autoreject_scale*fabs(pair_ptr->sigma))
								pair_ptr->rd_energy = MAXVALUE;

					} //count contributions

				} // if recalculate

				// sum all of the pairwise terms 
				potential += pair_ptr->rd_energy + pair_ptr->lrc;

			} // pair
		} // atom
	} // molecule

	// molecule self-energy for rd_crystal -> energy of molecule interacting with its periodic neighbors

	if (rd_crystal)
		for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
			for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				potential += rd_crystal_self(atom_ptr, cutoff);

	// calculate self LRC interaction
	if (rd_lrc)
		for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
			for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				potential += lj_lrc_self(atom_ptr, cutoff);

	return potential;

}



double System::lj_lrc_corr(Atom * atom_ptr, Pair * pair_ptr, double cutoff)
{

	double sig_cut  = 0,
	       sig3     = 0,
	       sig_cut3 = 0,
	       sig_cut9 = 0;

	// include the long-range correction.    I'm  not sure that I'm handling spectre pairs correctly
	// we can't use rd_excluded flag, since that disqualifies inter-molecular, but that DOES contribute to LRC
	// ALL OF THESE MUST BE TRUE TO PERFORM LRC CALCULATION
	if ((pair_ptr->epsilon != 0 && pair_ptr->sigma != 0) &&  //if these are zero, then we won't waste our time
		!(atom_ptr->spectre && pair_ptr->atom->spectre) &&  //i think we want to disqualify s-s pairs 
		!(pair_ptr->frozen) &&  //disqualify frozen pairs
		((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != pbc.volume)
		) { //LRC only changes if the volume change

		pair_ptr->last_volume = pbc.volume;
		sig_cut = fabs(pair_ptr->sigma) / cutoff;
		sig3 = fabs(pair_ptr->sigma);
		sig3 *= sig3 * sig3;
		sig_cut3 = sig_cut * sig_cut*sig_cut;
		sig_cut9 = sig_cut3 * sig_cut3*sig_cut3;

		if (cdvdw_sig_repulsion)
			return (4.0 / 9.0) * pi * pair_ptr->sigrep * sig3 * sig_cut9 / pbc.volume;
		else if (polarvdw) //only repulsion term, if polarvdw is on
			return (16.0 / 9.0) * pi * pair_ptr->epsilon * sig3 * sig_cut9 / pbc.volume;
		else //if polarvdw is off, do the usual thing
			return ((16.0 / 3.0)*pi * pair_ptr->epsilon * sig3)  *  ((1.0 / 3.0)*sig_cut9 - sig_cut3) / pbc.volume;
	}
	else return pair_ptr->lrc; //use stored value

}


double System::lj_lrc_self(Atom * atom_ptr, double cutoff)
{
	double sig_cut, sig3, sig_cut3, sig_cut9;

	if (((atom_ptr->sigma != 0) && (atom_ptr->epsilon != 0)) &&  // non-zero parameters
		!(atom_ptr->frozen) &&  // not frozen
		!(atom_ptr->spectre)                                      // not spectre 
		) {

		sig_cut = fabs(atom_ptr->sigma) / cutoff;
		sig3 = fabs(atom_ptr->sigma);
		sig3 *= sig3 * sig3;
		sig_cut3 = sig_cut * sig_cut*sig_cut;
		sig_cut9 = sig_cut3 * sig_cut3*sig_cut3;

		if (cdvdw_sig_repulsion)
			return (1.0 / 3.0) * pi * hBar / kB * au2invseconds * atom_ptr->omega * atom_ptr->polarizability * atom_ptr->polarizability / sig3 * sig_cut9 / pbc.volume;
		else if (polarvdw) //only repulsion term, if polarvdw is on
			return (16.0 / 9.0) * pi * atom_ptr->epsilon * sig3 * sig_cut9 / pbc.volume;
		else //if polarvdw is off, do the usual thing
			return ((16.0 / 3.0) * pi * atom_ptr->epsilon*sig3)  *  ((1.0 / 3.0) * sig_cut9 - sig_cut3) / pbc.volume;
	}

	return 0;
}



double System::lj_fh_corr(Molecule * molecule_ptr, Pair * pair_ptr, int order, double term12, double term6)
{
	double reduced_mass;
	double dE, d2E, d3E, d4E; //energy derivatives
	double corr;
	double ir = 1.0 / pair_ptr->rimg;
	double ir2 = ir * ir;
	double ir3 = ir2 * ir;
	double ir4 = ir3 * ir;

	if ((order != 2) && (order != 4))
		throw invalid_setting; //must be order 2 or 4

	reduced_mass = AMU2KG * molecule_ptr->mass*pair_ptr->molecule->mass /
		(molecule_ptr->mass + pair_ptr->molecule->mass);

	if (cdvdw_sig_repulsion) {
		dE = -6.0*pair_ptr->sigrep*(2.0*term12 - term6) * ir;
		d2E = 6.0*pair_ptr->sigrep*(26.0*term12 - 7.0*term6) * ir2;
	}
	else {
		dE = -24.0*pair_ptr->epsilon*(2.0*term12 - term6) * ir;
		d2E = 24.0*pair_ptr->epsilon*(26.0*term12 - 7.0*term6) * ir2;
	}

	//2nd order correction
	corr = M2A2 *
		(hBar2 / (24.0*kB*temperature*reduced_mass)) *
		(d2E + 2.0*dE / pair_ptr->rimg);

	if (order >= 4) {

		if (cdvdw_sig_repulsion) {
			d3E = -336.0*pair_ptr->sigrep*(6.0*term12 - term6) * ir3;
			d4E = 3024.0*pair_ptr->sigrep*(10.0*term12 - term6) * ir4;
		}
		else {
			d3E = -1344.0*pair_ptr->epsilon*(6.0*term12 - term6) * ir3;
			d4E = 12096.0*pair_ptr->epsilon*(10.0*term12 - term6) * ir4;
		}

		//4th order corection
		corr += M2A4 *
			(hBar4 / (1152.0*kB2*temperature*temperature*reduced_mass*reduced_mass)) *
			(15.0*dE*ir3 + 4.0*d3E*ir + d4E);
	}

	return corr;
}



double System::rd_crystal_self(Atom * aptr, double cutoff)
{

	double curr_pot, term12, term6;
	double sigma_over_r6, sigma_over_r12, sigma_over_r, r;
	int i[3], p, q;
	double a[3];
	curr_pot = 0;

	if (aptr->sigma == 0 && aptr->epsilon == 0) return 0; //skip if no LJ interaction

	sigma_over_r6 = 0;
	sigma_over_r12 = 0; //need to init these guys

	for (i[0] = -(rd_crystal_order - 1); i[0] <= rd_crystal_order - 1; i[0]++)
		for (i[1] = -(rd_crystal_order - 1); i[1] <= rd_crystal_order - 1; i[1]++)
			for (i[2] = -(rd_crystal_order - 1); i[2] <= rd_crystal_order - 1; i[2]++) {
				if (!i[0] && !i[1] && !i[2]) continue; //no (0,0,0)
				//calculate pair separation (atom with it's image)
				for (p = 0; p < 3; p++) {
					a[p] = 0;
					for (q = 0; q < 3; q++)
						a[p] += pbc.basis[q][p] * i[q];
				}
				r = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

				if (r > cutoff)
					continue; //too far away! will be included in LRC if enabled
				sigma_over_r = fabs(aptr->sigma) / r;
				sigma_over_r6 += 0.5*pow(sigma_over_r, 6); //multiply by 0.5 to get counting correct
				sigma_over_r12 += 0.5*pow(sigma_over_r, 12);
			}

	if (spectre) {
		term6 = 0;
		curr_pot = term12 = sigma_over_r12;
	}
	else {
		if (polarvdw)
			term6 = 0; //vdw calc'd by vdw.c
		else
			term6 = sigma_over_r6;

		if (aptr->sigma < 0.0)
			term12 = 0; //attractive only
		else
			term12 = sigma_over_r12;

		if (cdvdw_sig_repulsion)
			curr_pot = 0.75*hBar / kB * au2invseconds*aptr->omega*aptr->polarizability*aptr->polarizability / pow(aptr->sigma, 6)*term12; //C6*sig^6/r^12
		else if (polarvdw)
			curr_pot = 4.0*aptr->epsilon*term12;
		else
			curr_pot = 4.0*aptr->epsilon*(term12 - term6);
	}
	return curr_pot;
}



double System::lj_buffered_14_7()
{
	double potential = 0.0, potential_classical;
	double r_over_sigma, first_term, second_term;

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	Pair     * pair_ptr = nullptr;

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {
					pair_ptr->rd_energy = 0;

					// make sure we're not excluded or beyond the cutoff
					if (!((pair_ptr->rimg > pbc.cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
						r_over_sigma = pair_ptr->rimg / pair_ptr->sigma;
						first_term = pow(1.07 / (r_over_sigma + 0.07), 7);
						second_term = (1.12 / (pow(r_over_sigma, 7) + 0.12) - 2);
						potential_classical = pair_ptr->epsilon*first_term*second_term;
						pair_ptr->rd_energy += potential_classical;

						// cavity autoreject
						if (cavity_autoreject)
							if (pair_ptr->rimg < cavity_autoreject_scale * pair_ptr->sigma)
								pair_ptr->rd_energy = MAXVALUE;
					}
				}
				potential += pair_ptr->rd_energy;
			}
		}
	}

	return potential;
}



double System::lj_buffered_14_7_nopbc()
{
	double potential = 0,
		potential_classical = 0,
		r_over_sigma = 0,
		first_term = 0,
		second_term = 0;

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	Pair     * pair_ptr = nullptr;

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {
					pair_ptr->rd_energy = 0;

					// make sure we're not excluded or beyond the cutoff
					if (!((pair_ptr->rimg > pbc.cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
						r_over_sigma = pair_ptr->rimg / pair_ptr->sigma;
						first_term = pow(1.07 / (r_over_sigma + 0.07), 7);
						second_term = (1.12 / (pow(r_over_sigma, 7) + 0.12) - 2);
						potential_classical = pair_ptr->epsilon*first_term*second_term;
						pair_ptr->rd_energy += potential_classical;

						// cavity autoreject
						if (cavity_autoreject)
							if (pair_ptr->rimg < cavity_autoreject_scale*pair_ptr->sigma)
								pair_ptr->rd_energy = MAXVALUE;
					}
				}
				potential += pair_ptr->rd_energy;
			}
		}
	}

	return potential;
}






//   Coulombic Energy Functions   ///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// No ewald summation - regular accumulation of Coulombic terms without out consideration of PBC
// Only used by surface module 
double System::coulombic_nopbc(Molecule * molecules) {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	double pe = 0;
	double total_pe = 0;

	total_pe = 0;
	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
				if (!pair_ptr->es_excluded) {
					pe = atom_ptr->charge * pair_ptr->atom->charge / pair_ptr->r;
					total_pe += pe;
				}
			}
		}
	}

	return total_pe;
}



double System::coulombic_nopbc_gwp() {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	double pe = 0,
		total_pe = 0;
	double qi = 0,
		qj = 0,
		ai = 0,
		aj = 0,
		r = 0;

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				r = pair_ptr->rimg;
				qi = atom_ptr->charge;
				qj = pair_ptr->atom->charge;
				ai = atom_ptr->gwp_alpha;
				aj = pair_ptr->atom->gwp_alpha;

				if (atom_ptr->gwp_spin || pair_ptr->atom->gwp_spin) {
					pe = qi * qj *  erf(sqrt(3.0 / 2.0*(ai*ai + aj * aj)) * r) / r;
				}
				else {
					pe = qi * qj / r;
				}

				total_pe += pe;
			}
		}
	}

	return total_pe;
}



double System::coulombic_kinetic_gwp() {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	double     ai = 0,
		mass = 0,
		energy = 0;

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			if (atom_ptr->gwp_spin) {
				ai = atom_ptr->gwp_alpha;
				mass = atom_ptr->mass;
				energy += 9.0 * hBar*hBar / (8.0* (ai / METER2ANGSTROM) * (ai / METER2ANGSTROM) * (AMU2KG*mass)) / kB;
			}

		} // atom
	} // molecule

	return(energy);
}


// total ES energy term 
double System::coulombic() {

	double real = 0,
		reciprocal = 0,
		self = 0,
		potential = 0;

	// construct the relevant ewald terms 
	if (wolf)
		potential = coulombic_wolf();
	else {
		real = coulombic_real();
		reciprocal = coulombic_reciprocal();
		self = coulombic_self();

		// return the total electrostatic energy
		potential = real + reciprocal + self;
	}

	return potential;
}



double System::coulombic_wolf() {

	Molecule * mptr;
	Atom     * aptr;
	Pair     * pptr;

	double pot = 0,
		alpha = ewald_alpha,
		R = pbc.cutoff,
		iR = 1.0 / R,
		erfaRoverR = erf(alpha*R) / R,
		r = 0,
		ir = 0;

	for (mptr = molecules; mptr; mptr = mptr->next) {
		for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
			for (pptr = aptr->pairs; pptr; pptr = pptr->next) {

				if (pptr->recalculate_energy) {
					pptr->es_real_energy = 0;

					r = pptr->rimg;
					ir = 1.0 / r;
					if ((!pptr->frozen) && (!pptr->es_excluded) && (r < R)) {
						pptr->es_real_energy =
							aptr->charge * pptr->atom->charge * (ir - erfaRoverR - iR * iR*(R - r));

						// get feynman-hibbs contribution
						if (feynman_hibbs) {
							Output::err("COULOMBIC: FH + es_wolf is not implemented\n");
							throw incompatible_settings;
						}
					}  // r<cutoff
				} //recalculate

				pot += pptr->es_real_energy;

			} //pair
		} //atom
	} //molecule

	return(pot);
}


// real space sum 
double System::coulombic_real() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	Pair     * pair_ptr = nullptr;

	double alpha = ewald_alpha,
		r = 0,
		erfc_term = 0,
		gaussian_term = 0,
		potential = 0,
		potential_classical = 0;


	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {
					pair_ptr->es_real_energy = 0;

					if (!pair_ptr->frozen) {

						r = pair_ptr->rimg;
						if (!((r > pbc.cutoff) || pair_ptr->es_excluded)) { // unit cell part

							//calculate potential contribution
							erfc_term = erfc(alpha*r);
							gaussian_term = exp(-alpha * alpha*r*r);
							potential_classical = atom_ptr->charge * pair_ptr->atom->charge * erfc_term / r;
							//store for pair pointer, so we don't always have to recalculate
							pair_ptr->es_real_energy += potential_classical;

							if (feynman_hibbs)
								pair_ptr->es_real_energy += coulombic_real_FH(molecule_ptr, pair_ptr, gaussian_term, erfc_term);

						}
						else if (pair_ptr->es_excluded) // calculate the charge-to-screen interaction
							pair_ptr->es_self_intra_energy = atom_ptr->charge * pair_ptr->atom->charge * erf(alpha*pair_ptr->r) / pair_ptr->r;

					} // frozen 
				} // recalculate 

				// sum all of the pairwise terms
				potential += pair_ptr->es_real_energy - pair_ptr->es_self_intra_energy;

			} // pair
		} // atom
	} // molecule

	return potential;
}


// feynman-hibbs for real space
double System::coulombic_real_FH(Molecule * molecule_ptr, Pair *pair_ptr, double gaussian_term, double erfc_term) {

	double du = 0,
		d2u = 0,
		d3u = 0,
		d4u = 0; //derivatives of the pair term
	double fh_2nd_order = 0;
	double fh_4th_order = 0;
	double r = pair_ptr->rimg;
	double rr = r * r;
	double ir = 1.0 / r;
	double ir2 = ir * ir;
	double ir3 = ir * ir2;
	double ir4 = ir2 * ir2;
	double order = feynman_hibbs_order;
	double alpha = ewald_alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * alpha;
	double a4 = a3 * alpha;
	double reduced_mass = AMU2KG * molecule_ptr->mass*pair_ptr->molecule->mass / (molecule_ptr->mass + pair_ptr->molecule->mass);

	du = -2.0 * alpha * gaussian_term / (r*sqrt(pi)) - erfc_term * ir2;
	d2u = (4.0 / sqrt(pi))*gaussian_term*(a3 + 1.0*ir2) + 2.0*erfc_term*ir3;

	fh_2nd_order = (M2A2)  *  (hBar2 / (24.0 * kB * temperature * reduced_mass))  *  (d2u + 2.0*du / r);

	if (order >= 4) {

		d3u = (gaussian_term / sqrt(pi))*(-8.0*(a3*a2)*r - 8.0*(a3) / r - 12.0*alpha*ir3) - 6.0*erfc(alpha*r)*ir4;
		d4u = (gaussian_term / sqrt(pi))*(8.0*a3*a2 + 16.0*a3*a4*rr + 32.0*a3*ir2 + 48.0*ir4) + 24.0*erfc_term*(ir4*ir);
		fh_4th_order = M2A4 * (hBar4 / (1152.0*(kB*kB*temperature*temperature*reduced_mass*reduced_mass)))  *  (15.0*du*ir3 + 4.0*d3u / r + d4u);
	}
	else
		fh_4th_order = 0.0;

	return fh_2nd_order + fh_4th_order;
}


// fourier space sum 
double System::coulombic_reciprocal() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;

	int    kmax             = ewald_kmax,
	       l[3]             = { 0 };
	double alpha            = ewald_alpha,
	       k[3]             = { 0 },
	       k_squared        = 0,
	       position_product = 0,
	       SF_re            = 0,
	       SF_im            = 0, // structure factor 
	       potential        = 0;

	// perform the fourier sum over a hemisphere (skipping certain points to avoid overcounting the face) 
	for (l[0] = 0; l[0] <= kmax; l[0]++) {
		for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
			for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

				//if norm is out of the sphere, skip
				if (UsefulMath::iidotprod(l, l) > kmax*kmax)
					continue;

				// get the reciprocal lattice vectors 
				for (int p = 0; p < 3; p++) {
					k[p] = 0;
					for (int q = 0; q < 3; q++)
						k[p] += 2.0 * pi * pbc.reciprocal_basis[p][q] * l[q];
				}
				k_squared = UsefulMath::dddotprod(k, k);

				// structure factor 
				SF_re = 0;
				SF_im = 0;
				for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
					for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

						if (atom_ptr->frozen)
							continue; //skip frozen
						if (atom_ptr->charge == 0.0)
							continue; //skip if no charge

						// the inner product of the position vector and the k vector 
						position_product = UsefulMath::dddotprod(k, atom_ptr->pos);

						SF_re += atom_ptr->charge * cos(position_product);
						SF_im += atom_ptr->charge * sin(position_product);

					} // atom
				} // molecule

				potential += exp(-k_squared / (4.0*alpha*alpha)) / k_squared * (SF_re*SF_re + SF_im * SF_im);

			} // end for n
		} // end for m
	} // end for l

	potential *= 4.0 * pi / pbc.volume;

	return potential;
}



double System::coulombic_self() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	double     alpha = ewald_alpha,
		self_potential = 0.0;

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			if (atom_ptr->frozen)
				continue;
			atom_ptr->es_self_point_energy = alpha * atom_ptr->charge*atom_ptr->charge / sqrt(pi);
			self_potential -= atom_ptr->es_self_point_energy;
		}
	}

	return self_potential;
}





//   Axilrod Teller 3-body energy   /////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double System::axilrod_teller ()
// Copyright 2015 Adam Hogan
{
	double potential = 0.0, c9, rij, rik, rjk, cos_part;
	Molecule * molecule1, *molecule2, *molecule3;
	Atom     * atom1, *atom2, *atom3;
	Pair temp;
	Vector3D  ij, ik, jk, a, b;

	for (molecule1 = molecules; molecule1; molecule1 = molecule1->next) {
		for (molecule2 = molecules; molecule2; molecule2 = molecule2->next) {
			for (molecule3 = molecules; molecule3; molecule3 = molecule3->next) {
				int number_of_unique_molecules = 1;
				if (molecule2 != molecule1) {
					number_of_unique_molecules++;
				}
				if (molecule3 != molecule1 && molecule3 != molecule2) {
					number_of_unique_molecules++;
				}
				if (number_of_unique_molecules > 1) {
					for (atom1 = molecule1->atoms; atom1; atom1 = atom1->next) {
						for (atom2 = molecule2->atoms; atom2; atom2 = atom2->next) {
							for (atom3 = molecule3->atoms; atom3; atom3 = atom3->next) {
								if (atom1 != atom2 && atom1 != atom3 && atom2 != atom3) {

									double atom1_c9, atom2_c9, atom3_c9;
									atom1_c9 = atom1->c9;
									atom2_c9 = atom2->c9;
									atom3_c9 = atom3->c9;

									// Axilrod-Teller parameters in http://arxiv.org/pdf/1201.1532.pdf and http://dx.doi.org/10.1063/1.440310
									// Midzuno-Kihara approximation for c9 http://dx.doi.org/10.1143/JPSJ.11.1045
									if (midzuno_kihara_approx)
									{
										atom1_c9 = 3.0 / 4.0*atom1->polarizability*6.7483345*atom1->c6;
										atom2_c9 = 3.0 / 4.0*atom2->polarizability*6.7483345*atom2->c6;
										atom3_c9 = 3.0 / 4.0*atom3->polarizability*6.7483345*atom3->c6;
									}

									// Mixing rule, ep. 20 of http://dx.doi.org/10.1063/1.440310
									c9 = pow(pow(atom1->polarizability*6.7483345, 3)*pow(atom2->polarizability*6.7483345, 3)*pow(atom3->polarizability*6.7483345, 3), 1.0 / 3.0) *
										3.0 / (1.0 / (atom1_c9 / pow(atom1->polarizability*6.7483345, 3)) + 1.0 / (atom2_c9 / pow(atom2->polarizability*6.7483345, 3)) + 1.0 / (atom3_c9 / pow(atom3->polarizability*6.7483345, 3)));

									if (atom1->polarizability == 0.0 || atom2->polarizability == 0.0 || atom3->polarizability == 0.0)
										c9 = 0.0; // avoid division by zero

									c9 *= 0.0032539449 / (3.166811429*0.000001); // convert H*Bohr^9 to K*Angstrom^9

									// reset fake pair ptr
									temp.d_prev[0] = temp.d_prev[1] = temp.d_prev[2] = -999999999999.;
									// get minimum image distance
									minimum_image (atom1, atom2, &temp);
									rij = temp.rimg;
									ij.set (temp.dimg[0], temp.dimg[1], temp.dimg[2]);

									// reset fake pair ptr
									temp.d_prev[0] = temp.d_prev[1] = temp.d_prev[2] = -999999999999.;
									// get minimum image distance
									minimum_image (atom1, atom3, &temp);
									rik = temp.rimg;
									ik.set (temp.dimg[0], temp.dimg[1], temp.dimg[2]);

									// reset fake pair ptr
									temp.d_prev[0] = temp.d_prev[1] = temp.d_prev[2] = -999999999999.;
									// get minimum image distance
									minimum_image (atom2, atom3, &temp);
									rjk = temp.rimg;
									jk.set (temp.dimg[0], temp.dimg[1], temp.dimg[2]);

									cos_part = 3;

									a = -1.0 * ij;
									b = -1.0 * ik;

									cos_part *= (a * b) / (a.norm() * b.norm());

									a = ij;
									b = -1.0*jk;

									cos_part *= (a * b) / (a.norm() * b.norm());

									a = ik;
									b = jk;

									cos_part *= (a * b) / (a.norm() * b.norm());

									potential += c9 * ((1.0 + cos_part) / pow(rij*rik*rjk, 3));
								}
							}
						}
					}
				}
			}
		}
	}
	// We're counting each pair 6 times
	potential = potential / 6;
	return potential;
}






//   Silvera-Goldman H2 potential   /////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Silvera-Goldman parameters
// http://www.pnas.org/content/99/3/1129.full.pdf
const double ALPHA = 1.713;    // unitless
const double BETA = 1.5671;    // 1/a.u.
const double GAMMA = 0.00993;  // 1/a.u.^2
const double  C6 = 12.14;      // multipole term1 a.u.^6
const double  C8 = 215.2;      // multipole term2 a.u.^8
const double  C10 = 4813.9;    // multipole term3 a.u.^10
const double  C9 = 143.1;      // 3-body term a.u.^9
const double  RM = 8.321;      // position of max well depth (a.u.) times 1.28


double System::sg() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	Pair     * pair_ptr = nullptr;
	double     rimg = 0,
		r6 = 0,
		r8 = 0,
		r9 = 0,
		r10 = 0,
		r_rm = 0,
		repulsive_term = 0,
		multipole_term = 0,
		exponential_term = 0,
		first_r_diff_term = 0,
		second_r_diff_term = 0,
		first_derivative = 0,
		second_derivative = 0,
		potential_classical = 0,
		potential_fh_second_order = 0,
		potential = 0;


	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {

					pair_ptr->rd_energy = 0;
					rimg = pair_ptr->rimg;

					if (rimg < pbc.cutoff) {

						// convert units to Bohr radii 
						rimg /= AU2ANGSTROM;

						// classical pairwise part 
						repulsive_term = exp(ALPHA - BETA * rimg - GAMMA * rimg*rimg);

						r6 = pow(rimg, 6);
						r8 = pow(rimg, 8);
						r9 = pow(rimg, 9);
						r10 = pow(rimg, 10);
						multipole_term = C6 / r6 + C8 / r8 + C10 / r10 - C9 / r9;


						r_rm = RM / rimg;
						if (rimg < RM)
							exponential_term = exp(-pow((r_rm - 1.0), 2));
						else
							exponential_term = 1.0;

						potential_classical = (repulsive_term - multipole_term * exponential_term);
						pair_ptr->rd_energy += potential_classical;

						if (feynman_hibbs) {

							// FIRST DERIVATIVE 
							first_derivative = (-BETA - 2.0*GAMMA*rimg)*repulsive_term;
							first_derivative += (6.0*C6 / pow(rimg, 7) + 8.0*C8 / pow(rimg, 9) - 9.0*C9 / pow(rimg, 10) + 10.0*C10 / pow(rimg, 11))*exponential_term;
							first_r_diff_term = (r_rm*r_rm - r_rm) / rimg;
							first_derivative += -2.0*multipole_term*exponential_term*first_r_diff_term;

							// SECOND DERIVATIVE
							second_derivative = (pow((BETA + 2.0*GAMMA*rimg), 2) - 2.0*GAMMA)*repulsive_term;
							second_derivative += (-exponential_term)*(42.0*C6 / pow(rimg, 8) + 72.0*C8 / pow(rimg, 10) - 90.0*C9 / pow(rimg, 11) + 110.0*C10 / pow(rimg, 10));
							second_derivative += exponential_term * first_r_diff_term*(12.0*C6 / pow(rimg, 7) + 16.0*C8 / pow(rimg, 9) - 18.0*C9 / pow(rimg, 10) + 20.0*C10 / pow(rimg, 11));
							second_derivative += exponential_term * pow(first_r_diff_term, 2)*4.0*multipole_term;
							second_r_diff_term = (3.0*r_rm*r_rm - 2.0*r_rm) / (rimg*rimg);
							second_derivative += exponential_term * second_r_diff_term*2.0*multipole_term;

							potential_fh_second_order = pow(METER2ANGSTROM, 2)*(hBar*hBar / (24.0*kB*temperature*(AMU2KG*molecule_ptr->mass)))*(second_derivative + 2.0*first_derivative / rimg);
							pair_ptr->rd_energy += potential_fh_second_order;
						}

						// convert units from Hartrees back to Kelvin 
						pair_ptr->rd_energy *= HARTREE2KELVIN;

					}

				} // recalculate
			} // pair
		} // atom
	} // molecule

	potential = 0;
	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
				potential += pair_ptr->rd_energy;

	return potential;

}

// same as above, but no periodic boundary conditions 
double System::sg_nopbc(Molecule * molecules) {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	Pair     * pair_ptr = nullptr;
	double     potential = 0,
		multipole_term = 0,
		result = 0,
		exp_result = 0,
		r = 0, r6 = 0, r8 = 0, r10 = 0, r9 = 0, r_rm = 0, r_rm_2 = 0, r_exp = 0;


	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {

					r = pair_ptr->r / AU2ANGSTROM;


					r6 = pow(r, 6);
					r8 = pow(r, 8);
					r10 = pow(r, 10);
					r9 = pow(r, 9);

					multipole_term = C6 / r6 + C8 / r8 + C10 / r10 - C9 / r9;

					if (r < RM) {
						r_rm = RM / r;
						r_rm -= 1.0;
						r_rm_2 = pow(r_rm, 2);
						r_rm_2 *= -1.0;
						r_exp = exp(r_rm_2);

						multipole_term *= r_exp;

					}

					result = ALPHA - BETA * r - GAMMA * r*r;

					exp_result = exp(result);
					exp_result -= multipole_term;
					pair_ptr->rd_energy = HARTREE2KELVIN * exp_result;

				} // recalculate 
			} // pair 
		} // atom 
	} // molecule 

	potential = 0;
	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
				potential += pair_ptr->rd_energy;

	return potential;

}






//   Dispersion/Expansion   /////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double System::disp_expansion()
{
	double potential = 0.0;

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {

					// pair LRC
					if (rd_lrc)
						pair_ptr->lrc = disp_expansion_lrc(pair_ptr, pbc.cutoff);

					// make sure we're not excluded or beyond the cutoff
					if (!(pair_ptr->rd_excluded || pair_ptr->frozen)) {
						const double r = pair_ptr->rimg;
						const double r2 = r * r;
						const double r4 = r2 * r2;
						const double r6 = r4 * r2;
						const double r8 = r6 * r2;
						const double r10 = r8 * r2;

						double c6 = pair_ptr->c6;
						const double c8 = pair_ptr->c8;
						const double c10 = pair_ptr->c10;

						if (disp_expansion_mbvdw == 1)
							c6 = 0.0;

						double repulsion = 0.0;

						if (pair_ptr->epsilon != 0.0    &&    pair_ptr->sigma != 0.0)
							repulsion = 315.7750382111558307123944638 * exp(-pair_ptr->epsilon*(r - pair_ptr->sigma)); // K = 10^-3 H ~= 316 K

						if (damp_dispersion)
							pair_ptr->rd_energy = -tt_damping(6, pair_ptr->epsilon*r)*c6 / r6 - tt_damping(8, pair_ptr->epsilon*r)*c8 / r8 - tt_damping(10, pair_ptr->epsilon*r)*c10 / r10 + repulsion;
						else
							pair_ptr->rd_energy = -c6 / r6 - c8 / r8 - c10 / r10 + repulsion;

						if (cavity_autoreject)
						{
							if (r < cavity_autoreject_scale * pair_ptr->sigma)
								pair_ptr->rd_energy = MAXVALUE;
							if (cavity_autoreject_repulsion != 0.0   &&   repulsion > cavity_autoreject_repulsion)
								pair_ptr->rd_energy = MAXVALUE;
						}
					}

				}
				potential += pair_ptr->rd_energy + pair_ptr->lrc;
			}
		}
	}

	if (disp_expansion_mbvdw == 1)
	{
		thole_amatrix();
		potential += vdw();
	}

	// calculate self LRC interaction 
	if (rd_lrc)
	{
		for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		{
			for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			{
				atom_ptr->lrc_self = disp_expansion_lrc_self(atom_ptr, pbc.cutoff);
				potential += atom_ptr->lrc_self;
			}
		}
	}

	return potential;
}



double System::disp_expansion_lrc(Pair * pair_ptr, const double cutoff) // ignoring the exponential repulsion bit because it decays exponentially
{
	if (!(pair_ptr->frozen) &&  // disqualify frozen pairs
		((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != pbc.volume)) { // LRC only changes if the volume change

		pair_ptr->last_volume = pbc.volume;

		return -4.0 * pi * (pair_ptr->c6 / (3.0*cutoff*cutoff*cutoff) + pair_ptr->c8 / (5.0*cutoff*cutoff*cutoff*cutoff*cutoff) + pair_ptr->c10 / (7.0*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff)) / pbc.volume;
	}

	else return pair_ptr->lrc; // use stored value
}



double System::tt_damping(int n, double br)
{
	double sum = 0.0;

	for (int i = 0; i <= n; i++)
	{
		sum += pow(br, i) / UsefulMath::factorial(i);
	}

	const double result = 1.0 - exp(-br)*sum;

	if (result > 0.000000001)
		return result;
	else
		return 0.0; // This is so close to zero lets just call it zero to avoid rounding error and the simulation blowing up
}



double System::disp_expansion_lrc_self(Atom * atom_ptr, const double cutoff)
{
	if (!(atom_ptr->frozen) && // disqualify frozen atoms 
		((atom_ptr->lrc_self == 0.0) || atom_ptr->last_volume != pbc.volume)) { // LRC only changes if the volume change) 

		atom_ptr->last_volume = pbc.volume;

		if (extrapolate_disp_coeffs)
		{
			double c10;
			if (atom_ptr->c6 != 0.0&&atom_ptr->c8 != 0.0)
				c10 = 49.0 / 40.0*atom_ptr->c8*atom_ptr->c8 / atom_ptr->c6;
			else
				c10 = 0.0;

			return -4.0 * pi * (atom_ptr->c6 / (3.0*cutoff*cutoff*cutoff) + atom_ptr->c8 / (5.0*cutoff*cutoff*cutoff*cutoff*cutoff) + c10 / (7.0*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff)) / pbc.volume;
		}
		else
			return -4.0 * pi * (atom_ptr->c6 / (3.0*cutoff*cutoff*cutoff) + atom_ptr->c8 / (5.0*cutoff*cutoff*cutoff*cutoff*cutoff) + atom_ptr->c10 / (7.0*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff)) / pbc.volume;
	}

	return atom_ptr->lrc_self; // use stored value
}






//   Dreiding Potential   ///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//  @2010, Jonathan Belof
//  Space Research Group
//  Department of Chemistry
//  University of South Florida

const double DREIDING_GAMMA = 12.0;



double System::dreiding() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	Pair     * pair_ptr = nullptr;
	double     gamma = DREIDING_GAMMA,
		r_over_sigma = 0,
		termexp = 0,
		term6 = 0,
		potential = 0,
		potential_classical = 0;

#ifdef XXX
	double first_derivative, second_derivative, third_derivative, fourth_derivative;
	double potential, potential_classical, potential_fh_second_order, potential_fh_fourth_order;
	double reduced_mass;
	double sig3, sig_cut, sig_cut3, sig_cut9;
#endif // XXX 


	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {

					pair_ptr->rd_energy = 0;

					// make sure we're not excluded or beyond the cutoff 
					if (!((pair_ptr->rimg > pbc.cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {

						r_over_sigma = pair_ptr->rimg / pair_ptr->sigma;

						// the DREIDING potential
						term6 = pow(r_over_sigma, -6);
						term6 *= gamma / (gamma - 6.0);

						if (pair_ptr->attractive_only)
							termexp = 0;
						else {
							if (pair_ptr->rimg < 0.4*pair_ptr->sigma)
								termexp = MAXVALUE;
							else {
								termexp = exp(gamma*(1.0 - r_over_sigma));
								termexp *= (6.0 / (gamma - 6.0));
							}
						}
						potential_classical = pair_ptr->epsilon*(termexp - term6);

						pair_ptr->rd_energy += potential_classical;

						// XXX 
						// need to do fh for dreiding pbc->volume
#ifdef XXX
						if (system->feynman_hibbs) {

							reduced_mass = AMU2KG * molecule_ptr->mass*pair_ptr->molecule->mass / (molecule_ptr->mass + pair_ptr->molecule->mass);

							// FIRST DERIVATIVE
							first_derivative = -24.0*pair_ptr->epsilon*(2.0*term12 - term6) / pair_ptr->rimg;

							// SECOND DERIVATIVE
							second_derivative = 24.0*pair_ptr->epsilon*(26.0*term12 - 7.0*term6) / pow(pair_ptr->rimg, 2);

							potential_fh_second_order = pow(METER2ANGSTROM, 2)*(HBAR*HBAR / (24.0*KB*system->temperature*reduced_mass))*(second_derivative + 2.0*first_derivative / pair_ptr->rimg);
							pair_ptr->rd_energy += potential_fh_second_order;

							if (system->feynman_hibbs_order >= 4) {

								// THIRD DERIVATIVE
								third_derivative = -1344.0*pair_ptr->epsilon*(6.0*term12 - term6) / pow(pair_ptr->rimg, 3);

								// FOURTH DERIVATIVE
								fourth_derivative = 12096.0*pair_ptr->epsilon*(10.0*term12 - term6) / pow(pair_ptr->rimg, 4);

								potential_fh_fourth_order = pow(METER2ANGSTROM, 4)*(pow(HBAR, 4) / (1152.0*pow(KB*system->temperature*reduced_mass, 2)))*(15.0*first_derivative / pow(pair_ptr->rimg, 3) + 4.0*third_derivative / pair_ptr->rimg + fourth_derivative);
								pair_ptr->rd_energy += potential_fh_fourth_order;

							}

						}

#endif // XXX 
						// cause an autoreject on insertions closer than a certain amount 
						if (cavity_autoreject) {
							if (pair_ptr->rimg < cavity_autoreject_scale * pair_ptr->sigma)
								pair_ptr->rd_energy = MAXVALUE;
						}

					}

					// XXX need to derive lrc for dreiding
#ifdef XXX
	// include the long-range correction
					if (!(pair_ptr->rd_excluded || pair_ptr->frozen) && (pair_ptr->lrc == 0.0) && system->rd_lrc) {

						sig_cut = fabs(pair_ptr->sigma) / system->pbc->cutoff;
						sig3 = pow(fabs(pair_ptr->sigma), 3);
						sig_cut3 = pow(sig_cut, 3);
						sig_cut9 = pow(sig_cut, 9);
						pair_ptr->lrc = ((-8.0 / 3.0)*M_PI*pair_ptr->epsilon*sig3)* sig_cut3 / system->pbc->volume;

					}
#endif // XXX

				} // if recalculate

				// sum all of the pairwise terms
				potential += pair_ptr->rd_energy + pair_ptr->lrc;

			} // pair
		} // atom
	} // molecule


	return(potential);

}


// same as above, but no periodic boundary conditions 
double System::dreiding_nopbc(Molecule *molecules) {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	Pair     * pair_ptr = nullptr;
	double     gamma = DREIDING_GAMMA,
		r_over_sigma = 0,
		termexp = 0,
		term6 = 0,
		potential = 0;



	for (molecule_ptr = molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				// make sure we're not excluded or beyond the cutoff 
				if (!pair_ptr->rd_excluded) {

					r_over_sigma = pair_ptr->r / pair_ptr->sigma;

					// the DREIDING potential 
					term6 = pow(r_over_sigma, -6);
					term6 *= gamma / (gamma - 6.0);

					if (pair_ptr->attractive_only)
						termexp = 0;
					else {
						if (pair_ptr->rimg < 0.35*pair_ptr->sigma)
							termexp = MAXVALUE;
						else {
							termexp = exp(gamma*(1.0 - r_over_sigma));
							termexp *= (6.0 / (gamma - 6.0));
						}
					}
					potential += pair_ptr->epsilon*(termexp - term6);

				}

			} // pair
		} // atom
	} // molecule

	return potential;

}






//   EXP Repulsion   ////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double System::exp_repulsion()
{

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;
	double     r = 0,
		term = 0,
		potential = 0,
		potential_classical = 0,
		cutoff = 0;
	int        i[3] = { 0 };
	double     a[3] = { 0 };

	//set the cutoff
	if (rd_crystal)
		cutoff = 2.0 * pbc.cutoff * ((double)rd_crystal_order - 0.5);
	else
		cutoff = pbc.cutoff;

	potential = 0;
	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->recalculate_energy) {

					pair_ptr->rd_energy = 0;

					// pair LRC 
					if (rd_lrc) pair_ptr->lrc = exp_lrc_corr(atom_ptr, pair_ptr, cutoff);

					// to include a contribution, we require
					if ((pair_ptr->rimg - SMALL_dR < cutoff)  //inside cutoff?
						&& (!pair_ptr->rd_excluded || rd_crystal) //either not excluded OR rd_crystal is ON
						&& !pair_ptr->frozen) //not frozen
					{

						//loop over unit cells
						if (rd_crystal) {
							term = 0;
							for (i[0] = -(rd_crystal_order); i[0] <= rd_crystal_order; i[0]++)
								for (i[1] = -(rd_crystal_order); i[1] <= rd_crystal_order; i[1]++)
									for (i[2] = -(rd_crystal_order); i[2] <= rd_crystal_order; i[2]++) {
										if (!i[0] && !i[1] && !i[2] && pair_ptr->rd_excluded)
											continue; //no i=j=k=0 for excluded pairs (intra-molecular)
										//calculate pair separation (atom with it's image)
										for (int p = 0; p < 3; p++) {
											a[p] = 0;
											for (int q = 0; q < 3; q++)
												a[p] += pbc.basis[q][p] * i[q];
											a[p] += atom_ptr->pos[p] - pair_ptr->atom->pos[p];
										}
										r = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

										if (r + SMALL_dR > cutoff)	continue;
										term += exp(-r / (2.0*pair_ptr->epsilon));
									}
						}
						else //otherwise, calculate as normal
							term = exp(-pair_ptr->rimg / (2.0*pair_ptr->epsilon));

						potential_classical = pair_ptr->sigma*term;
						pair_ptr->rd_energy += potential_classical;

						if (feynman_hibbs)
							pair_ptr->rd_energy +=
							exp_fh_corr(molecule_ptr, pair_ptr, feynman_hibbs_order, potential_classical);

					} //count contributions

				} // if recalculate

				// sum all of the pairwise terms
				potential += pair_ptr->rd_energy + pair_ptr->lrc;

			} // pair 
		} // atom
	} // molecule

	// molecule self-energy for rd_crystal -> energy of molecule interacting with its periodic neighbors 

	if (rd_crystal)
		for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
			for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				potential += exp_crystal_self(atom_ptr, cutoff);

	// calculate self LRC interaction
	if (rd_lrc)
		for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
			for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				potential += exp_lrc_self(atom_ptr, cutoff);

	return(potential);

}



double System::exp_lrc_corr(Atom * atom_ptr, Pair * pair_ptr, double cutoff)
{

	double eps = pair_ptr->epsilon;
	double rover2e = cutoff / (2.0*eps);

	// include the long-range correction 
	// I'm  not sure that I'm handling spectre pairs correctly
	// we can't use rd_excluded flag, since that disqualifies inter-molecular, but that DOES contribute to LRC
	// ALL OF THESE MUST BE TRUE TO PERFORM LRC CALCULATION 
	if ((pair_ptr->epsilon != 0 && pair_ptr->sigma != 0) &&  //if these are zero, then we won't waste our time
		!(atom_ptr->spectre && pair_ptr->atom->spectre) && //i think we want to disqualify s-s pairs 
		!(pair_ptr->frozen) &&  //disqualify frozen pairs
		((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != pbc.volume)) { //LRC only changes if the volume change

		pair_ptr->last_volume = pbc.volume;

		return (8.0*pi)*exp(1. - rover2e)*(cutoff*cutoff + 4.0*eps*cutoff + 8.0*eps*eps)*pair_ptr->sigma / pbc.volume;

	}
	else return pair_ptr->lrc; //use stored value

}



double System::exp_fh_corr(Molecule * molecule_ptr, Pair * pair_ptr, int order, double pot)
{
	double reduced_mass = 0,
		dE = 0,
		d2E = 0,
		d3E = 0,
		d4E = 0, //energy derivatives
		corr = 0,
		ir = 1.0 / pair_ptr->rimg,
		ir2 = ir * ir,
		ir3 = ir2 * ir;

	if ((order != 2) && (order != 4))
		throw invalid_setting; //must be order 2 or 4

	reduced_mass = AMU2KG * molecule_ptr->mass*pair_ptr->molecule->mass / (molecule_ptr->mass + pair_ptr->molecule->mass);

	dE = -pot / (2.0*pair_ptr->epsilon);
	d2E = dE / (2.0*pair_ptr->epsilon);

	//2nd order correction
	corr = M2A2 *
		(hBar2 / (24.0*kB*temperature*reduced_mass)) *
		(d2E + 2.0*dE / pair_ptr->rimg);

	if (order >= 4) {

		d3E = -d2E / (2.0*pair_ptr->epsilon);
		d4E = d3E / (2.0*pair_ptr->epsilon);

		//4th order corection
		corr += M2A4 *
			(hBar4 / (1152.0*kB2*temperature*temperature*reduced_mass*reduced_mass)) *
			(15.0*dE*ir3 + 4.0*d3E*ir + d4E);
	}

	return corr;
}



double System::exp_crystal_self(Atom * aptr, double cutoff)
{
	double term = 0,
		r = 0,
		eps = aptr->epsilon;
	int    i[3] = { 0 };
	double a[3] = { 0 };

	if (aptr->sigma == 0 || aptr->epsilon == 0) return 0; //skip if no LJ interaction

	for (i[0] = -(rd_crystal_order); i[0] <= rd_crystal_order; i[0]++)
		for (i[1] = -(rd_crystal_order); i[1] <= rd_crystal_order; i[1]++)
			for (i[2] = -(rd_crystal_order); i[2] <= rd_crystal_order; i[2]++) {
				if (!i[0] && !i[1] && !i[2])
					continue; //no (0,0,0)
					//calculate pair separation (atom with it's image)
				for (int p = 0; p < 3; p++) {
					a[p] = 0;
					for (int q = 0; q < 3; q++)
						a[p] += pbc.basis[q][p] * i[q];
				}
				r = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

				if (r > cutoff) continue; //too far away! will be included in LRC if enabled
				term += 0.5*exp(-r / (2.0*eps));  //multiply by 0.5 to get counting correct
			}

	return aptr->sigma * term;
}



double System::exp_lrc_self(Atom * atom_ptr, double cutoff)
{
	double eps = atom_ptr->epsilon,
		rover2e = cutoff / (2.0*eps);

	if (((atom_ptr->sigma != 0) && (atom_ptr->epsilon != 0)) && //non-zero parameters
		!(atom_ptr->frozen) && //not frozen
		!(atom_ptr->spectre)) { //not spectre 

		return (8.0*pi)*exp(1. - rover2e)*(cutoff*cutoff + 4.0*eps*cutoff + 8.0*eps*eps)*atom_ptr->sigma / pbc.volume;
	}
	return 0;
}






//   Polarization Energy   //////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void print_matrix(int N, double **matrix) {

	int i, j;

	printf("\n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			printf("%.3f ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}



void zero_out(Molecule * m) {

	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	int p;

	//zero out the electric field for each site
	for (mptr = m; mptr; mptr = mptr->next) {
		for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
			for (p = 0; p < 3; p++) {
				aptr->ef_static[p] = 0.0;
				aptr->ef_static_self[p] = 0.0;
			}
		}
	}

	return;
}


// get the induction energy 
double System::polar() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;

	int        num_iterations = 0;
	double     potential = 0;
	char       linebuf[maxLine] = { 0 };

	// take measures to let N fluctuate
	if ((ensemble == ENSEMBLE_UVT || ensemble == ENSEMBLE_REPLAY) && !polar_zodid)
		thole_resize_matrices();

	// get the A matrix
	if (!polar_zodid) {
		thole_amatrix();
		if (polarizability_tensor) {
			Output::out("POLAR: A matrix:\n");
			print_matrix(3 * ((int)checkpoint->thole_N_atom), A_matrix);
		}
	}

	// find the dipoles

	if (polar_ewald_full) {
		//do a full-ewald polarization treatment
		ewald_full();

	}
	else if (polar_iterative) {
		//solve the self-consistent problem
		thole_field(); //calc e-field
		num_iterations = thole_iterative(); //calc dipoles

		nodestats->polarization_iterations = (double)num_iterations; //statistics
		observables->dipole_rrms = get_dipole_rrms();

		if (iterator_failed) {
			switch (ensemble) {
			case ENSEMBLE_UVT:
			case ENSEMBLE_NVT:
			case ENSEMBLE_NVE:
			case ENSEMBLE_NPT:
				sprintf(linebuf, "POLAR: polarization iterative solver convergence failure on mc step %d.\n", step);
				Output::err(linebuf);
				break;
			case ENSEMBLE_REPLAY:
				sprintf(linebuf, "POLAR: polarization iterative solver convergence failure on configuration %d.\n", step);
				Output::err(linebuf);
				break;
			case ENSEMBLE_SURF:
			case ENSEMBLE_SURF_FIT:
			case ENSEMBLE_TE:
			default:
				sprintf(linebuf, "POLAR: polarization iterative solver failed to reach convergence.\n");
				Output::err(linebuf);
			}
		}

	}
	else {
		//do matrix inversion
		thole_field(); //calc e-field
		thole_bmatrix(); //matrix inversion
		thole_bmatrix_dipoles(); //get dipoles

		// output the 3x3 molecular polarizability tensor 
		if (polarizability_tensor) {
			Output::out("POLAR: B matrix:\n");
			print_matrix(3 * ((int)checkpoint->thole_N_atom), B_matrix);
			thole_polarizability_tensor();
			throw exception_ok;
		}
	}

	// calculate the polarization energy as 1/2 mu*E
	potential = 0;
	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			potential += UsefulMath::dddotprod(atom_ptr->mu, atom_ptr->ef_static);
			if (polar_palmo)
				potential += UsefulMath::dddotprod(atom_ptr->mu, atom_ptr->ef_induced_change);
		}
	}
	potential *= -0.5;

#ifdef DEBUG
	fprintf(stderr, "mu MOLECULE ATOM * DIPOLES * STATIC * INDUCED * pot/atom -0.5*mu*E_s\n");
	for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			fprintf(stderr, "mu %4d %4d * %8.5lf %8.5lf %8.5lf * %8.5lf %8.5lf %8.5lf * %8.5lf %8.5lf %8.5lf * %lf %lf\n",
				molecule_ptr->id, atom_ptr->id,
				atom_ptr->mu[0], atom_ptr->mu[1], atom_ptr->mu[2],
				atom_ptr->ef_static[0], atom_ptr->ef_static[1], atom_ptr->ef_static[2],
				atom_ptr->ef_induced[0], atom_ptr->ef_induced[1], atom_ptr->ef_induced[2],
				potential / system->natoms, -0.5*atom_ptr->mu[0] * atom_ptr->ef_static[0]);
		}
	}
#endif

	return potential;
}


// RRMS of dipoles
double System::get_dipole_rrms() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;

	double     NAtoms = 0,
		dipole_rrms = 0;



	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			if (std::isfinite(atom_ptr->dipole_rrms))
				dipole_rrms += atom_ptr->dipole_rrms;
			NAtoms++;
		}
	}
	return dipole_rrms / NAtoms;
}


// calculate the dipole field tensor 
void System::thole_amatrix() {

	int     ii, jj;
	int     NAtoms = natoms;
	Pair  * pair_ptr = nullptr;
	double  damp1 = 0, damp2 = 0, wdamp1 = 0, wdamp2 = 0, v = 0, s = 0,
		r = 0, r2 = 0, ir3 = 0, ir5 = 0, ir = 0,

		rcut = pbc.cutoff,
		rcut2 = rcut * rcut,
		rcut3 = rcut2 * rcut,

		l = polar_damp,
		l2 = l * l,
		l3 = l2 * l,

		explr = 0, //exp(-l*r)
		explrcut = exp(-l * rcut);

	zero_out_amatrix(NAtoms);

	// set the diagonal blocks 
	for (int i = 0; i < NAtoms; i++) {
		ii = i * 3;
		for (int p = 0; p < 3; p++) {
			if (atom_array[i]->polarizability != 0.0)
				A_matrix[ii + p][ii + p] = 1.0 / atom_array[i]->polarizability;
			else
				A_matrix[ii + p][ii + p] = MAXVALUE;
		}
	}

	// calculate each Tij tensor component for each dipole pair
	for (int i = 0; i < (NAtoms - 1); i++) {
		ii = i * 3;
		pair_ptr = atom_array[i]->pairs;
		for (int j = (i + 1); j < NAtoms; j++, pair_ptr = pair_ptr->next) {
			jj = j * 3;

			r = pair_ptr->rimg;
			r2 = r * r;

			// inverse displacements
			if (pair_ptr->rimg == 0.)
				ir3 = ir5 = MAXVALUE;
			else {
				ir = 1.0 / r;
				ir3 = ir * ir*ir;
				ir5 = ir3 * ir*ir;
			}

			//evaluate damping factors
			switch (damp_type) {
			case DAMPING_OFF:
				if (pair_ptr->es_excluded)
					damp1 = damp2 = wdamp1 = wdamp2 = 0.0;
				else
					damp1 = damp2 = wdamp1 = wdamp2 = 1.0;
				break;
			case DAMPING_LINEAR:
				s = l * pow(atom_array[i]->polarizability*atom_array[j]->polarizability, 1.0 / 6.0);
				v = r / s;
				if (r < s) {
					damp1 = (4.0 - 3.0*v)*v*v*v;
					damp2 = v * v*v*v;
				}
				else {
					damp1 = damp2 = 1.0;
				}
				break;
			case DAMPING_EXPONENTIAL:
				explr = exp(-l * r);
				damp1 = 1.0 - explr * (0.5*l2*r2 + l * r + 1.0);
				damp2 = damp1 - explr * (l3*r2*r / 6.0);
				if (polar_wolf_full) { //subtract off damped interaction at r_cutoff
					wdamp1 = 1.0 - explrcut * (0.5*l2*rcut2 + l * rcut + 1.0);
					wdamp2 = wdamp1 - explrcut * (l3*rcut3 / 6.0);
				}
				break;
			default:
				Output::err("error: something unexpected happened in thole_matrix.c");
			}

			// build the tensor
			for (int p = 0; p < 3; p++) {
				for (int q = 0; q < 3; q++) {

					A_matrix[ii + p][jj + q] = -3.0*pair_ptr->dimg[p] * pair_ptr->dimg[q] * damp2*ir5;
					if (polar_wolf_full)
						A_matrix[ii + p][jj + q] -= -3.0*pair_ptr->dimg[p] * pair_ptr->dimg[q] * wdamp2*ir*ir / rcut3;

					// additional diagonal term
					if (p == q) {
						A_matrix[ii + p][jj + q] += damp1 * ir3;
						if (polar_wolf_full)
							A_matrix[ii + p][jj + q] -= wdamp1 / (rcut3);
					}
				}
			}

			// set the lower half of the tensor component 
			for (int p = 0; p < 3; p++)
				for (int q = 0; q < 3; q++)
					A_matrix[jj + p][ii + q] = A_matrix[ii + p][jj + q];

		} // end j 
	} // end i 

	return;
}



void System::zero_out_amatrix (int NAtoms) {

	// zero out the matrix 
	for (int i = 0; i < 3 * NAtoms; i++)
		for (int j = 0; j < 3 * NAtoms; j++)
			A_matrix[i][j] = 0;
	return;
}

//do full polarization calculation using ewald
//see nymand and linse jcp 112 6152 (2000)
void System::ewald_full() {

	//int max_iter=system->polar_max_iter;  (unused variable)
	int keep_iterating, iteration_counter;

	//calculate static e-field
	zero_out(molecules);
	recip_term();
	real_term();

	//calculate induced e-field
	init_dipoles_ewald();

	keep_iterating = 1;
	iteration_counter = 0;
	while (keep_iterating) {

		if ((iteration_counter >= (int)(MAX_ITERATION_COUNT)) && polar_precision) {
			iterator_failed = 1;
			return;
		}

		//set induced field to zero
		clear_ef_induced();

		//calculate induced field
		induced_real_term();
		induced_recip_term();
		induced_corr_term();

		if (polar_rrms || polar_precision > 0)
			calc_dipole_rrms();

		//recalculate dipoles using new induced field
		new_dipoles(iteration_counter);

		keep_iterating = are_we_done_yet(iteration_counter);

		if (polar_palmo && !keep_iterating) //if last iteration
			ewald_palmo_contraction();

		iteration_counter++;
	}

	return;
}

//we deviate from drexel's treatment, and instead do a trig identity to get from a pairwise sum to two atomwise rums
//or ignore drexel, and derive this term from eq (29) in nymand and linse
void System::recip_term() {

	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	int        l[3] = { 0 },
		kmax = ewald_kmax;
	double     ea = polar_ewald_alpha, //actually sqrt(ea)
		k[3] = { 0 },
		k2 = 0,
		kweight[3] = { 0 },
		float1 = 0,
		float2 = 0;


	//k-space sum (symmetry for k -> -k, so we sum over hemisphere, avoiding double-counting on the face)
	for (l[0] = 0; l[0] <= kmax; l[0]++) {
		for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
			for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

				// if |l|^2 > kmax^2, then it doesn't contribute (outside the k-space cutoff)
				if (UsefulMath::iidotprod(l, l) > kmax*kmax) continue;

				for (int p = 0; p < 3; p++) {
					k[p] = 0;
					for (int q = 0; q < 3; q++)
						k[p] += 2.0 * pi * pbc.reciprocal_basis[p][q] * l[q];
				}
				k2 = UsefulMath::dddotprod(k, k);

				kweight[0] = k[0] / k2 * exp(-k2 / (4.0*ea*ea));
				kweight[1] = k[1] / k2 * exp(-k2 / (4.0*ea*ea));
				kweight[2] = k[2] / k2 * exp(-k2 / (4.0*ea*ea));

				float1 = float2 = 0;
				for (mptr = molecules; mptr; mptr = mptr->next)
					for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
						float1 += aptr->charge * cos(UsefulMath::dddotprod(k, aptr->pos));
						float2 += aptr->charge * sin(UsefulMath::dddotprod(k, aptr->pos));
					}

				for (mptr = molecules; mptr; mptr = mptr->next)
					for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
						for (int p = 0; p < 3; p++) {
							aptr->ef_static[p] += kweight[p] * sin(UsefulMath::dddotprod(k, aptr->pos)) * float1;
							aptr->ef_static[p] -= kweight[p] * cos(UsefulMath::dddotprod(k, aptr->pos)) * float2;
						}
					}

			} //l2
		} //l1
	} //l0

	for (mptr = molecules; mptr; mptr = mptr->next) {
		for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
			for (int p = 0; p < 3; p++) {
				//factor of 2 more, since we only summed over hemisphere
				aptr->ef_static[p] *= 8.0*pi / pbc.volume;
			}
		}
	}

	return;
}



void System::real_term() {

	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	Pair     * pptr = nullptr;

	double     r = 0,
		r2 = 0,
		factor = 0,
		a = polar_ewald_alpha; //some ambiguity between ea and ea^2 across the literature;


	for (mptr = molecules; mptr; mptr = mptr->next) {
		for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
			for (pptr = aptr->pairs; pptr; pptr = pptr->next) { //for each pair
				if (pptr->frozen) continue; //if the pair is frozen (i.e. MOF-MOF interaction) it doesn't contribute to polar
				r = pptr->rimg;
				if ((r > pbc.cutoff) || (r == 0.0)) continue; //if outside cutoff sphere (not sure why r==0 ever) -> skip
				r2 = r * r;
				if (pptr->es_excluded) {
					//need to subtract self-term (interaction between a site and a neighbor's screening charge (on the same molecule)
					factor = (2.0 * a * OneOverSqrtPi * exp(-a * a*r2) * r - erf(a*r)) / (r*r2);
					for (int p = 0; p < 3; p++) {
						aptr->ef_static[p] += factor * pptr->atom->charge * pptr->dimg[p];
						pptr->atom->ef_static[p] -= factor * aptr->charge * pptr->dimg[p];
					}
				} //excluded
				else { //not excluded

					factor = (2.0 * a * OneOverSqrtPi * exp(-a * a*r2) * r + erfc(a*r)) / (r2*r);
					for (int p = 0; p < 3; p++) { // for each dim, add e-field contribution for the pair
						aptr->ef_static[p] += factor * pptr->atom->charge * pptr->dimg[p];
						pptr->atom->ef_static[p] -= factor * aptr->charge * pptr->dimg[p];
					}
				} //excluded else
			} //ptr
		} //aptr
	} //mptr

	return;
}


//set zeroth iteration dipoles
void System::init_dipoles_ewald() {
	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;

	for (mptr = molecules; mptr; mptr = mptr->next)
		for (aptr = mptr->atoms; aptr; aptr = aptr->next)
			for (int p = 0; p < 3; p++) {
				aptr->old_mu[p] = 0;
				aptr->new_mu[p] = aptr->mu[p] = aptr->polarizability*aptr->ef_static[p];
			}

	return;
}


//reset the ef_induced values
void System::clear_ef_induced() {
	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;


	for (mptr = molecules; mptr; mptr = mptr->next)
		for (aptr = mptr->atoms; aptr; aptr = aptr->next)
			for (int p = 0; p < 3; p++)
				aptr->ef_induced[p] = 0;

	return;
}



void System::induced_recip_term() {
	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr,
		** aarray = nullptr;
	int        NAtoms = 0,
		l[3] = { 0 },
		kmax = ewald_kmax;
	double     Psin = 0,
		Pcos = 0,
		kweight = 0,
		a = polar_ewald_alpha,
		dotprod1 = 0,
		dotprod2 = 0,
		k[3] = { 0 },
		k2 = 0;

	//make atom array
	for (mptr = molecules; mptr; mptr = mptr->next) {
		for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
			SafeOps::realloc(aarray, sizeof(Atom *) * (NAtoms + 1), __LINE__, __FILE__);
			aarray[NAtoms] = aptr;
			NAtoms++;
		}
	}

	//k-space sum (symmetry for k -> -k, so we sum over hemisphere, avoiding double-counting on the face)
	for (l[0] = 0; l[0] <= kmax; l[0]++)
		for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++)
			for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

				// if |l|^2 > kmax^2, then it doesn't contribute (outside the k-space cutoff)
				if (UsefulMath::iidotprod(l, l) > kmax*kmax) continue;

				for (int p = 0; p < 3; p++) {
					k[p] = 0;
					for (int q = 0; q < 3; q++)
						k[p] += 2.0 * pi * pbc.reciprocal_basis[p][q] * l[q];
				}
				k2 = UsefulMath::dddotprod(k, k);

				for (int p = 0; p < 3; p++)
					kweight = 8.0 * pi / pbc.volume * exp(-k2 / (4.0*a*a)) / k2 * k[p];

				//calculate Pcos, Psin for this k-point
				Pcos = Psin = 0;
				for (int j = 0; j < NAtoms; j++) {
					dotprod1 = UsefulMath::dddotprod(k, aarray[j]->mu);
					dotprod2 = UsefulMath::dddotprod(k, aarray[j]->pos);
					Pcos += dotprod1 * cos(dotprod2);
					Psin += dotprod1 * sin(dotprod2);
				}

				//calculate ef_induced over atom array
				for (int i = 0; i < NAtoms; i++) {
					//for each cartesian dimension
					for (int p = 0; p < 3; p++) {

						dotprod1 = UsefulMath::dddotprod(k, aarray[i]->pos);
						aarray[i]->ef_induced[p] += kweight * (-sin(dotprod1)*Psin - cos(dotprod1)*Pcos);

					} //dim
				} //ef_incuded over atom array
			} //kspace	

	free(aarray);

	return;
}



void System::induced_real_term() {

	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	Pair     * pptr = nullptr;
	double     erfcar = 0,
		expa2r2 = 0,
		T = 0, //dipole-interaction tensor component
		s1 = 0,
		s2 = 0, //common term (s_2 via eq 10. JCP 133 243101)
		a = polar_ewald_alpha, //ewald damping
		l = polar_damp, //polar damping
		r = 0, ir = 0, ir3 = 0, ir5 = 0;



	for (mptr = molecules; mptr; mptr = mptr->next) {
		for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
			for (pptr = aptr->pairs; pptr; pptr = pptr->next) {
				if (aptr->polarizability == 0 || pptr->atom->polarizability == 0)
					continue; //don't waste CPU time
				if (pptr->rimg > pbc.cutoff)
					continue; //if outside cutoff sphere skip
				//some things we'll need
				r = pptr->rimg;
				ir = 1.0 / r; ir3 = ir * ir*ir; ir5 = ir * ir*ir3;
				erfcar = erfc(a*r);
				expa2r2 = exp(-a * a*r*r);

				//E_static_realspace_i = sum(i!=j) d_xi d_xj erfc(a*r)/r u_j 
				s2 = erfcar + 2.0*a*r*OneOverSqrtPi * expa2r2 + 4.0*a*a*a*r*r*r / 3.0*OneOverSqrtPi*expa2r2 - damp_factor(l*r, 3);

				for (int p = 0; p < 3; p++) {
					for (int q = p; q < 3; q++) {  //it's symmetric!

						if (p == q)
							s1 = erfcar + 2.0*a*r*OneOverSqrtPi * expa2r2 - damp_factor(l*r, 2);
						else
							s1 = 0;

						//real-space dipole interaction tensor
						T = 3.0 * pptr->dimg[p] * pptr->dimg[q] * s2 * ir5 - s1 * ir3;

						aptr->ef_induced[p] += T * pptr->atom->mu[q];
						pptr->atom->ef_induced[p] += T * aptr->mu[q];

						if (p != q) {
							aptr->ef_induced[q] += T * pptr->atom->mu[p];
							pptr->atom->ef_induced[q] += T * aptr->mu[p];
						}

					} //loop over q dim
				} //loop over p dim
			} //pptr loop
		} //aptr loop
	} //mptr loop

	return;
}

//damping term (see e.g. Souaille et al., Comp. Phys. Comm. 180 276-301) below eq (9).
//signs are intentionally reversed (different convention)
double System::damp_factor(double t, int i) {

	double temp = 1.0 + t + 0.5 * t*t;

	if (i == 3)
		temp += t * t*t / 6.0;

	return temp * exp(-t);
}



void System::induced_corr_term() {
	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;

	double     a = polar_ewald_alpha,
		totalmu[3] = { 0 };

	for (int p = 0; p < 3; p++)
		totalmu[p] = 0;

	for (mptr = molecules; mptr; mptr = mptr->next)
		for (aptr = mptr->atoms; aptr; aptr = aptr->next)
			for (int p = 0; p < 3; p++)
				totalmu[p] += aptr->mu[p];

	//other term
	for (mptr = molecules; mptr; mptr = mptr->next)
		for (aptr = mptr->atoms; aptr; aptr = aptr->next)
			for (int p = 0; p < 3; p++)

				aptr->ef_induced[p] += -4.0 * pi / (3.0*pbc.volume) * totalmu[p] + 4.0*a*a*a / (3.0*SqrtPi) * aptr->mu[p];

	return;
}



void System::calc_dipole_rrms() {

	Atom ** aa = atom_array;
	double  carry = 0;

	// get the dipole RRMS 
	for (int i = 0; i < natoms; i++) {
		// mean square difference 
		aa[i]->dipole_rrms = 0;
		for (int p = 0; p < 3; p++) {
			carry = aa[i]->new_mu[p] - aa[i]->old_mu[p];
			aa[i]->dipole_rrms += carry * carry;
		}

		// normalize
		aa[i]->dipole_rrms /= UsefulMath::dddotprod(aa[i]->new_mu, aa[i]->new_mu);
		aa[i]->dipole_rrms = sqrt(aa[i]->dipole_rrms);
		if (!std::isfinite(aa[i]->dipole_rrms))
			aa[i]->dipole_rrms = 0;
	}

#ifdef DEBUG
	double totalrms = 0;
	for (int i = 0; i < system->natoms; i++) {
		totalrms += aa[i]->dipole_rrms;
	}
	fprintf(stderr, "TOTAL DIPOLE RRMS %lf\n", totalrms);
#endif

	return;
}



void System::new_dipoles(int count) {

	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;

	for (mptr = molecules; mptr; mptr = mptr->next) {
		for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
			for (int p = 0; p < 3; p++) {

				//set dipoles
				aptr->old_mu[p] = aptr->mu[p];
				if (polar_sor) {
					aptr->new_mu[p] = aptr->polarizability*(aptr->ef_static[p] + aptr->ef_induced[p]);
					aptr->new_mu[p] = aptr->mu[p] = polar_gamma * aptr->new_mu[p] + (1.0 - polar_gamma)*aptr->old_mu[p];
				}
				else if (polar_esor) {
					aptr->new_mu[p] = aptr->polarizability*(aptr->ef_static[p] + aptr->ef_induced[p]);
					aptr->new_mu[p] = aptr->mu[p] = (1.0 - exp(-polar_gamma * (count + 1)))*aptr->new_mu[p] +
						exp(-polar_gamma * (count + 1))*aptr->old_mu[p];
				}
				else {
					//if no sor, still need new_mu for polar_palmo
					aptr->mu[p] = aptr->new_mu[p] = aptr->polarizability*(aptr->ef_static[p] + aptr->ef_induced[p]);
				}

			}
		}
	}

	return;
}



int System::are_we_done_yet(int iteration_counter) {

	Atom ** aa = atom_array;
	int     NAtoms = natoms;
	double  allowed_sqerr = 0,
		error = 0;

	if (polar_precision == 0.0) {	// ... by fixed iteration ...
		if (iteration_counter != polar_max_iter)
			return 1;
	}

	else { // ... or by dipole precision 
		allowed_sqerr = polar_precision * polar_precision * DEBYE2SKA * DEBYE2SKA;
		for (int i = 0; i < NAtoms; i++) { //check the change in each dipole component
			for (int p = 0; p < 3; p++) {
				error = aa[i]->new_mu[p] - aa[i]->old_mu[p];
				if (error*error > allowed_sqerr)
					return 1; //we broke tolerance
			}
		}
	}

	return 0;
}



void System::ewald_palmo_contraction() {

	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;

	//set induced field to zero
	clear_ef_induced();

	//calculate induced field
	induced_real_term();
	induced_recip_term();
	induced_corr_term();

	for (mptr = molecules; mptr; mptr = mptr->next) {
		for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
			if (aptr->polarizability == 0) continue;
			for (int p = 0; p < 3; p++) {
				aptr->ef_induced_change[p] = //current induced - last induced (backed out from dipole values)
					aptr->ef_induced[p] - (aptr->new_mu[p] / aptr->polarizability - aptr->ef_static[p]);
			}
		}
	}

	return;
}


// calculate the field with periodic boundaries
void System::thole_field() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;

	// zero the field vectors 
	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			for (int p = 0; p < 3; p++) {
				atom_ptr->ef_static[p] = 0;
				atom_ptr->ef_static_self[p] = 0;
			}

		}
	}

	// calculate the electrostatic field 
	if (polar_ewald)
		ewald_estatic();
	else if (polar_wolf || polar_wolf_full)
		thole_field_wolf();
	else
		thole_field_nopbc();

}


// calculate the field without ewald summation/wolf 
void System::thole_field_nopbc() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	Pair     * pair_ptr = nullptr;
	double     r = 0;

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (pair_ptr->frozen)
					continue;
				if (molecule_ptr == pair_ptr->molecule)
					continue; //don't let molecules polarize themselves

				r = pair_ptr->rimg;

				//inclusive near the cutoff
				if ((r - SMALL_dR < pbc.cutoff) && (r != 0.)) {

					for (int p = 0; p < 3; p++) {
						atom_ptr->ef_static[p] += pair_ptr->atom->charge*pair_ptr->dimg[p] / (r*r*r);
						pair_ptr->atom->ef_static[p] -= atom_ptr->charge*pair_ptr->dimg[p] / (r*r*r);
					}

				} // cutoff 

			} // pair
		} // atom
	} // molecule

	return;
}


// calc field using wolf sum (JCP 124 234104 (2006) equation 19
void System::thole_field_wolf() {

	Molecule  * molecule_ptr = nullptr;
	Atom      * atom_ptr = nullptr;
	Pair      * pair_ptr = nullptr;

	double      r = 0,
		rr = 0, // 1/r (reciprocal of r)
		R = pbc.cutoff,
		rR = 1. / R;
	//used for polar_wolf_alpha (aka polar_wolf_damp)
	double      a = polar_wolf_alpha,
		erR = erfc(a*R), //complementary error functions
		cutoffterm = (erR*rR*rR + 2.0*a*OneOverSqrtPi*exp(-a * a*R*R)*rR),
		bigmess = 0;

	//init lookup table if needed
	if (polar_wolf_alpha_lookup && !(polar_wolf_alpha_table))
		polar_wolf_alpha_table = polar_wolf_alpha_lookup_init();

	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if (molecule_ptr == pair_ptr->molecule)
					continue; //don't let molecules polarize themselves
				if (pair_ptr->frozen)
					continue; //don't let the MOF polarize itself

				r = pair_ptr->rimg;

				if ((r - SMALL_dR < pbc.cutoff) && (r != 0.)) {
					rr = 1. / r;

					//we will need this shit if wolf alpha != 0 
					if ((a != 0) & polar_wolf_alpha_lookup)
						bigmess = polar_wolf_alpha_getval(r);
					else if (a != 0) //no lookup  
						bigmess = (erfc(a*r)*rr*rr + 2.0*a*OneOverSqrtPi*exp(-a * a*r*r)*rr);

					for (int p = 0; p < 3; p++) {
						//see JCP 124 (234104)
						if (a == 0) {
							atom_ptr->ef_static[p] += (pair_ptr->atom->charge)*(rr*rr - rR * rR)*pair_ptr->dimg[p] * rr;
							pair_ptr->atom->ef_static[p] -= (atom_ptr->charge)*(rr*rr - rR * rR)*pair_ptr->dimg[p] * rr;
						}
						else {
							atom_ptr->ef_static[p] += pair_ptr->atom->charge*(bigmess - cutoffterm)*pair_ptr->dimg[p] * rr;
							pair_ptr->atom->ef_static[p] -= atom_ptr->charge*(bigmess - cutoffterm)*pair_ptr->dimg[p] * rr;
						}

					}

				} // no lookup table
			} // pair
		} // atom
	} // molecule

	return;
}

//only calculate the static e-field via ewald
//see http://www.pages.drexel.edu/~cfa22/msim/node50.html
void System::ewald_estatic() {

	//calculate static e-field
	// no need to zero out dipoles; this is done in polar.c
	recip_term();
	real_term();

	return;
}



// the point of this is to store all the polar_wolf_alpha calculations in a big array, and then just look shit up
// that way we don't need to calculate erfc's and exp's over and over and over
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double * System::polar_wolf_alpha_lookup_init () {
	double * rval = nullptr;
	int      i = 0;
	double   a = polar_wolf_alpha,
		r = 0,
		rr = 0;

	polar_wolf_alpha_table_max = (int)(ceil(polar_wolf_alpha_lookup_cutoff)) * 1000;
	SafeOps::malloc(rval, sizeof(double) * polar_wolf_alpha_table_max, __LINE__, __FILE__);

	for (i = 1; i < polar_wolf_alpha_table_max; i++) {
		r = (double)i / 1000.;
		rr = 1.0 / r;
		rval[i] = erfc(a*r)*rr*rr + 2.0*a*OneOverSqrtPi*exp(-a * a*r*r)*rr;
	}
	//store the zeroth component without blowing up
	rval[0] = rval[1];


	return rval;
}

double System::polar_wolf_alpha_getval(double r) {

	int i = (int)(r * 1000);

	if (i >= polar_wolf_alpha_table_max)
		return 0.0; //answer will be zero if cutoff is large enough

	return polar_wolf_alpha_table[i];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// iterative solver of the dipole field tensor returns the number of iterations required 
int System::thole_iterative() {

	int     NAtoms = natoms,
		iteration_counter = 0,
		keep_iterating = 1;
	Atom ** aa = atom_array;
	int   * ranked_array = nullptr;




	// array for ranking
	SafeOps::calloc(ranked_array, NAtoms, sizeof(int), __LINE__, __FILE__);
	for (int i = 0; i < NAtoms; i++)
		ranked_array[i] = i;

	//set all dipoles to alpha*E_static * polar_gamma
	init_dipoles();

	// if ZODID is enabled, then stop here and just return the alpha*E dipoles
	if (polar_zodid) {
		free(ranked_array);
		return(0);
	}

	// iterative solver of the dipole field equations 
	keep_iterating = 1;
	iteration_counter = 0;
	while (keep_iterating) {
		iteration_counter++;

		// divergence detection
		// if we fail to converge, then return dipoles as alpha*E
		if (iteration_counter >= MAX_ITERATION_COUNT && polar_precision) {
			for (int i = 0; i < NAtoms; i++)
				for (int p = 0; p < 3; p++) {
					aa[i]->mu[p] = aa[i]->polarizability * (aa[i]->ef_static[p] + aa[i]->ef_static_self[p]);
					aa[i]->ef_induced_change[p] = 0.0; //so we don't break palmo
				}
			//set convergence failure flag
			iterator_failed = 1;

			free(ranked_array);
			return(iteration_counter);
		}

		//zero out induced e-field
		for (int i = 0; i < natoms; i++)
			for (int p = 0; p < 3; p++)
				aa[i]->ef_induced[p] = 0;

		//save the current dipole information if we want to calculate precision (or if needed for relaxation)
		if (polar_rrms || polar_precision > 0 || polar_sor || polar_esor) {
			for (int i = 0; i < NAtoms; i++)
				for (int p = 0; p < 3; p++)
					aa[i]->old_mu[p] = aa[i]->mu[p];
		}

		// contract the dipoles with the field tensor (gauss-seidel/gs-ranked optional)
		contract_dipoles(ranked_array);

		if (polar_rrms || polar_precision > 0)
			calc_dipole_rrms();

		// determine if we are done... 
		keep_iterating = are_we_done_yet(iteration_counter);

		// if we would be finished, we want to contract once more to get the next induced field for palmo
		if (polar_palmo && !keep_iterating)
			palmo_contraction(ranked_array);

		//new gs_ranking if needed
		if (polar_gs_ranked && keep_iterating)
			update_ranking(ranked_array);

		// save the dipoles for the next pass
		for (int i = 0; i < NAtoms; i++) {
			for (int p = 0; p < 3; p++) {
				// allow for different successive over-relaxation schemes
				if (polar_sor)
					aa[i]->mu[p] = polar_gamma * aa[i]->new_mu[p] + (1.0 - polar_gamma)*aa[i]->old_mu[p];
				else if (polar_esor)
					aa[i]->mu[p] = (1.0 - exp(-polar_gamma * iteration_counter))*aa[i]->new_mu[p] + exp(-polar_gamma * iteration_counter)*aa[i]->old_mu[p];
				else
					aa[i]->mu[p] = aa[i]->new_mu[p];
			}
		}

	} //end iterate
	free(ranked_array);

	// return the iteration count
	return(iteration_counter);
}


// set them to alpha*E_static
void System::init_dipoles() {

	Atom ** aa = atom_array;

	for (int i = 0; i < natoms; i++) {
		for (int p = 0; p < 3; p++) {
			aa[i]->mu[p] = aa[i]->polarizability*(aa[i]->ef_static[p] + aa[i]->ef_static_self[p]);
			// should improve convergence since mu's typically grow as induced fields are added in
			if (!polar_sor && !polar_esor)
				aa[i]->mu[p] *= polar_gamma;
		}
	}
	return;
}



void System::contract_dipoles(int * ranked_array) {
	int     ii = 0,
		jj = 0,
		index = 0;
	Atom ** aa = atom_array;

	for (int i = 0; i < natoms; i++) {
		index = ranked_array[i]; //do them in the order of the ranked index
		ii = index * 3;
		if (aa[index]->polarizability == 0) { //if not polar
			//aa[index]->ef_induced[p] is already 0
			aa[index]->new_mu[0] = aa[index]->new_mu[1] = aa[index]->new_mu[2] = 0; //might be redundant?
			aa[index]->mu[0] = aa[index]->mu[1] = aa[index]->mu[2] = 0; //might be redundant?
			continue;
		}
		for (int j = 0; j < natoms; j++) {
			jj = j * 3;
			if (index != j)
				for (int p = 0; p < 3; p++)
					aa[index]->ef_induced[p] -= UsefulMath::dddotprod((A_matrix[ii + p] + jj), aa[j]->mu);
		} // end j 

		// dipole is the sum of the static and induced parts
		for (int p = 0; p < 3; p++) {
			aa[index]->new_mu[p] = aa[index]->polarizability*(aa[index]->ef_static[p] + aa[index]->ef_static_self[p] + aa[index]->ef_induced[p]);

			// Gauss-Seidel
			if (polar_gs || polar_gs_ranked)
				aa[index]->mu[p] = aa[index]->new_mu[p];
		}

	} // end matrix multiply

	return;
}



void System::palmo_contraction(int * ranked_array) {

	int     ii = 0,
		jj = 0,
		index = 0,
		NAtoms = natoms;
	Atom ** aa = atom_array;

	// calculate change in induced field due to this iteration
	for (int i = 0; i < NAtoms; i++) {
		index = ranked_array[i];
		ii = index * 3;

		for (int p = 0; p < 3; p++)
			aa[index]->ef_induced_change[p] = -aa[index]->ef_induced[p];

		for (int j = 0; j < NAtoms; j++) {
			jj = j * 3;
			if (index != j)
				for (int p = 0; p < 3; p++)
					aa[index]->ef_induced_change[p] -= UsefulMath::dddotprod(A_matrix[ii + p] + jj, aa[j]->mu);
		}
	}

	return;
}



void System::update_ranking(int * ranked_array) {

	int     sorted = 0,
		tmp = 0,
		NAtoms = natoms;
	Atom ** aa = atom_array;

	// rank the dipoles by bubble sort
	if (polar_gs_ranked) {
		for (int i = 0; i < NAtoms; i++) {
			sorted = 1;
			for (int j = 0; j < (NAtoms - 1); j++) {
				if (aa[ranked_array[j]]->rank_metric < aa[ranked_array[j + 1]]->rank_metric) {
					sorted = 0;
					tmp = ranked_array[j];
					ranked_array[j] = ranked_array[j + 1];
					ranked_array[j + 1] = tmp;
				}
			}
			if (sorted)
				break;
		}
	}

	return;
}


// invert the A matrix
void System::thole_bmatrix() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr = nullptr;
	int        NAtoms = 0;

	// count the number of atoms
	for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			++NAtoms;

	UsefulMath::invert_matrix(3 * NAtoms, A_matrix, B_matrix);
}


// get the dipoles by vector matrix multiply
void System::thole_bmatrix_dipoles() {

	int      ii = 0,
		NAtoms = natoms;
	double * mu_array = nullptr,
		*field_array = nullptr;

	// allocate working arrays
	SafeOps::calloc(mu_array, 3 * NAtoms, sizeof(double), __LINE__, __FILE__);
	SafeOps::calloc(field_array, 3 * NAtoms, sizeof(double), __LINE__, __FILE__);

	// copy the field in
	for (int i = 0; i < NAtoms; i++) {
		ii = i * 3;
		for (int p = 0; p < 3; p++)
			field_array[ii + p] = atom_array[i]->ef_static[p] + atom_array[i]->ef_static_self[p];
	}

	// multiply the supervector with the B matrix
	for (int i = 0; i < 3 * NAtoms; i++)
		for (int j = 0; j < 3 * NAtoms; j++)
			mu_array[i] += B_matrix[i][j] * field_array[j];

	// copy the dipoles out
	for (int i = 0; i < NAtoms; i++) {
		ii = i * 3;
		for (int p = 0; p < 3; p++)
			atom_array[i]->mu[p] = mu_array[ii + p];
	}

	// free the working arrays
	free(mu_array);
	free(field_array);

}


// calculate the molecular polarizability tensor from the B matrix 
void System::thole_polarizability_tensor() {
	char linebuf[maxLine];
	int ii, jj, NAtoms;

	double isotropic;

	NAtoms = checkpoint->thole_N_atom;

	// clear the polarizability tensor
	for (int p = 0; p < 3; p++)
		for (int q = 0; q < 3; q++)
			C_matrix[p][q] = 0;

	// sum the block terms for the 3x3 molecular tensor
	for (int p = 0; p < 3; p++) {
		for (int q = 0; q < 3; q++) {
			for (int i = 0; i < NAtoms; i++) {
				ii = i * 3;
				for (int j = 0; j < NAtoms; j++) {
					jj = j * 3;
					C_matrix[p][q] += B_matrix[ii + p][jj + q];
				}
			}
		}
	}

	// get the isotropic term
	isotropic = 0;
	for (int p = 0; p < 3; p++)
		isotropic += C_matrix[p][p];
	isotropic /= 3.0;

	Output::out("POLARIZATION: polarizability tensor (A^3):\n");
	Output::out("##########################\n");
	for (int p = 0; p < 3; p++) {
		for (int q = 0; q < 3; q++) {
			sprintf(linebuf, "%.4f ", C_matrix[p][q]);
			Output::out(linebuf);
		}
		Output::out("\n");
	}

	Output::out("##########################\n");
	sprintf(linebuf, "isotropic = %.4f\n", isotropic);
	Output::out(linebuf);
	sprintf(linebuf, "XX/ZZ = %.4f\n", C_matrix[0][0] / C_matrix[2][2]);
	Output::out(linebuf);

}