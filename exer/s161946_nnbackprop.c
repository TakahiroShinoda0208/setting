/* M.Nielsen April 2002 mniel@cbs.dtu.dk */

/* 
Copyright (C) 2008-2015 Danish Technical University

This suite of programs and library routine is free software. You can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

In other words, you are free to modify, copy, or redistribute this
source code and its documentation in any way you like, but you must
distribute all derivative versions as free software under the same
terms that I've provided my code to you (i.e. the GNU General Public
License). This precludes any use of the code in proprietary or
commercial software unless your source code is made freely available.

If you wish to use the code under a different Open Source license
that's not compatible with the GPL (like the Artistic License, BSD
license, or the Netscape Public License), please contact me
(Morten Nielsen, mniel@cbs.dtu.dk) for permission.

Incorporation into commercial software under non-GPL terms is possible
by obtaining a specially licensed version from The Danish Technical University.
Contact Morten Nielsen (mniel@cbs.dtu.dk) to arrange licensing terms.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

static double squash( a ) double a; { return (1.0 / (1.0 + exp( -a ))); }

int	p_verbose;
int	p_nhid;
FILENAME	p_syn;
int     p_seed;
int	p_ncycles;
float	p_eta;
FILENAME	p_trainpred;
FILENAME        p_testpred;
FILENAME	p_testfile;
int     p_ntest;
float	p_bplim;
float	p_weight;
int	p_dtype;

PARAM   param[] = {
	"-v", VSWITCH	p_verbose,	"Verbose mode", "0",
	"-nh", VINT	p_nhid, "Number of hidden neurons", "2",
	"-syn", VFNAME	p_syn,	"Name of synaps file", "syn.dat",
        "-s", VINT      p_seed, "Seed [-1] Default, [0] Time [>0] Specific seed", "-1",
	"-nc", VINT     p_ncycles, "Number of iterations", "500",
	"-eta", VFLOAT	p_eta, "Eta for weight update", "0.05",
	"-ol", VFNAME	p_trainpred, "File for training data prediction output", "trainpred.out",
	"-ot", VFNAME   p_testpred, "File for test data prediction output", "testpred.out",
	"-tf", VFNAME	p_testfile, "File with test input", "",
	"-nt", VINT     p_ntest, "Test interval", "10",
	"-bl", VFLOAT	p_bplim, "Limit for backpropagation", "0.00001",
	"-w",  VFLOAT	p_weight, "Initial value for weight", "0.1",
	"-dtype", VINT 	p_dtype, "Dump type [0] Error [1] Pearson", "0",
	0
};

float	*scan_input( char *line, int *n )

{
	float	*v;

	/* Split input line on blanks */
	v = fvector_split( line, n );

	return( v );
}

void    setnet(float **x, SYNAPS * syn)

{
        int     l;

	/* Initialize bias input neuron */
        for (l = 0; l < syn->nlay; l++)
                x[l][syn->nnlay[l]] = 1.0;
}

void	forward( float **x, SYNAPS *s, float *inp )

{

	int	i,l,k;
	float	a;

	setnet( x, s );

	for ( i=0; i<s->nnlay[0]; i++ ) 
		x[0][i] = inp[i];

	for ( l = 1; l < s->nlay; l++ ) 
	for ( k = 0; k < s->nnlay[l]; k++ ) {
		a = 0.0;
		for (i = 0; i < s->nnlay[l-1]+1; i++ ) {
			a += x[l-1][i] * s->wgt[l-1][k][i];
		}
		x[l][k] = squash( a );
	}

}

void	backprop( float **x, SYNAPS *s, float *inp, float t, float *dw, float **dv )


{
	int	j,k;
	float	***wgt;
	float	Hj, delta, d_o;
	float	gprime_o, gprime_h;
	float	tmp1;
	float	w0j, Ik;

	wgt = s->wgt;

	/* back propagation in output layer */
	/* The backprop only works for 3 layer networks with one output neuron! */
	/* x[0][i] is the input */
	/* x[1][j] is the output from the hidden layer - H_j */
	/* x[2][0] is the output from the network */

	/* From lecture notes

	delta = ( O - t )* g'(o)
	dE/dw_j = delta * H_j 
	dE_dv_jk = g'(h_j) * I_k * delta * w_j

	*/

	/* 
	d_o is the difference between output and target
	gprime_o is g'(o)

	*/

	/* d_o = XXXXX; */
	/* gprime_o = XXX; */

	/* delta = XXXX; */

	d_o = (x[2][0] - t); //inputがおそらくmeasure
	gprime_o = (1-x[2][0])*x[2][0];
	delta = d_o*gprime_o;

	/* Calculate change to weights between hidden and output layer */

	for ( j=0; j<=s->nnlay[1]; j++ ) {

              //dw[j] = XXXX;
              dw[j] = delta * x[1][j];
                
	}

	/* Update weight between hidden and output layer */

	for ( j=0; j<=s->nnlay[1]; j++ )
              //wgt[1][0][j] -= XXXXX; ?? substract??
              wgt[1][0][j] -= p_eta*dw[j];

	/* back propagation in input layer */

	/* 
	Hj is the output from the hidden layer
	gprime_h is g'(h_j)
	*/

	for ( j = 0; j<s->nnlay[1]; j++ ) {

		/* Hj = XXXX; */
		/* w0j = XXXX; */
		/* gprime_h = XXXX; */

		Hj = x[1][j];
		w0j = wgt[1][0][j];
		gprime_h = (1-Hj)*Hj;

		tmp1 = w0j * gprime_h * delta;

		for ( k = 0; k<=s->nnlay[0]; k++ ) {

                      //Ik = XXXX;
                      Ik = x[0][k];

			dv[j][k] = tmp1 * Ik;

		}
	}

	/* Update weight between input and hidden layer */

	for ( k=0; k<=s->nnlay[0]; k++ )	
	for ( j=0; j<s->nnlay[1]; j++ )
              wgt[0][j][k] -= p_eta*dv[j][k];

}

float	***init_wgt( SYNAPS *syn )

{
	float	***wgt;
	int	l, k, i;

	wgt = (float ***) malloc((unsigned) (syn->nlay - 1) * sizeof(float **));

	if ( ! wgt ) {
		printf( "Allocation failure 1 in synaps_wmtx_alloc\n");
                exit(1);
        } 

	for (l = 0; l < syn->nlay - 1; l++) {

		wgt[l] = fmatrix(0, syn->nnlay[l + 1] - 1, 0, syn->nnlay[l]);
               
                if (!wgt[l]) {
                        printf( "Allocation failure 2 in synaps_wmtx_alloc\n");
                        exit(1);
                }

		for (k = 0; k < syn->nnlay[l + 1]; k++)
                for (i = 0; i < syn->nnlay[l] + 1; i++) 
			wgt[l][k][i] = p_weight * ( 0.5 - drand48());
	}

	return( wgt );
}

float	**init_input( char *filename, int *ncolumn, int *ndata )

{
	LONGLINELIST	*linelist, *ln;
	int	n;
	int	nc, nc0;
	float	**inp;

	linelist = longlinelist_read( filename );

	if ( ! linelist ) {
		printf( "Error. Cannot read input data from file %s. Exit\n", filename );
                exit( 1 );
        }

	for ( n=0,ln = linelist; ln; ln=ln->next, n++ );

	inp = (float **) malloc((unsigned) ( n ) * sizeof(float * ));

	for ( n=0,ln = linelist; ln; ln=ln->next, n++ ) {

		inp[n] = scan_input( ln->line, &nc );

		if ( n == 0 )
			nc0 = nc;

		if ( nc != nc0 ) {
			printf( "Error. Inconsistent number of inputs %i %i in %s\n", 
				nc, nc0, ln->line );
			exit( 1 );
		}
	}

	longlinelist_free( linelist );

	*ncolumn = nc;
	*ndata = n;

	return( inp );
}

void	syn_print( SYNAPS *syn, int cycle, char *filename, char *comment )

{
	FILE	*fp;
	int	l,i,k,nc;

	if ( ( fp = fopen( filename, "w" ) ) == NULL ) {
	        printf( "Cannot open file %s\n", filename );
                exit( 1 );
        }

	fprintf( fp, "TESTRUNID %s\n", comment );

	for ( l=0; l<syn->nlay; l++ )
		fprintf( fp, "%8i  LAYER: %4i\n", syn->nnlay[l], l+1 );

	fprintf( fp, "%8i   :ILEARN\n", cycle );

	nc = 0;

	for ( l = 0; l < syn->nlay - 1; l++) {

		for (k = 0; k < syn->nnlay[l + 1]; k++)
                for (i = 0; i < syn->nnlay[l] + 1; i++) {
                        fprintf(fp, "%13f", syn->wgt[l][k][i]);

			nc++;

			if ( nc == 5 ) {
				fprintf(fp, "\n" );
				nc = 0;
			}
		}

		
        }

	if ( nc != 0 )
		fprintf(fp, "\n" );

	fclose( fp );
}

void	print_prediction( char *filename, float *p, float *target, int n )

{
	int	i;
	FILE	*fp;

        if ( ( fp = fopen( filename, "w" ) ) == NULL ) {
                printf( "Cannot open file %s\n", filename );
                exit( 1 );
        }

	for ( i=0;i<n; i++ )
		fprintf( fp, "%i Pred %f Target %f\n", i+1, p[i], target[i] );

	fclose( fp );
}

main( int argc, char *argv[] )

{

	SYNAPS	*syn;
	int	i;
	float	**x;
	float	**inp;
	float	*target, *p;
	int	nout, nin, n, nc;
	int	ncy;
	LINE	comment;
	float	**t_inp;
	int	t_n, t_nc;
	float	*t_target, *t_p;
	float   t_err, best_t_err;
        float   t_pcc, best_t_pcc;
        float   err, pcc;
	int	*order, ix;
	float	**dv, *dw;

	pparse( &argc, &argv, param, 1, "inputfile" );

	if ( p_nhid < 1 ) {
		printf( "Error. NNLIN can only run with number of hidden neurons > 0. Exit\n" );
		exit( 1 );
	}

	if ( p_seed >= 0 )
		setseed( p_seed );

	/* Read input from file argv[1] */
	/* nc is number of inputs per line and n is number of input lines */
	inp = init_input( argv[1], &nc, &n );

	/* Number of output neurons is 1 */
	/* Number of input neurons is nc - 1 */
	nout = 1;
	nin = nc - 1;

	/* target is the measured target value */
	target = fvector( 0, n-1 );

	for ( i=0; i<n; i++ )
		target[i] = inp[i][nin];

	/* If test file is defined, read test data */
        if ( strlen( p_testfile ) > 0 ) {

		/* read test data */
                t_inp = init_input( p_testfile, &t_nc, &t_n );

		/* allocate vectors to store test targets and prediction values */
                t_target = fvector( 0, t_n-1 );
                t_p = fvector( 0, t_n-1 );

		/* Test for consistency between test and training data */
                if ( t_nc != nc ) {
                        printf( "Test data and training data are inconsistent %i %i\n", nc, t_nc );
                        exit( 1 );
                }
                
		/* Store test target values */
                for ( i=0; i<t_n; i++ )
                        t_target[i] = t_inp[i][nin];
        }

	if ((syn = synaps_alloc()) == NULL) {
                printf( "Error. Cannot allocate SYNAPS\n");
                exit(1);
        }

	/* define synapse structure */
	/* 
	Number of layers is 3 always!
	*/

	syn->nlay = 3;
	syn->nnlay = ivector( 0, syn->nlay-1 );

	syn->nnlay[0] = nin;
	syn->nnlay[1] = p_nhid;
	syn->nnlay[2] = nout;

	syn->maxnper = ( nin > p_nhid ? nin : p_nhid );

	printf( "# Ninput %i. Noutput %i\n", nin, nout );

	for ( i=0; i<syn->nlay; i++ )
		printf( "# Nlayer %i NNlay %i\n", i, syn->nnlay[i] );

	/* 
	define x matrix
	define vector and matrix to store changes in synapses
	*/

	x = fmatrix( 0, 2, 0, syn->maxnper );
	dw = fvector( 0, syn->nnlay[1] );
	dv = fmatrix( 0, syn->nnlay[1], 0, syn->nnlay[0] );

	/* Initialize synapses */
	syn->wgt = init_wgt( syn );

	/* define vector to store network output */
	p = fvector( 0, n-1 );

	best_t_err = 9999.9;
	best_t_pcc = -9999.9;

	/* define vector with order of training samples */
	order = ivector_ramp( 0, n-1 );

	for ( ncy=0;ncy<p_ncycles; ncy++ ) {

		/* rerandomize order of training samples */
		ivector_rerandomize( order, 0, n-1 );

		for ( i = 0; i < n; i++ ) {
			
			/* select training sample */
			ix = order[i];

			/* do feed forward */
			forward( x, syn, inp[ix] );

			/* store prediction */
			p[ix] = x[2][0];

			/* Check size of error */
			/* If error greater than p_bplim to backpropagation */
                        //SQR(誤差)が 0.00001 以上であれば、
			if ( SQR( p[ix] - target[ix] ) > p_bplim ) 
				backprop( x, syn, inp[ix], target[ix], dw, dv );

		}

		/* Test set performance is evaluated every p_ntest cyles */
		if ( ( ncy % p_ntest ) == 0 ) {

			/* Get train set pearsons correlation and mean squared error */
			pcc = fvector_xycorr( n, p, target );
                	err = fvector_xyerror( n, p, target );

			/* Get test set predictions */
                        for ( i = 0; i < t_n; i++ ) {
                                                     
                                forward( x, syn, t_inp[i] );
                                                            
                                t_p[i] = x[2][0];        
                        }                           
                        
			/* Get test set pearsons correlation and mean squared error */                            
                        t_err = fvector_xyerror( t_n, t_p, t_target );
                        t_pcc = fvector_xycorr( t_n, t_p, t_target );   
                                                                          
			/* Check if test error (for p_dtype == 0) or test pcc (for p_dtype != 0 ) has improved */ 
                        if ( ( p_dtype == 0 && t_err < best_t_err ) || ( ( p_dtype > 0 ) && t_pcc > best_t_pcc ) ) {

                                best_t_err = t_err;
                                best_t_pcc  = t_pcc;
                                     
                                sprintf( comment, "Ncycles %i Train Pearson %f Err %f Test Pearson %f Err %f",
                                        ncy, pcc, err, t_pcc, t_err );
                                
				/* Dump synapses to file */       
                                syn_print( syn, ncy, p_syn, comment );
                                     
                        }            

			printf( "Train epoch %i N Train %i Pearson %f Err %f Test Pearson %f Err %f\n", 
				ncy, n, pcc, err, t_pcc, t_err );

                }   
		
	}

	/* Read in best synapses from file */
        syn = synaps_read_file( p_syn );
                                    
	/* If file for train set prediction is defined, run train predicitons and save to file */ 
        if ( strlen( p_trainpred ) ) {
                for ( i = 0; i < n; i++ ) {
                        forward( x, syn, inp[i] );
                                     
                        p[i] = x[2][0];
                }                    
                                     
                print_prediction( p_trainpred, p, target, n );
        }                            
        
	/* If file for test set prediction is defined, run test predicitons and save to file */                             
        if ( strlen( p_testpred ) ) {
                for ( i = 0; i < t_n; i++ ) {
                                     
                        forward( x, syn, t_inp[i] );
                                     
                        t_p[i] = x[2][0];
                }                    
                                     
                print_prediction( p_testpred, t_p, t_target, t_n );
        }              

	exit( 0 );
}
