/* mniel@cbs.dtu.dk November 2008 */
/* Code strongly inspired by the umdhmm-v1.02 code by Tapas Kanungo, kanungo@cfar.umd.edu */

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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"

FILENAME	p_alphabet;
int		p_printhmm;
int		p_printmat;

PARAM   param[] = {
	"-a", VFILENAME	p_alphabet, "Symbol alphabet", "123456",
	"-phmm", VSWITCH p_printhmm, "Print the HMM", "0",
	"-pmat", VSWITCH p_printmat, "Print Viterbi matrix", "0",
        0
};

typedef struct {
        int N;          /* number of states;  Q={0,1,...,N-1} */
        int M;          /* number of observation symbols; V={0,1,...,M-1}*/
        double  **A;    /* A[0..N-1][0..N-1]. a[i][j] is the transition prob
                           of going from state i at time t to state j
                           at time t+1 */
        double  **B;    /* B[0..N-1][0..M-1]. b[j][k] is the probability of
                           of observing symbol k in state j */
        double  *pi;    /* pi[0..N-1] pi[i] is the initial state distribution. */
} HMM;

HMM	*hmm_alloc()

{

	HMM	*n;

	if ( ( n = ( HMM * ) malloc ( sizeof(HMM))) != NULL ) {
		n->A = NULL;
		n->B = NULL;
		n->pi = NULL;
		n->M = -99;
		n->N = -99;
        }

        return( n );
}

void	hmm_free( HMM *hmm )

{

	dmatrix_free( hmm->A, 0, hmm->N-1, 0, hmm->N-1);
        dmatrix_free( hmm->B, 0, hmm->N-1, 0, hmm->M-1);
        dvector_free( hmm->pi, 0, hmm->N-1);

	free( hmm );
}

HMM	*hmm_read( char *filename )

{
	HMM	*hmm;
	int i, j;
	FILE	*fp;

	if ( (fp = fopen( filename, "r" ) ) == NULL ) {
		printf( "Error. Cannot open file %s. Exit\n", filename );
		exit( 1 );
	}	

	hmm = hmm_alloc();

	if ( ! hmm ) {
		printf( "Error. Cannot allocate HMM. Exit\n" );
		exit( 1 );
	}

	if ( fscanf( fp, "M= %d\n", &(hmm->M)) != 1 ) {
		printf( "Error. Wrong File format M= \n" );
		exit( 1 );
	}

	if ( fscanf( fp, "N= %d\n", &(hmm->N)) != 1 ) {
		printf( "Error. Wrong File format N= \n" );
		exit( 1 );
	}

	fscanf( fp, "A:\n" );

	hmm->A = dmatrix( 0, hmm->N - 1, 0, hmm->N -1 );
	for ( i = 0; i<hmm->N; i++) {
                for (j = 0; j < hmm->N; j++) {
                        fscanf(fp, "%lf", &(hmm->A[i][j]));
                }
                fscanf(fp,"\n");
        }

	fscanf( fp, "B:\n" );

        hmm->B = dmatrix( 0, hmm->N - 1, 0, hmm->M -1 );
        for ( i = 0; i<hmm->N; i++) {
                for (j = 0; j < hmm->M; j++) {
                        fscanf(fp, "%lf", &(hmm->B[i][j]));
                }
                fscanf(fp,"\n");
        }

	fscanf( fp, "pi:\n" );

        hmm->pi = dvector( 0, hmm->N - 1 );
	for ( i = 0; i<hmm->N; i++)
		fscanf(fp, "%lf", &(hmm->pi[i]));

	fclose( fp );

	return( hmm );
}

void	hmm_print( HMM *hmm )

{
	int	i,j;

	printf( "M= %d\n", hmm->M );
	printf( "N= %d\n", hmm->N );

	printf( "A:\n" );
	for ( i=0; i<hmm->N; i++ ) {
		for ( j=0; j<hmm->N-1; j++ )
			printf( "%6.4f ", hmm->A[i][j] );

		printf( "%6.4f\n", hmm->A[i][hmm->N - 1] );
	}

	printf( "B:\n" );
	for ( i=0; i<hmm->N; i++ ) {
                for ( j=0; j<hmm->M-1; j++ )
                        printf( "%6.4f ", hmm->B[i][j] );

                printf( "%6.4f\n", hmm->B[i][hmm->M - 1] );
        }

	printf( "pi:\n" );
	for ( i=0; i<hmm->N-1; i++ ) 
		printf( "%6.4f ", hmm->pi[i] );
	printf( "%6.4f\n", hmm->pi[hmm->N - 1] );

}

char	*sequence_read( char *filename, int *T )

{
	LINELIST        *linelist, *ln;
        int     i;
        char    tvec[10000], *cvec;

	linelist = linelist_read( filename );

        if ( ! linelist ) {
                printf( "Error. Cannot read Sequence from file %s\n", filename );
                exit( 1 );
        }

        *T = 0;

        for ( ln=linelist; ln; ln=ln->next ) {

		for ( i=0; i<strlen(ln->line); i++ ) 
			tvec[(*T)++] = (ln->line)[i];

		if ( *T >= 10000 ) {
			printf( "Error. Sequence too long %i maxsize 10000\n", *T );
			exit( 1 );
		}
	}

	cvec = cvector( 0, *T );

	strncpy( cvec, tvec, *T );

	cvec[*T] = 0;

	return( cvec );
}

int	*seq2index( char *seq, int len )

{

	int	*ivec, i, ix;

	ivec = ivector( 0, len-1 );

	for ( i=0; i<len; i++ ) {
		ix = strpos( p_alphabet, seq[i] );

                if ( ix < 0 ) {
                	printf( "Error. Symbol %c not in alphabet %s\n", seq[i], p_alphabet );
                        exit( 1 );
                }

 	  	ivec[i] = ix;
	}

	return( ivec );
}

double	viterbi( HMM *hmm, int T, int *O, double **delta, int **psi, int *q )

{
        int     i, j;   /* state indices */
        int     t;      /* time index */

        int     maxvalind;
        double  maxval, val, pprob;

	/* The viterbi recursion

		delta[t][j] = p(t,j) * max ( a_ji * delta[t-1][i] )

	where p(t,j) == hmm->B[j][O[t]] is the probability that the state j emits the symbol O[t] and
	a_ji == hmm->A[i][j] is the transition probability from state i to j. 

	psi[t][j] holds the state index of the selected transition from time t-1 to t into state j

	*/

	/* 1. Initialization  */

        for (i = 0; i < hmm->N; i++) {
              //delta[0][i] = XXXXXXX;
              //psi : 初期値 , B[i][O[0]] : iのサイコロを使用して、symbol O[0]が出る確率
              delta[0][i] = hmm->pi[i]*(hmm->B[i][O[0]]);
              psi[0][i] = 0;
        }

	 /* 2. Recursion */
        //T : number of symbol read
        //N : number of state
        for (t = 1; t < T; t++) {
                for (j = 0; j < hmm->N; j++) {
                        maxval = 0.0;
                        maxvalind = 0;
                        //for ( XXXXXX ) {
                        int k; //変数kの使用
                        for ( k = 0; k < hmm->N; k++ ) {
                              //val = XXXX ;
                              //t-1番目の状態のstate kの値 k→j(全てのjに変わる変化を記録しておくため)
                              val = (hmm->B[j][O[t]])*delta[t-1][k]*hmm->A[k][j];
                              
                              if (val > maxval) {
                                      /* maxval = XXXXX; */
                                      /* maxvalind = XXXXX; */
                                      maxval = val;
                                      maxvalind = k; 
                                      
                                }
                        }

                        //delta[t][j] = XXXX;
                        delta[t][j] = maxval;
                        psi[t][j] = maxvalind;

                }
        }

	/* 3. Termination */
	/* Find highest scoring cell in Viterbi matrix */
	/* q is a string that holds the selected states */

        pprob = 0.0;
        q[T-1] = 0;
        for (i = 0; i < hmm->N; i++) {
                if (delta[T-1][i] > pprob) {
                      //pprob = XXXXXXX;
                      pprob = delta[T-1][i];
                      q[T-1] = i;
                }
        }

        /* 4. Path (state sequence) backtracking */

        for (t = T - 2; t >= 0; t--)
                q[t] = psi[t+1][q[t+1]];

	return( pprob );

}

void	sequence_print( int T, int *O )

{
	int	i;

	for ( i=0; i<T; i++ )
		printf( "%c ", p_alphabet[O[i]] );

	printf( "\n" );

}

main(int argc, char **argv)

{
	int 	T; 
	HMM  	*hmm;
	char	*seq;	/* Observed sequence */
	int	*O;	/* observation index sequence O[0..T-1] */
	int	*q;	/* state sequence q[0..T-1] */
	double **delta;
	int	**psi;
	double 	prob; 
	int	i,j;

	pparse( &argc, &argv, param, 2, "hmm sequence" );

	/* Read HMM model */

	hmm = hmm_read( argv[1] );

	if (hmm == NULL) {
		printf( "Error: Cannot reda HMM from file %s\n", argv[1]);
		exit (1);
	}

	printf( "# HMM read from file %s\n", argv[1] );

	if ( p_printhmm )
		hmm_print( hmm );

	/* Read observations */
	/* T holds the number of symbols read */

	seq = sequence_read(argv[2], &T );

	if ( seq == NULL) {
		printf( "Error: Cannot read Sequence data from %s\n", argv[2]);
		exit (1);
	}

	printf( "# Sequence read from file %s. Number of symbols %i\n",
		argv[2], T );

	/* Transform observation symbols to index values defined by the alphabet */
	/* Here 1 -> 0, 2 -> 1 etc */
	/* This is just like we transformed amino acids into number between 0 and 19 */

	O = seq2index( seq, T );

	/* Allocate vector for state path */
	q = ivector( 0, T-1 );

	/* Allocate matrix for Viterbi */
	delta = dmatrix( 0, T-1, 0, hmm->N - 1);

	/* Allocate vector for backtracking */
	psi = imatrix( 0, T-1, 0, hmm->N - 1);

	printf("# ------------------------------------\n");
	printf("# Viterbi using direct probabilities\n");

	/* Find most likely path using Viterbi */
	prob = viterbi( hmm, T, O, delta, psi, q); 

	/* Print Viterbi matrix (for debugging) */
	if ( p_printmat ) {
		for ( i=0; i< hmm->N; i++ ) {
                	printf( "%i", i );
                	for ( j=0; j<T; j++ ) {
                        	printf( " %f", log(delta[j][i])/log(10 ) );
                	}
                	printf( "\n" );
        	}
	}

	printf( "Viterbi  MLE log prob = %E\n", log(prob));

	printf( "# Sequence: " );
	sequence_print( T, O );

	printf( "Optimal state sequence:\n");

	sequence_print( T, q);

	ivector_free( q, 0, T-1);
	ivector_free( O, 0, T-1);
	imatrix_free( psi, 0, T-1, 0, hmm->N - 1);
	dmatrix_free( delta, 0, T-1, 0, hmm->N - 1);
	hmm_free( hmm );
}
