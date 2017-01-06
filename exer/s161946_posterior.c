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
int		p_state;

PARAM   param[] = {
	"-a", VFILENAME	p_alphabet, "Symbol alphabet", "123456",
	"-phmm", VSWITCH p_printhmm, "Print the HMM", "0",
	"-s", VINT	p_state, "State for posterior decoding", "0",
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
	int 	i, j;
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

#define MAXSIZE 10000

char	*sequence_read( char *filename, int *T )

{
	LINELIST        *linelist, *ln;
        int     i;
        char    tvec[MAXSIZE], *cvec;

	linelist = linelist_read( filename );

        if ( ! linelist ) {
                printf( "Error. Cannot read Sequence from file %s\n", filename );
                exit( 1 );
        }

        *T = 0;

        for ( ln=linelist; ln; ln=ln->next ) {

		for ( i=0; i<strlen(ln->line); i++ ) 
			tvec[(*T)++] = (ln->line)[i];

		if ( *T >= MAXSIZE ) {
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

void	sequence_print( int T, int *O )

{
	int	i;

	for ( i=0; i<T; i++ )
		printf( "%c ", p_alphabet[O[i]] );

	printf( "\n" );

}

void	forward( HMM *hmm, int T, int *O, double **alpha )

{
        int     i, j;  
        int     t; 
        double  sum;

	/* Forward recursion is

		alpha[t+1][j] = p(t+1,j) * sum ( alpha[t][i] * a_ij )

	where p(t+1,j) == hmm->B[j][O[t+1]] is the probability that the state j emits the symbol O[t+1] and
        a_ij == hmm->A[i][j] is the transition probability from state i to j.

        */

	/* 1. Initialization  */

        for (i = 0; i < hmm->N; i++) 
		alpha[0][i] = hmm->pi[i]* hmm->B[i][O[0]];

	 /* 2. Recursion */

        for (t = 1; t < T; t++) {
                for (j = 0; j < hmm->N; j++) {
			sum = 0.0;
                        for (i = 0; i < hmm->N; i++){
                              //sum += XXXXXX;
                              sum+= alpha[t-1][i]*hmm->A[i][j];
                        }
                        
                        
                        //alpha[t][j] = sum * XXXX;
                        alpha[t][j] = sum * hmm->B[j][O[t]];
                        
                }
        }

}

void backward( HMM *hmm, int T, int *O, double **beta )

{

	int     i, j, t;
        double  sum;

	/* Backward recursion is

		beta[t][i] = sum ( p(t+1,j) * beta[t+1][j] * a_ij )

	where p(t+1,j) == hmm->B[j][O[t+1]] is the probability that the state i emits the symbol O[t+1] and
        a_ij == hmm->A[i][j] is the transition probability from state i to j

        */

        /* 1. Initialization */

        for (i = 0; i < hmm->N; i++)
                beta[T-1][i] = 1.0;

        /* 2. Recursion */

        for (t = T - 2; t >= 0; t-- ) {

                for (i = 0; i < hmm->N; i++) {
                        sum = 0.0;
                        for (j = 0; j < hmm->N; j++)
                        {
                              //sum += XXXX;
                              sum += hmm->B[j][O[t+1]]*beta[t+1][j]*(hmm->A[i][j]);
                        }
                        
                        //beta[t][i] = XXXX;                        
                        beta[t][i] = sum;
                        
                }
                
        }
}


main(int argc, char **argv)

{
	int 	T; 
	HMM  	*hmm;
	char	*seq;
	int	*O;	
	double **alpha, **beta;
	double 	px, fki, bki, logpost, logpx; 
	int	i;

	pparse( &argc, &argv, param, 2, "hmm sequence" );

	hmm = hmm_read( argv[1] );

	if (hmm == NULL) {
		printf( "Error: Cannot reda HMM from file %s\n", argv[1]);
		exit (1);
	}

	printf( "# HMM read from file %s\n", argv[1] );

	if ( p_printhmm )
		hmm_print( hmm );

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

	alpha = dmatrix( 0, T-1, 0, hmm->N - 1);
	beta = dmatrix( 0, T-1, 0, hmm->N - 1);

	printf( "# Sequence: " );
	sequence_print( T, O );

	/* Calculate P(x) */

	forward( hmm, T, O, alpha ); 

	/* 
		P(x) = sum alpha[T-1][i]

		where the sum of over all states

	*/

	/* px is the probability of observing the sequence given the model */
	px = 0.0;
	for ( i=0; i < hmm->N; i++ )
        {
              //px += XXXXX;
              px += alpha[T-1][i];
        }
        

#define PLOW 1.0e-40

	if ( px > PLOW )
		logpx = log(px);
	else
		logpx = log(PLOW);

	printf( "# log(P(x)) %f\n", logpx );

	/* Calculate backward matrix Beta */

	backward( hmm, T, O, beta );

	for ( i=0; i<T; i++ ) {

		fki = alpha[i][p_state];
		bki = beta[i][p_state];

		/* Post = log(fki*bki/px) */

                //logpost = XXX + XXX - logpx;
                logpost = log(fki) + log(bki) - logpx;
                


                printf( "Posterior %3i %2i %c %7.3f %7.3f %6.4f\n", i, O[i], p_alphabet[O[i]], log(fki), log(bki), exp(logpost) );
        }

	ivector_free( O, 0, T-1);
	dmatrix_free( alpha, 0, T-1, 0, hmm->N - 1);
	dmatrix_free( beta, 0, T-1, 0, hmm->N - 1);
	hmm_free( hmm );
}
