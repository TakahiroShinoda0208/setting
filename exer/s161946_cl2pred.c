/* M.Nielsen April 2002 mniel@cbs.dtu.dk */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

FILENAME	p_matrix;
WORD    	p_weights;

PARAM   param[] = {
	"-mat", VFNAME p_matrix, "Matrix file", "my.mat",
	"-w", VWORD     p_weights, "Weigths on positions", "",
	0
};

main( int argc, char *argv[] )

{
	int		c;
	int		nn;
	int		len, mlen;
	char		*alphabet;
	float		s, max_s;
	float		**matrix;
	int		i,j;
	WORD		core;
	int		best_pos;
	float		*weights=NULL;
	int		nw;
	PEPLIST		*peplist, *pl;

	pparse( &argc, &argv, param, 1, "input_file" );

	/* Read PSSM, and the scoring amino acids alphabet, and the motif length */
	matrix = hmm_matrix_read_l( p_matrix, &alphabet, &mlen );

	/* Read the peptide list */
	peplist = peplist_read( argv[1] );

	/* Process the position specific weights */
	if ( strlen( p_weights ) > 0 ) {

		weights = fvector_split( p_weights, &nw );

		printf( "# Number of weigths found in list %i\n", nw );

		if ( nw != mlen ) {
			printf( "Error. nw %i != l %i. Exit\n",
				nw, mlen );
			exit( 1 );
		}

	}

	nn = 0;

	for ( pl=peplist; pl; pl=pl->next ) {

		len = pl->len;

		if ( len < mlen ) {
			printf( "# Length of %s is %i shorter than motif length %i\n",
				pl->pep, len, mlen );
			continue;
		}

		max_s = -9999;

		/* Loop over all submers of length mlen */
		/* to find highest scoring submer */
		for ( i=0; i<len-mlen+1; i++ ) {

			/* Find score of submer */
			s = 0.0;

			for ( j=0; j<mlen; j++ ) {

                              //長さlenのペプチドを長さmlenで切り出して、loopを回して 一番高いスコアを探している。
                              // mappingで少しずつ進める感じ

				/* Find the amino acid type at position i+j */
				c = strpos( alphabet, pl->pep[i+j] );

				/* the syntax 
					s += mat * ( weights ? weights[j] : 1.0 ) 

				is a compact way of writing

				if ( weights != NULL )
					s += mat * weights[j];
				else
					s += mat * 1;
				*/

				if ( c >= 0 ) {
                                      //s += p_matrix[j][c] * ( weights ? weights[j] : 1.0 );
                                      s += matrix[j][c] * ( weights ? weights[j] : 1.0 );
				}
			}

			if ( s > max_s ) {
				max_s = s;
				best_pos = i;
			}

		}	

		strncpy( core, pl->pep+best_pos, mlen );
		core[mlen] = 0;
		printf( "%4i %4i %12s %8.4f %s\n", nn, best_pos, core, max_s, pl->line );

		nn++;

	}

	printf( "# Number of peptides %i\n", nn );

	exit( 0 );
}
