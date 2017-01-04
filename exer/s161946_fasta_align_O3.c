/* M.Nielsen May 2008, mniel@cbs.dtu.dk */

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

float   p_gapf; 
float   p_gapn;
int     p_minall;
FILENAME	p_blmat;
int	p_verbose;
int	p_show;

PARAM   param[] =
{
	"-gf", VFLOAT p_gapf, "penalty for first gap", "11",
	"-gn", VFLOAT p_gapn, "penalty for next gap", "1",
	"-mal", VINT p_minall, "Minimum alignment lenght", "3",
	"-m", VFNAME	p_blmat, "Scoring matrix", "$ALGOHOME/data/BLOSUM50",
	"-v", VSWITCH	p_verbose, "Verbose mode", "0",
	"-show", VSWITCH p_show, "Show matrices", "0",
	0
};

/* Smith Watermann alignment algorithm */

float	sw_alignment( float **m,		/* Scoring matrix, called d in lecture */
	int l1,			/* Length of query sequence */
	int l2,			/* Length of database sequence */
	float fg,		/* Penalty for first gap */
	float ng,		/* Penalty for each of subsequent gaps */
	float **S,		/* match scores, called D is lecture */
	int *firsti,		/* Offset in query sequence */
	int *firstj,		/* Offset in database sequence */
	char *qal,		/* Query alignment */
	char *dal, 		/* DB algnemnt */
	int *alen,		/* alignemnt length */
	char *qseq,		/* Query sequence */
	char *dseq		/* Database sequence */
	)

{
	int	**E; /* Path matrix called e in lecture */
	int     i, j;
	float   temp1, temp2, temp;
	float   sij, dij, qij;
	float	score;
	int	keep_going;
	int	best;
	int	k;
	int	best_i, best_j;

	score = 0;
	/* firsti and firstj will store the start location of the alignment */
	(*firsti) = -1;
	(*firstj) = -1;

	E = imatrix(0, l1, 0, l2);

	/* Loop over Query sequence */
	for (i = l1 - 1; i >= 0; i--) {

		for (j = l2 - 1; j >= 0; j--) {

/* HERE YOU MUST FILL OUT THE MISSING CODE (XXXXX)*/

			/* We must find the best of five possible scores
				1) match
				2) gap opening query 
				3) gap extension query
				4) gap opening in db
				5) gap extension in db

			*/

			/* Try match state */

                      //sij = XXXXX
                      sij = S[i+1][j+1] + m[i][j]; 
                      

			/* Try gap in Query sequence (insertion in database) */

                      //temp1 = XXXXX /* Open a gap */
                      temp1 = S[i][j+1] - fg;                      

			temp2 = 0;
			best_j = -9; /* best_j keeps track of the lenght of the gap */
			/* best_j > 0 will indicate a gap extension, and best_j = -9 
			that a gap opening was used  */
			for ( k=j+2; k<l2; k++ ) { /* Gap extension */

				temp = S[i][k] - fg - (k-j-1)*ng;

				if ( temp > temp2 ) {
					temp2 = temp;
					best_j = k;
				}
			}

			if ( temp1 >= temp2 ) { /* gap opening best */
				qij = temp1;
				best_j = -9;
			}
			else { /* gap extension best */
				qij = temp2; 
			}

			/* Try gap in database sequence (insertion in Query) */

			//temp1 = XXXXX; /* Open a gap */
                        temp1 = S[i+1][j] - fg;
                        
			temp2 = 0;
			best_i = -9; /* best_i keeps track of the lenght of the gap */
			for ( k=i+2; k<l1; k++ ) { /* Gap extension */

				/* And some more code selection the best cap extension */
                              temp = S[k][j] - fg - (k-i-1)*ng;

                              if ( temp > temp2 ) {
                                    temp2 = temp;
                                    best_i = k;
                              }
			}

			/* Some code determining if gap opening or gap extension is best. 
			Set dij to the best of these two and set best_i according (see above) */

			if ( temp1 >= temp2 ) { /* gap opening best */
                              dij = temp1;
                              best_i = -9;
			}
			else { /* gap extension best */
                              dij = temp2; 
			}
                        

 			/* eij is the backtrack direction matrix
                       		eij = 0 stop back tracking
                        	eij = 1 match
                        	eij = 2 gap-opening database
                        	eij = 3 gap-extension database
                        	eij = 4 gap-opening query
                        	eij = 5 gap-extension query
                	*/

			if ( sij >= qij && sij >= dij ) { /* Match is best */

				E[i][j] = 1;

			}
			else if ( qij > dij ) { /* Gap in query (insertion in DB) is best */

				sij = qij;
				if ( best_j > 0 ) /* Gap extension */
					E[i][j] = 5;
				else
					E[i][j] = 4;

			}
			else { /* Gap in database (insertion in Query) is best */

                              sij = dij;
				if ( best_i > 0 ) /* Gap extension */
					E[i][j] = 3;
				else
					E[i][j] = 2;

			}

			if (sij > score) {
				score = sij;
				(*firsti) = i;
				(*firstj) = j;
			}

			if ( sij <= 0 ) {
				sij = 0.0;
				E[i][j] = 0;
			}

			S[i][j] = sij;
		}
	}

	if ( p_show ) {
		printf( "# S-matrix\n" );

             	printf( "  " );
       		for ( j=0; j<l2; j++ )
                	printf( "%5c ", dseq[j] );
        	printf( "\n" );
                       
        	for ( i=0; i<=l1; i++ ) {
                	printf( "%c", ( i<l1 ? qseq[i] : ' ' ));
                	for ( j=0; j<=l2; j++ )
                        	printf( " %5.2f", S[i][j] );
                	printf( "\n" );
        	}      

		printf( "# Eij-matrix\n" );
                           
		printf( "  " );
        	for ( j=0; j<l2; j++ )
                	printf( "%2c ", dseq[j] );
        	printf( "\n" );
                          
        	for ( i=0; i<=l1; i++ ) {
                	printf( "%c", ( i<l1 ? qseq[i] : ' ' ));
                	for ( j=0; j<=l2; j++ )
                        	printf( " %2i", E[i][j] );
                	printf( "\n" );
        	}     
	}

	/* Do back tracking */

	if (*firsti < 0 || *firstj < 0 ) 
		printf( "No alignment found. Exit\n" );

	(*alen) = 0;

	i = *firsti;
	j = *firstj;
	qal[(*alen)] = qseq[i];
	dal[(*alen)] = dseq[j];
	i++;
	j++;

	(*alen)++;

	keep_going = 1;

#define SMALL 1e-10

	while ((i < l1 ) && (j < l2 ) && keep_going ) {


		if ( E[i][j] == 0 ) {
			keep_going = 0;
		} 
		else if ( E[i][j] == 1 ) { /* Match */

			qal[(*alen)] = qseq[i];
			dal[(*alen)] = dseq[j];
			i++;
			j++;
			(*alen)++;

		}
		else if ( E[i][j] == 4 ) { /* gap opening in Query */

			qal[(*alen)] = '-';
			dal[(*alen)] = dseq[j];
			j++;
			(*alen)++;

		}
		else if (  E[i][j] == 5 ) { /* gap extension in Query */

                      best = j+2;
                      temp = S[i][best] - fg - ( best-j-1) * ng - S[i][j];

                      while ( temp*temp > SMALL ) {
                            best++;
                            temp = S[i][best] - fg - (best-j-1) * ng - S[i][j];
                            //printf("%f\n",S[i][best]);
                            
                      }

                      for ( k=j; k<best; k++ ) {
                            qal[(*alen)] = '-';
                            dal[(*alen)] = dseq[j];
                            j++;
                            (*alen)++;
                      }
                      
		}
		else if ( E[i][j] == 2 ) { /* gap opening in Database */

			qal[(*alen)] = qseq[i];
			dal[(*alen)] = '-';
			i++;
			(*alen)++;

		}
		else if (  E[i][j] == 3 ) { /* gap extension in Database */

			/* Write some code to find the best gap extension (see above code for query) */

			best = i+2;
			temp = S[best][j] - fg - ( best-i-1) * ng - S[i][j];

			while ( temp*temp > SMALL ) {
				best++;
				temp = S[best][j] - fg - ( best-i-1) * ng - S[i][j];
			}

			for ( k=i; k<best; k++ ) {
				qal[(*alen)] = qseq[i];
				dal[(*alen)] = '-';
				i++;
				(*alen)++;
			}
                        
                        
                        
		}
	}
        
/* THIS IS THE END */
        
	qal[(*alen)] = 0;
	dal[(*alen)] = 0;
        
	imatrix_free( E, 0, l1, 0, l2);
        
	return( score );
}

ALN    *align(float **m, FSALIST *q, FSALIST *d, float gapf, float gapn )

{
	/* m is the score matrix called d in the lecture */
	float **sco;
	float   score;
	int     firsti, firstj;
	int     alen;
	ALN    *new = NULL;
	int     i;
	char	qpal[15000], dpal[15000];

	sco = fmatrix(0, q->len, 0, d->len); /* D matrix from lecture */

	score = sw_alignment(m, q->len, d->len, gapf, gapn, sco, &firsti, &firstj, 
		qpal, dpal, &alen, q->seq, d->seq );

	fmatrix_free(sco, 0, q->len, 0, d->len);

	if (alen < p_minall)
		return( NULL );

	/* Allocate alignment structure to store alignment */

	if ((new = aln_alloc()) == NULL) {
		printf("Cannot alloc new ALN\n");
		exit(1);
	}

	new->alen = alen;

	new->score = score;

	new->mlen = 0;
	new->ngap = 0;
	new->nid = 0;

	for (i = 0; i < new->alen; i++) {
		if (qpal[i] != '-' && dpal[i] != '-') {
			new->mlen++;
		} else {
			new->ngap++;
		}
		if (qpal[i] == dpal[i])
			new->nid++;
	}

	new->qof = firsti;
	new->qal = cvector(0, strlen(qpal));
	strcpy(new->qal, qpal);

	new->dof = firstj;
	new->dal = cvector(0, strlen(dpal));
	strcpy(new->dal, dpal);

	strcpy(new->qname, q->name);
	new->qlen = q->len;

	strcpy(new->dname, d->name);
	new->dlen = d->len;

	strcpy(new->type, "SW_ALN");

	new->rscore = -new->score;

	return (new);

}

float	**score_mat( FSALIST *q, FSALIST *d, float **blmat )

{
	float	**scomat;
	int	i,j;
	int	ix, jx;

	/* Allocate the matrix */
	scomat = fmatrix( 0, q->len-1, 0, d->len-1 );

	for ( i=0; i<q->len; i++ ) {

		ix = q->i[i];

		if ( ix < 0 ) {
			printf( "Error. Unconventional amino acid i query sequence %c %s\n", q->seq[i], q->name );
			exit( 1 );
		}

		if ( ix > 19 )
			continue;

		for ( j=0; j<d->len; j++ ) {

			jx = d->i[j];

			if ( jx < 0 ) {
				printf( "Error. Unconventional amino acid i query sequence %c %s\n", q->seq[i], q->name );
				exit( 1 );
			}

			if ( jx > 19 )
				continue;

			scomat[i][j] = blmat[ix][jx];
		}
	}

	return( scomat );

}

main(int argc, char *argv[])

{
	FSALIST		*q_fsa, *db_fsa, *d;
	float   	gapf, gapn;
	ALN		*new;
	float		**blmat, **scomat;
	char		*alphabet;

	pparse(&argc, &argv, param, 2, "fsa1 db");

	/* Read Blosum substutution scoring matrix from file */

	blmat = read_blosummat( p_blmat, &alphabet );

        if ( ! blmat ) {
                printf( "Error. Cannot read BLOSUM matrix from file %s. Exit\n", p_blmat );
                exit( 1 );
        }


	/* Read query FASTA file */
	if ( ( q_fsa = fsalist_read( argv[1] ) ) == NULL ) {
		printf("Cannot read fasta file %s\n", argv[1] );
		exit(1);
	}

	q_fsa = fsalist_check_names( q_fsa );

	/* Assign Blosum alphabet order to variable q_fsa->i */
	
	fsalist_iassign_profile_order( q_fsa );

	/* Read database FASTA file */

	if ( ( db_fsa = fsalist_read( argv[2] ) ) == NULL ) {
                printf("Cannot read fasta file %s\n", argv[2] );
                exit(1);
        }

	db_fsa = fsalist_check_names( db_fsa );

	/* Assign Blosum alphabet order to variable ->i in all faste entries in db_fsa */

	fsalist_iassign_profile_order( db_fsa );

	gapf = p_gapf; //penalty for opening gap
	gapn = p_gapn; //penalty for extension gap

	printf("# Gap penalties. fgap: %f. ngap: %f\n", gapf, gapn);

	for ( d = db_fsa; d; d=d->next ) { //pointerを次に進めるときは->next

		/* Make scoring matrix as described in exercise 1 */
		/* This is the d matrix from the lecture */

		/* The score_mat routine returns a matrix of size q_fsa->len, d->len with the blosum scores for matching i against j */
		scomat = score_mat( q_fsa, d, blmat );

		/* Calculate alignment */

		new = align( scomat, q_fsa, d, gapf, gapn);

		if ( new )  /* Print alignment to screen */
			 aln_write_single( new );

		/* Free dymanically allocated memory */

		fmatrix_free( scomat, 0, q_fsa->len-1, 0, d->len-1 );
		aln_free( new );
	}

	exit(0);
}
