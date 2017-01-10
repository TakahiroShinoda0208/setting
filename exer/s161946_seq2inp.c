/* M.Nielsen April 2008 mniel@cbs.dtu.dk */

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
#include <string.h>
#include <stdlib.h>
#include "utils.h"

FILENAME	p_blm;
float	p_max;
float	p_min;
float	p_blnorm;
int	p_bl;
int	p_strict;

PARAM   param[] = {
	"-blm", VFNAME p_blm, "Blosum matrix", "$ALGOHOME/data/BLOSUM50",
	"-max",	VFLOAT	p_max, "Max value in sparse encoding", "0.9",
	"-min", VFLOAT  p_min, "Min value in sparse encoding", "0.05",
	"-bln", VFLOAT p_blnorm, "Normalizing factor for blosum score", "5.0",
	"-bl", VSWITCH p_bl, "Use Blosum encoding [default is sparse]", "0",
	"-s", VSWITCH p_strict, "Strict mode. Stop if non standart aa [Default ignore]", "0",
	0
};

main( int argc, char *argv[] )

{
	PEPLIST		*peplist, *pl;
	int		i,j;
	float		**blm;
	int		c;
	int		nlines;
	int		len;
	char		*alphabet;
	char		aa;

	pparse( &argc, &argv, param, 1, "input_file" );

	/* Read peptide list */
	peplist = peplist_read( argv[1] );

	nlines = 0;

	if ( p_bl ) {
		/* Read Blosum matrix */
		blm = read_blosummat( p_blm, &alphabet );
		printf( "# Blosum matrix %s initialized\n", p_blm );

		/* Normalize Blosum matrix */
		for ( i=0; i<strlen(alphabet); i++ )
		for ( j=0; j<strlen(alphabet); j++ )
			blm[i][j] /= ( p_blnorm > 0 ? p_blnorm : 1.0 );
	}
	else {
		/* Define amino acids alphabet */
		alphabet = cvector( 0, strlen( "ARNDCQEGHILKMFPSTWYV" ) );
		strcpy( alphabet, "ARNDCQEGHILKMFPSTWYV" );
	}

	for ( pl=peplist; pl; pl=pl->next ) {

		/* 
		Each peplist element is a structure of the form (see peputils.h file in utils)
		typedef struct peplist  {
        		struct  peplist *next;
        		WORD    pep;
        		int     len;
        		float   score; 
        		int     nn;
        		int     *iv;
        		LINE    line;
		} PEPLIST;

		the score variable contain the peptide affinity
		*/
	
		nlines++;

		len = pl->len;

		if ( len != peplist->len ) {
			printf( "Inconsistent peptide length %s %i != %i. Exit\n", pl->pep, len, peplist->len );
			exit( 1 );
		}

		for ( i=0; i<len; i++ ) {

			aa = pl->pep[i];

			/* Find position of amino acid in alphabet */
			c = strpos( alphabet, aa );

                        //if ( c < 0 && aa != 'X' ) {
                        if ( c < 0 && aa != 'X' ) {
                              printf( "Error. Sequence letter unknown %c %s\n",
                                      aa, pl->pep );
                              exit( 1 );
			}

			/* Encode amino acid */
			/* For sparse encoding the amino acids are encoded as 1 p_max, and 19 p_min */
			/* For Blosum encoding the amino acids are encoded as the normalized blosum vector */
			/* For sparse encoding X's are coded as 20 p_min, for Blosum encoding X are coded as 20 0's */
			
			for ( j=0; j<20; j++ ) {

                              if ( ! p_bl ) { //sparse encoding used
                                      if ( j == c ) {
                                            //printf( "%5.3f ", XXXXXXXXX );
                                            printf( "%5.3f ", p_max ); //これまで1で埋めてた部分
                                      }
                                      else{
                                            //printf( "%5.3f ", XXXXXXXXX );
                                            printf( "%5.3f ", p_min ); //これまで0で埋めてた部分
                                      }
                                      
				}
                              else { //Blosum encoding used
                                      if ( c < 0 ){
                                            printf( "%5.3f ", 0.0 );
                                      }
                                      else{
                                            printf( "%5.3f ", blm[c][j] );
                                      }
                                      
				}
			}
		}

		printf( "%f\n", pl->score );

	}

	printf( "# Number of lines %i\n", nlines );

	exit( 0 );
}
