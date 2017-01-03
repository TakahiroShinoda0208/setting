/* M.Nielsen July 2008, mniel@cbs.dtu.dk */

/* 

Copyright (C) 2008 Danish Technical University

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


/* First you inlude the standard libraries you might need */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

FILENAME     p_blmat;

PARAM   param[] =
{
	"-blm", VFNAME	p_blmat, "File with BLOSUM substitution scores", "$ALGOHOME/data/BLOSUM50",
	0
};

/* Note, $HOME is a variable that point to your home directory and is replaced by
/home/people/stud0XX when the code is executed. */

/* Now follows subroutines */

/* Routine to read Blosum matrix file */

float	**read_blosummat( char *filename, char **alphabet )

{
	LINELIST	*linelist, *ln;
	float		**m;
	int		nc, l, i;
	char		**wvec;

	linelist = linelist_read( filename );

	if ( ! linelist ) {
		printf( "Error. Cannot read from file %s\n", filename );
		exit( 1 );
	}

        (*alphabet) = cvector( 0, 20 );
        m = fmatrix( 0, 19, 0 ,19 );
	l = 0;
          
#define MFORMAT "%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"

	for ( ln=linelist; ln; ln=ln->next ) {

		if ( strlen( ln->line ) <= 1 ) continue;
               
                if ( strncmp( ln->line, "#", 1 ) == 0 ) continue;
               
                if ( strncmp( ln->line, "   A", 4 ) == 0 ) {

			wvec = split( ln->line, &nc );
               
                        for ( i=0;i<20;i++ )
                                (*alphabet)[i] = wvec[i][0];
               
                        (*alphabet)[20] = 0;

                        continue;
               
                }
	
		if ( l < 20 && sscanf( ln->line, MFORMAT,
                        &m[l][0], &m[l][1], &m[l][2], &m[l][3], &m[l][4],
                        &m[l][5], &m[l][6], &m[l][7], &m[l][8], &m[l][9],
                        &m[l][10], &m[l][11], &m[l][12], &m[l][13], &m[l][14],
                        &m[l][15], &m[l][16], &m[l][17], &m[l][18], &m[l][19] ) != 20 ) {
               
                        printf( "Wrong line format %s", ln->line );
                        exit( 1 );
                }
               
                l++;
        }      

	linelist_free( linelist );
               
        printf( "# Read_realblosum done. Alphabet %s\n", *alphabet );

	return( m );
}

main(int argc, char *argv[])

{

/* First declair your variables */

	FSALIST		*fsa1, *fsa2;
	float		**m;
	float		**blmat;
	char		*alphabet;
	int		l1, l2, i, j;
	int		i1, i2;
	int		ix1, ix2;
	float		score;

/* Parse the command line options */

	pparse(&argc, &argv, param, 2, "fsa1 fsa2");

/* Start coding */

/* Read first fasta file */ 

	if ( ( fsa1 = fsalist_read( argv[1] ) ) == NULL ) {
		printf("Cannot read fasta file %s\n", argv[1] );
		exit(1);
	}

/* Read second fasta file */

	if ( ( fsa2 = fsalist_read( argv[2] ) ) == NULL ) {
		printf("Cannot read fasta file %s\n", argv[1] );
		exit(1);
	}

/* Read Blosum matrix */

	blmat = read_blosummat( p_blmat, &alphabet );

	l1 = fsa1->len;
	l2 = fsa2->len;

/* Allocate matrix to store scoring matrix */

	m = fmatrix( 0, l1-1, 0, l2-1 );

/* Build matrix */

/* Here you must code the missing part. (The part with XXXX) */

/* An important function to use could be 

strpos( string, character )

This routine returns the position of the character in the string. 
Ex strpos( "ARNDCQEGHILKMFPSTWYV", "N" ) = 2 

*/

	for ( i=0; i<l1; i++ ) {
              
              int h1;
              h1 = strpos(alphabet, fsa1->seq[i]);
              
              for ( j=0; j<l2; j++ ) {
                    
                    int h2;
                    h2 = strpos(alphabet, fsa1->seq[j]);
                    m[i][j] = blmat[h1][h2];
                    
              }
              
 	}

/* Here your coding ends */

	for ( j=0;j<l2;j++ ) 
		printf( "\t%i%c", j+1, fsa2->seq[j] );
	printf( "\n" );

	for ( i=0;i<l1;i++ ) {
		printf( "%i%c", i+1, fsa1->seq[i] );
		for ( j=0; j<l2; j++ ) 
			printf( "\t%3.0f", m[i][j] );
		printf( "\n" );
	}

	exit(0);
}
