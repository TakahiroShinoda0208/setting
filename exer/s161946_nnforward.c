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

static double squash( a ) double a; { return (1.0 / (1.0 + exp( -a ))); }

int	p_ptarget;
int	p_synfile;
int	p_verbose;

PARAM   param[] = {
	"-pt", 	VINT	p_ptarget,	"Print target value", "1",
	"-v", VSWITCH	p_verbose,	"Verbose mode", "0",
	"-s", VSWITCH	p_synfile, 	"Use single synfile [default is synlist]", "0",
	0
};

/* GLOBAL variables */

float	*inp;

int	scan_input_fp( FILE *fp, int n )

{
	int 	i;
	char	ch;
	LINE	line;

	while ( ( ch = fgetc( fp ) ) && ungetc( ch, fp ) && ch == '#' ) { /* skip comment line */
                fgets( line, sizeof(line), fp );
                if ( p_verbose )
                        printf( "# %s", line );
        }

	for ( i = 0; i <= n; i++) /* Read input line with n elements */
        	if ( fscanf(fp, "%f", &inp[i]) != 1 ) 
			return( 0 );

	return( 1 );
}

void    setnet(float **x, SYNAPS * syn)

{
        int     l;
        for (l = 0; l < syn->nlay; l++) /* Initialize the extra input value for each layer (the bias) to 1 */
              x[l][syn->nnlay[l]] = 1.0; //s->nnlay[1] is the number of neurons in the hidden layer set bias input as 1
}

/* synaps data structure and routines to read synaps files are defined in utils/nnutils.{c,h}

The format of the synaps file is

TESTRUNID Ncycles 90 Train Pearson 0.931930 Err 0.008701 Test Pearson 0.933660 Err 0.008496
     180  LAYER:    1
       2  LAYER:    2
       1  LAYER:    3
      90   :ILEARN
     0.021232     0.295686     0.110727     0.566801    -0.045417
     0.108592     0.472350     0.125460     0.068841    -0.450012
    -0.057411    -0.122409    -0.171865    -0.651800     0.677134
    -0.051108     0.115779     0.125024    -0.643785     0.115704
     0.167212     0.556199     0.069067     0.378024     0.011984
     ...

The first line is a comment line, 
the second line (and all lines having the keywork "LAYER:" at position 2) defines the network structure
the line with the Keyword ":ILEARN" give the training cycle when the synaps file was saved. 
The rest of the file contains the synaps weights. The weights are ordered as

wgt[0][0][0] .. wgt[0][0][nin] wgt[0][1][0] .. wgt[0][1][nin] etc

The data structure is 

typedef struct synaps {
        struct synaps *next;
        LINE    nin, nout;
        FILENAME wgtfile;
        int     nlay, *nnlay; 
        float ***wgt;
        int     maxnper;
} SYNAPS;

As always this is a linked list. Each element in the list has the following data

s->nlay: number of layers in network. 0 is input layer, 1 is hidden etc.
s->nnlay: is a vector with the number of neurons in each layer. That is s->nnlay[0] is the number of
inputs to the network, s->nnlay[1] is the number of neurons in the hidden layer etc.
s->wgt: is a 3 dimentional matrix with the synaps weights. The first index refers to the FROM layer
the second index to the neuron in the TO layer, and the third index refers to the FROM neuron. That

s->wgt[0][0][1] refers to the synaps weight between the 0 and 1 layers between the 1 neuron 
in the 0 layers and the 0 neuron on the 1 layer.

s->maxnper: Contains the size of the largest layer in the network

The other variables are strings containing information about the file types etc. These are unimportant.

*/

main( int argc, char *argv[] )

{

	SYNAPS	*syn, *s;
	FILE    *fp;
	int	i, l, k;
	float	**x;
	float	a, *p;
	int	nlines, ns, nout, nin;
	int	fc, ff;
	int	synnlay;

	pparse( &argc, &argv, param, 2, "synapslist/synapsfile inputfile" );

	if ( p_synfile ) { /* Read single synaps file */
		syn = synaps_read_file( argv[1] );
		printf( "# Synaps file read from %s\n", argv[1] );
	}
	else /* Read synaps list */
		syn = synaps_read( argv[1] );

	/* Print out network structure for first network */
	for ( i=0; i<syn->nlay; i++ ) 
		printf( "# Nlayer %i NNlay %i\n", i, syn->nnlay[i] );

	/* Find maximum number of layers in ensemble */
	for ( synnlay=-99, s=syn; s; s=s->next ) 
		synnlay = ( s->nlay > synnlay ? s->nlay : synnlay );

	/* Open file for input */
        if ( ( fp = stream_input( argv[2], &fc, &ff ) ) == NULL ) { 
                printf( "Error Cannot open file %s for input. Exit\n", argv[2] );
                exit( 1 );
        }

	/* Numner of output neurons is given by the size of the last layer */
	nout = syn->nnlay[syn->nlay-1];

	/* Numner of input neurons is given by the size of the first layer */
	nin = syn->nnlay[0];


	/* Allocate matrix for feed forward */
	/* The matrix x contains the out from each layer. That is
		x[0][0], .. ,x[0][nin-1] contains the input, x[0][nin] = 1 (for the bias) 
		x[1][0], .., x[1][s->nnlay[1]-1] contains the output from the second layer and x[1][s->nnlay[1]] = 1

		Note

                x[1][0] = squash( sum_i x[0][i] * wgt[0][1][i] )

		where squash is the sigmoid function

	*/
	x = fmatrix( 0, synnlay-1, 0, syn->maxnper ); 

	inp = fvector( 0, nin+nout-1 ); /* Allocate vector for input */

	p = fvector( 0, nout-1 ); /* Allocate vector for output. Note network can have more than one output value */ 

	printf( "# Ninput %i. Noutput %i\n", nin, nout );

	nlines = 0;

	while ( scan_input_fp( fp, nin+nout-1 ) ) { /* Read an input line from file fp */

		for ( i=0; i<nout; i++ ) 
			p[i] = 0.0;

		for ( ns=0, s=syn; s; s=s->next ) { /* Run each Network/synaps */

			setnet( x, s ); /* Initialize bias input to 1 */

			/* Set x[0][i] to input values */
			for ( i=0; i<s->nnlay[0]; i++ ) 
				x[0][i] = inp[i]; 

			for (  l = 1; l < s->nlay; l++ ) { /* loop over all layers */

				for ( k = 0; k < s->nnlay[l]; k++ ) { /* loop over neurons in layer l */

					a = 0.0;
					/* Add inputs from layer l-1 */
					/* This is just the sum of the weight * the output value from the
					previous layer */
					for (i = 0; i < s->nnlay[l-1]+1; i++ ) { 
                                              //a += XXXXXXXX;
                                              a += (x[l-1][i]) * (s->wgt[l-1][k][i]); //l-1 layer from i cell to k cell.
                                              
					}

					x[l][k] = squash( a ); /* Squash input */

				}
			}

			for ( i=0; i<nout; i++ ) /* Save output values from syapse s */
                              //p[i] += XXXXX;
                              p[i] += x[s->nlay-1][i];
                        
                        
			ns++;

		}

		for ( i=0; i<nout; i++ ) {
			p[i] = ( ns> 0 ? p[i]/ns : 0.0 );

			printf( "%8.6f ", p[i] );

			if ( p_ptarget ) 
				printf( "Target= %8.6f ", inp[syn->nnlay[0]+i] );

		}

		printf( "\n" );

		nlines++;
	}

        stream_close( fp, fc, argv[2] );

	printf( "# Nlines %i processed thorugh %i Neural Networks\n", nlines, ns );

	exit( 1 );
}
