/* Ole Lund Jan 2002 lund@cbs.dtu.dk */
/* M. Nielsen Feb 2002 mniel@cbs.dtu.dk */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

int	p_l;
int	p_niter;
int	p_seed;
FILENAME	p_matfile;
FILENAME	p_corefile;
WORD	p_alphabet;
int	p_printcore;
float	p_tstart;
float	p_tfinal;
int	p_nt;
float   p_wlowcount;
int	p_seqw;
FILENAME        p_blrealmat;
int	p_verbose;
float	p_nocountaalo;
int	p_p1allowed;
FILENAME p_p1letters;
float	p_factor;

PARAM   param[] = {
	"-l", VINT	p_l, "Motif length", "9",
	"-i", VINT	p_niter, "Number of iterations per train example", "50",
	"-s", VINT	p_seed, "Seed [-1] Default, [0] Time [>0] Specific seed", "-1",
	"-m", VFNAME	p_matfile, "File for scoring matrix", "--",
	"-alp", VWORD	p_alphabet, "Amino acid alphabet", "ARNDCQEGHILKMFPSTWYV",
	"-pc", VSWITCH p_printcore, "Print core of alignment", "0",
	"-c", VFNAME	p_corefile, "File for core", "core",
	"-ts", VFLOAT	p_tstart, "Start temperature", "3.0",
	"-te", VFLOAT	p_tfinal, "Final temperature", "0.0001",
	"-nt", VINT	p_nt, "Number of temperature steps", "10",
	"-wlc", VFLOAT  p_wlowcount, "Weigth on low count", "50.0",
	"-sw", VSWITCH	p_seqw, "Include heuristic sequence weighting", "0",
        "-blf", VFNAME p_blrealmat, "Real value blosum matrix filename", "$ALGOHOME/data/blosum62.freq_rownorm",
	"-v", VSWITCH	p_verbose, "Verbose mode", "0",
	"-noc", VFLOAT	p_nocountaalo, "Log-odds score asignment to no count elements", "0.0",
	"-p1a", VSWITCH p_p1allowed, "Allow only Hydrophobic AA at P1", "0",
	"-p1l", VFNAME p_p1letters, "Allowed letters at P1 [Only with -p1a option]", "ILVMFYW",
	"-factor", VFLOAT p_factor, "Scale factor for matrix 2/log2", "2.8854",
	0
};

float	*set_bgfreq_aa( int alen )
//background frequency
{
	float	*aabgfrec;

	aabgfrec = fvector(0, alen-1);

	aabgfrec[0]  = 0.074;
	aabgfrec[1]  = 0.052;
	aabgfrec[2]  = 0.045;
	aabgfrec[3]  = 0.054;
	aabgfrec[4]  = 0.025;
	aabgfrec[5]  = 0.034;
	aabgfrec[6]  = 0.054;
	aabgfrec[7]  = 0.074;
	aabgfrec[8]  = 0.026;
	aabgfrec[9]  = 0.068;
	aabgfrec[10] = 0.099;
	aabgfrec[11] = 0.058;
	aabgfrec[12] = 0.025;
	aabgfrec[13] = 0.047;
	aabgfrec[14] = 0.039;
	aabgfrec[15] = 0.057;
	aabgfrec[16] = 0.051;
	aabgfrec[17] = 0.013;
	aabgfrec[18] = 0.032;
	aabgfrec[19] = 0.073;

	return( aabgfrec );
}

typedef struct seqlist {
        struct seqlist *next;
        char    *pep;
        int     len;
        int     *i;
        float   w;
        int     clen;
        int     offset;
        float   aff;
} SEQLIST;

SEQLIST *seqlist_alloc()

{
        SEQLIST *n;

        if ( ( n = ( SEQLIST * ) malloc ( sizeof( SEQLIST ))) != NULL ) {
		n->next = NULL;
                n->pep = NULL;
                n->len = 0;
                n->i = NULL;
                n->w = 1.0;
                n->clen = 0;
                n->offset = -99;
		n->aff = -99;
	}

	return( n );
}

SEQLIST *seqlist_read( char *filename ) 

{

	SEQLIST         *first, *last, *new;
	LINELIST	*linelist, *ln;
	float           aff;
        WORD            pep;
	int		n;
	
	linelist = linelist_read( filename );

	first = NULL;
	n = 0;

	for ( ln=linelist; ln; ln=ln->next ) {

		if ( strlen( ln->line ) < 1 )
                        continue;

                aff = -99.9;

                if ( sscanf( ln->line, "%s %f", pep, &aff ) < 1 ) {
                        printf( "Error. Wrong line format. %s\n", ln->line );
                        exit( 1 );
                }

                if ( ( new = seqlist_alloc()) == NULL ) {
                        printf( "Error. Cannot allocate seqlist. Exit\n" );
                        exit( 1 );
                }

                new->len = strlen(pep);
                new->aff = aff;

                new->pep = cvector( 0, new->len );
                strcpy( new->pep, pep );

               	if ( first == NULL )
                        first = new;
                else
                        last->next = new;

                last = new;

		n++;

	}

	linelist_free( linelist );

	printf( "# Read %i elements on seqlist %s\n", n, filename );

        return( first );

}

void	seqlist_free( SEQLIST	*seqlist )

{
	SEQLIST	*s, *n;

	for ( s=seqlist; s; s=n ) {

		n = s->next;

		if ( s->pep )
			cvector_free( s->pep, 0, s->len );

		free( s );
	}
}

SEQLIST	**seqlist_table( int n )
                
{               
        SEQLIST **table; 
        int     i;
                
        if ( ( table = ( SEQLIST ** ) malloc(( unsigned )
                n * sizeof( SEQLIST * ) ) ) == NULL ) {
                printf( "Error. Cannot allocate SEQLIST TABLE %i\n", n );
                exit( 1 );
        }       
                
        for ( i=0;i<n;i++ )
                table[i] = NULL;
        
        return( table );
}       

float   sequence_weighting( SEQLIST *peplist, char *alphabet, int alen )

/* This code is taken from the pep2mat program */

{
        SEQLIST *pl;
        int     *r, **s;
        int     len;
        float   nseq;
        int     i,j, iv, o;
        float   w;

        len = peplist->clen;

        r = ivector( 0, len-1 );
        s = imatrix( 0, len-1, 0, alen-1 );

        /* s is a matrix containing the number occurrence of each amino acids at each position */
        for ( pl=peplist; pl; pl=pl->next ) {

		/* o is the offset value for the given peptide */
		/* That is the position where the core starts */
		o = pl->offset;

                for ( i=0; i<len; i++ ) {

			/* Find the encoding the the amino acids at posiiton i+o */
                        iv = pl->i[i+o];

                        s[i][iv]++; /* This is the sam as s[i][iv] = s[i][iv] + 1; */

                }

        }

        /* r is a vector containing the number of different amino acids at each position */
        for ( i=0; i<len; i++ ) {

                for ( j=0; j<alen; j++ ) {

                        r[i] += ( s[i][j] > 0 );
                }
        }

        for ( pl=peplist; pl; pl=pl->next ) {

                pl->w = 0.0;
		o = pl->offset;

                for ( i=0; i<len; i++ ) {

			/* Find the encoding the the amino acids at posiiton i+o */
                        iv = pl->i[i+o];

                        w = 1.0/( r[i] * s[i][iv] );
                        pl->w += w;
                }
        }

	/* The effective number of sequence is the average of r over the peptide positions */
        nseq = 0;
        for ( i=0; i<len; i++ ) 
                nseq += r[i];

        nseq /= len;

        return( nseq );
}

float   cal_kldsum( SEQLIST *peplist, float **aafrec, float **aalo,  
        float **aacounts, float *aabgfrec, float **aaraw, float **blf,
        float wlowcount, int seqw, int alen )

{

	int     i, j, k, len, ix;
        int     o;
        SEQLIST *pep;
        float   kldsum, odds;
        float   neff, sum, g;
	int	nseq;

	len = peplist->clen;

	/* Initialize count matrix */
	/* aacounts contain the sequence weigthed amino acid counts */
	/* aaraw contain the raw amino acid counts */
        for (j=0; j<len; j++) {
                for (k=0;k<alen;k++) {
                        aacounts[j][k]=0;
                        aaraw[j][k]=0;
                }
        }

        /* Estimate sequences weigths */
        if ( seqw )
                neff = sequence_weighting( peplist, p_alphabet, alen ); 
	else {
                for ( neff=0, pep=peplist; pep; pep=pep->next,neff+=1.0 )
                        pep->w = 1.0;
        }

	if ( p_verbose )
		printf( "# Number of sequences nseq %f\n", neff );

        /* Now find the amino acid counts at each peptide position */

	sum = 0;
	nseq = 0;
        for ( pep=peplist; pep; pep=pep->next ) {

		o = pep->offset;

                for ( i=0; i<len; i++ ) {

			/* Find the amino acids type at position i+o */
                        ix = pep->i[i+o];

                        aacounts[i][ix] += pep->w;
			aaraw[i][ix] += 1.0;

                }

		sum += pep->w;
		nseq++;

        }

	/* Normalize the amino acids counts to frequencies */
        for ( i=0; i<len; i++ ) {

                /* aacounts[i][j] /= ( sum > 0 ? sum : 1.0 ); is a compact syntax for dividing with 1.0
                if sum is <= or dividing with sum otherwise */

                for ( j=0; j<alen; j++ ) {
                        aacounts[i][j] /= ( sum > 0 ? sum : 1.0 );
			aaraw[i][j] /= (nseq > 0 ? nseq : 1.0 );
		}
        }

       	/* Calculate pseudo frequency g, and combined frequence (a*f+b*g)/(a+b), where a is neff-1 and b is
        (the command line parameter) weight on prior */

        for ( i=0; i<len; i++ ) {

		sum = 0.0;
                for ( j=0; j<alen; j++ ) {

                        /* g_j = sum_k f_k * q(j|k) */
                        /* f_k =  aacounts[i][k]; */
                        /* q(j|k) = blf[k][j] */

                        g = 0.0;
                        for ( k=0; k<alen; k++ )
                                g += aacounts[i][k] * blf[k][j];

                        aafrec[i][j] = ((neff-1) * aacounts[i][j] + p_wlowcount * g)/( (neff-1) + p_wlowcount ); 

			sum += aafrec[i][j];

                }

		for ( j=0; j<alen; j++ )
			aafrec[i][j] /= ( sum > 0.0 ? sum : 1.0 );

        }

#define MINODDS 0.0001

	/* Calculate odds and log odds ratios */
	for (j=0; j<len; j++)
        for (k=0; k<alen; k++) {

                if ( aafrec[j][k] == 0.0 ) {
                        aalo[j][k] = 0.0;
                }
                else {
                      odds = (aafrec[j][k])/aabgfrec[k]; //PSSM初期値
                        aalo[j][k] = ( odds > MINODDS ? log(odds) : log(MINODDS) ); 
                        aalo[j][k] *= p_factor;
                }
        }

        /* Calculate Kullbach Leibler distance sum */
	/* kldsum = sum_position sum_amino_acids number_aa_pos * aalo_aa_pos
         　総ポジション数* 総アミノ酸数* position pasでaa が観測された数 * 観測されたアミノ酸の重み付けlog score
           where number_aa_pos is the number of time amino acids aa is observed at position pos
	and aalo_aa_pos is the log_odds score of amino acid aa at position pos */
        


        kldsum =0.0;
        for (j=0; j<len; j++)
        for (k=0; k<alen; k++)
                kldsum += aaraw[j][k]*aalo[j][k];

        kldsum /= 2.0;

	return( kldsum );
}

main( int argc, char *argv[] )

{
    	FILE 		*fp;
	int		fc;
	SEQLIST		*peplist, *pl, **pep_table, *opl, *npl;
	int		*startposbest;
	float		*aabgfrec;
	float		de;
	int		nseq;
	int		i,j,k,acc,is, it, nw;
	int		iter, ns, niter;
	float 		**aacounts;
	float 		**aafrec;
	float 		**aalo;
	float           **aaraw;
	float		**blf;
	char		*bl_alphabet;
	int		o, old;
	float		kldsum, kldsumbest;
	float		t, dt;
	float		eo, en;
	float		tmp;
	int		rem;
	int		found;
	int		ix, alen;
	WORD		core;

	/* Parse command line options */
	pparse( &argc, &argv, param, 1, "pep_file" );

	if ( p_seed >= 0 )
		setseed( p_seed );

	/* Read peptide entries */
	peplist = seqlist_read( argv[1] );

	/* Count entries and remove too short peptides and peptides without P1 amino acids */
        // P1は プログラムの開始上必要なので、十分な長さがあるamino acidでも、P1がない場合は除去
	for ( nseq=0, opl=NULL, pl=peplist; pl; pl=npl ) {

		npl = pl->next;

		rem = 0;

		if ( pl->len < p_l ) {
			printf( "# Sequence in shorter than pattern length %s %i. Remove\n",
				pl->pep, p_l );
			rem = 1;
		}
		else if ( p_p1allowed ) {
			/* check is peptide contains at least one allowed amino acids to fit the P1 position */
                      // P1があるかどうか判定
			found = 0;
			for ( i=0; i<=pl->len-p_l && !found; i++ )
				if ( strpos( p_p1letters, pl->pep[i] ) >= 0 )
					found = 1;

			if ( ! found ) {
				printf( "# Sequence does not contain allowed P1 aa %s %s. Remove\n",
					pl->pep,  p_p1letters );
				rem = 1;
			}
		}

		if ( rem ) {
			if ( pl==peplist ) {
                                peplist = npl;
                                pl->next = NULL;
                                seqlist_free( pl );
                        }
                        else {
                                opl->next = npl;
                                pl->next = NULL;
                                seqlist_free( pl );
                        }
                }
                else {
                        nseq++;
                        opl = pl;
                }

	}

	printf( "# Number of sequences nseq %i\n", nseq );

        if ( peplist == NULL ) {
                printf( "No PEPLIST left in list. Stop\n" );
                exit( 0 );
        }

	/* Make table of peptides. This allows for fast access to each peptide entry */

	pep_table = seqlist_table( nseq );

	for ( i=0, pl=peplist; pl; pl=pl->next, i++ ) {
		pl->clen = p_l; /* assign the core length to clen */
		pl->i = ivector( 0, pl->len-1 ); /* allocate the vector i to contain the amino acids translated to integers */
		for ( j=0; j<pl->len; j++ ) {

			ix = strpos( p_alphabet, pl->pep[j] );

			if ( ix < 0 ) { 

				printf( "Error. Amino acid %c not part of alphabet %s\n",
					pl->pep[j], p_alphabet );
				exit( 1 );
			}

                        //リストi番目のペプチドj 番目のアミノ酸はp_alphabetでは
                        //ix番目である。
			pl->i[j] = ix; /* save the amino acid type */
                        
		}

		pep_table[i] = pl;
	}

	alen = strlen( p_alphabet );

	/* Background frequencies */
	aabgfrec = set_bgfreq_aa(alen);

	printf( "# Back ground frequences\n" );
	for ( i=0; i<alen; i++ )
		printf( "# %i %c %8.5f\n", i, p_alphabet[i], aabgfrec[i] );

	/* Assign random start positions */
	/* Allocate vector for optimal start positions */
	startposbest = ivector(0, nseq-1);

	for ( i=0,pl=peplist; pl; pl=pl->next,i++ ) {
		if ( p_p1allowed ) {
                      //乱数でランダムにスタート位置を割り振る
                        while ( ( o = (int)( drand48() * (pl->len-p_l+1) ) ) >=0 &&
                                ( strpos( p_p1letters, pl->pep[o] ) < 0 ) ) {
                        }
			pl->offset= o;
		}
		else
			pl->offset = (int)( drand48() * ( pl->len - p_l ) );

		startposbest[i] = pl->offset;
	}
        //BLOSUM matrix(substitution matrixの読み込み)
	blf = read_realblosum( p_blrealmat, &bl_alphabet );
		
       	if ( strcmp( bl_alphabet, p_alphabet ) ) {
               	printf( "Error. Bl alphabet %s not equal to p_alphabet %s\n",
                       	bl_alphabet, p_alphabet );
               	exit( 1 );
       	}

	printf( "# Motif length %i. Alen %i\n", p_l, alen );

	dt = ( p_tstart - p_tfinal)/(p_nt > 0 ? p_nt : 1 );

	printf( "# Number of temperature steps %i dt %f\n", p_nt, dt );

	/* Allocate matrix aa counts and frequences */
        aacounts = fmatrix(0, p_l-1, 0, alen-1);
	aaraw = fmatrix(0, p_l-1, 0, alen-1);
        aafrec = fmatrix(0, p_l-1, 0, alen-1);
        aalo = fmatrix(0, p_l-1, 0, alen-1);


        //aaloの初期値設定

	kldsum = cal_kldsum( peplist, aafrec, aalo, 
		aacounts, aabgfrec,
		aaraw, blf, p_wlowcount, p_seqw, alen );

	printf("# nt Ini T Ini Iteration: Ini Kullbach Leibler distance sum %f\n",
		kldsum);

	kldsumbest = kldsum;

	niter = nseq * p_niter; //effective number of seq * set number of iterator

	for ( it=0, t=p_tstart; it <= p_nt; it++, t-=dt ) {

	for ( iter=0;iter<niter;iter++) {

		if ( p_verbose )
			printf("# nt %i T %f Iteration: %d Kullbach Leibler distance sum %f\n", 
				it, t, iter+1, kldsum);

		/* 
		*	Select a random sequence
			Note the is not use is selecting peptides of length p_l since these are
			by default corectly aligned 
		*/
		while ( ( is = (int)(drand48()*nseq)) &&
                        ( (pep_table[is])->len == p_l ));

		pl = pep_table[is];  //plはランダムに選ばれたseq一つ
                
		/* 
		*	Select random peptide offset. If p_p1allowed is define, only select offsets
		*	where the P1 amino acids is one the allowed amino acids
		*/

		if ( p_p1allowed ) {
			while ( ( o = (int)( drand48() * (pl->len-p_l+1) ) ) >=0 && 
				( strpos( p_p1letters, pl->pep[o] ) < 0 ) ) {
			}

		}
		else {
			o = (int)( drand48() * (pl->len-p_l+1) );
		}
		
		/* Safe old offset */
		old = pl->offset;

		/* de = E_new - E_old, and E is the score of the peptide core to
		the PSSM (aalo), and k is a sum over peptide positions. Remember the 
                use of pl->i[k] to get the amino acids type at posiiton k in the peptide */

		de = 0;
                //p_l = motif length
		for (k=0; k<p_l;k++) {

                      //de += XXX;
                      de +=  aalo[k][pl->i[k]];
                      //aaloはPSSMで、Enewを求めるときはEoldから計算したPSSMを使う。
                      //これを繰り返すことで、練習問題でやったように重み付けがモデルにfitするよう変化していく。
                      //aaloの初期値は上で設定している。↑↑↑
                      
		}

		if ( de >= 0 ) {
			acc = 1;
		}
		else {
			/* Note we are trying to maximize the fitness function */
			if ( exp( de/t  ) > drand48() ) 
				acc = 1;
			else
				acc = 0;
		}

                //for文回ってないから、ここで計算できる Eって一つのpeptideだけなんじゃないの?
		if ( acc ) {

			pl->offset = o;

			kldsum = cal_kldsum( peplist, aafrec, aalo,
                		aacounts, aabgfrec,
                		aaraw, blf, p_wlowcount, p_seqw, alen );

			if ( kldsum > kldsumbest ) {
				kldsumbest = kldsum;
				startposbest[is] = (pep_table[is])->offset;
			}

			if ( p_verbose )
				printf( "# T %f Itr %i kldsum %f klsdum best %f\n", 
					t, iter, kldsum, kldsumbest );
		}


	} /* end iteration loop */

	} /* end temperature loop */
        
	kldsum = cal_kldsum( peplist, aafrec, aalo,
        		aacounts, aabgfrec,
        		aaraw, blf, p_wlowcount, p_seqw, alen );

	printf("# Iteration: Final Kullbach Leibler distance sum %f %f\n", kldsum, kldsumbest);

	/* 
	*	Print logodds  matrix
	*/

	if ( ( fp = stream_output( p_matfile, &fc, 0)) == NULL) {
		printf( " Error. Can't open %s\n", p_matfile);
		exit(1);
	}

	fprintf(fp, "\nLast position-specific scoring matrix computed       \n");

	fprintf(fp,"      " );
	for (k=0;k<alen; k++)
		fprintf(fp, "%7c ", p_alphabet[k] );
	fprintf(fp, "\n");

	for (j=0;j<p_l;j++) {
		fprintf(fp, "%3d A ",j+1);
		for (k=0;k<alen;k++) 
			fprintf(fp, "%7.3f ", aalo[j][k]);
		fprintf(fp, "\n");
	}

	stream_close( fp, fc, p_matfile );

	printf( "# Log-odds matrix written to file %s\n", p_matfile );

	/* 
	*	Print alignment
	*/

	/* Print core */

	if ( p_printcore ) {

		if ( ( fp = stream_output( p_corefile, &fc, 0 ) ) == NULL) {
			printf( "Error. Can't open %s\n", p_corefile);
			exit(1);
		}

		for ( pl=peplist; pl; pl=pl->next ) {
			strncpy( core, pl->pep + pl->offset, p_l );
			core[p_l] = 0;
			fprintf(fp, "%s %f\n", core, pl->w);
		}

		stream_close( fp, fc, p_corefile );
	}
		
	exit( 0 );
}
