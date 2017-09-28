
/* 
	tree-dependent iteration   
    algorithm A+ when group-to-group, C when group-to-singleSeqence 
	                 OR
    algorithm A+
*/

#include "mltaln.h"

#define DEBUG 0
#define RECORD 0

static void Writeoption2( FILE *fp, int cycle, double cut )
{
	fprintf( fp, "%dth cycle\n", cycle );
    fprintf( fp, "marginal score to search : current score * (100-%d) / 100\n", (int)cut );
}

static void Writeoptions( FILE *fp )
{
	fprintf( fp, "Tree-dependent-iteration\n" );
    if( scoremtx == 1 )
        fprintf( fp, "Dayhoff( ... )\n" );
    else if( scoremtx == -1 )
        fprintf( fp, "DNA\n" );
    else if( scoremtx == 2 )
        fprintf( fp, "Miyata-Yasunaga\n" );
	else
		fprintf( fp, "JTT %dPAM\n", pamN );

	if( scoremtx == 0 )
    	fprintf( fp, "Gap Penalty = %+5.3f, %5.2f, %+5.3f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
	else
		fprintf( fp, "Gap Penalty = %+5.3f\n", (double)penalty/1000 );
		

    if( scmtd == 3 )
        fprintf( fp, "score of rnd or sco\n" );
    else if( scmtd == 4 )
        fprintf( fp, "score = sigma( score for a pair of homologous amino acids ) / ( number of amino acids pairs )\n" );
    else if( scmtd == 5 )
        fprintf( fp, "score : SP\n" );
    if( mix )
        fprintf( fp, "?\n" );
    else
    { 
        if( weight == 2 )
            fprintf( fp, "weighted rationale-1,  geta2 = %f\n", geta2 );
        else if( weight == 3 )
            fprintf( fp, "weighted like ClustalW," );
        else if( weight == 4 )
            fprintf( fp, "weighted rationale-2,  geta2 = %f\n", geta2 );
        else
            fprintf( fp, "unweighted\n" );
    }
    if( weight && utree )
        fprintf( fp, "using tree defined by the file hat2.\n" );
    if( weight && !utree )
        fprintf( fp, "using temporary tree.\n" );

	if( treemethod == 'n' )
		fprintf( fp, "Tree is calculated with Neighbor-Joining Method.\n" );
	else if( treemethod == 's' )
		fprintf( fp, "Tree is calculated with simplified UPG Method.\n" );
	else if( treemethod == 'x' )
		fprintf( fp, "Tree is calculated with simplified UPG Method and UPG Method.\n" );
	else if( treemethod == 'a' )
		fprintf( fp, "Tree is calculated with UPG Method.\n" );
		
	if( alg == 'C' ) 
		fprintf( fp, "Algorithm A+ / C\n" );
	else if( alg == 'S' ) 
		fprintf( fp, "Algorithm S \n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+ \n" );
	else if( alg == 'a' ) 
		fprintf( fp, "Algorithm A \n" );
	else 
		fprintf( fp, "Algorithm ? \n" );

	if( use_fft )
	{
		fprintf( fp, "Skip not diagonalregions by FFT\n" );
		if( scoremtx == -1 )
		{
			fprintf( fp, "Basis : 4 nucleotides\n" );
		}
		else
		{
			if( fftscore )
				fprintf( fp, "Basis : Polarity and Volume\n" );
			else
				fprintf( fp, "Basis : 20 amino acids\n" );
		}
		fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
		fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
	}
}


int TreeDependentIteration( int locnjob, char name[M][B], int nlen[M], 
                             char **aseq, char **bseq, int ***topol, double **len, 
                             int alloclen )
{
	int i, j, k, l, iterate, ii;
	int lin, ldf; 
	int clus1, clus2;
	int s1, s2;
	static char name1[M][B], name2[M][B];
	static Node *stopol;
	static double *effarr = NULL;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	static double **mtx = NULL;
	static int **node = NULL;
	static int *branchnode = NULL;
	static double **branchWeight = NULL;
	static char **mseq1, **mseq2;
	static float ***history;
	FILE *trap;
	double tscore, mscore;
	int identity;
	int converged;
	int oscillating;
#if 0
	char pair[njob][njob];
#else
	static char **pair;
#endif
#if DEBUG + RECORD
	double score_for_check0, score_for_check1;
	static double **effmtx = NULL;
	extern double score_calc0();
#endif
	extern double intergroup_score( char **, char **, double *, double *, int, int );
	static char *indication1, *indication2;

	Writeoptions( trap_g );
	fflush( trap_g );

	if( effarr == NULL ) /* locnjob == njob ni kagiru */
	{
		indication1 = AllocateCharVec( njob*3+50 );
		indication2 = AllocateCharVec( njob*3+50 );
		effarr = AllocateDoubleVec( locnjob );
		effarr1 = AllocateDoubleVec( locnjob );
		effarr2 = AllocateDoubleVec( locnjob );
		mseq1 = AllocateCharMtx( locnjob, 1 );
		mseq2 = AllocateCharMtx( locnjob, 1 );
		mtx = AllocateDoubleMtx( locnjob, locnjob );
		node = AllocateIntMtx( locnjob, locnjob );
		branchnode = AllocateIntVec( locnjob );
		branchWeight = AllocateDoubleMtx( locnjob, 2 );
		history = AllocateFloatCub( 100, locnjob, 2 );
		stopol = (Node *)calloc( locnjob * 2, sizeof( Node ) );
#if 0
#else
		pair = AllocateCharMtx( locnjob, locnjob );
#endif
	}
#if DEBUG + RECORD
	if( !effmtx ) effmtx = AllocateDoubleMtx( locnjob, locnjob );
	for( i=0; i<locnjob; i++ ) for( j=0; j<locnjob; j++ ) effmtx[i][j] = 1.0;
#endif

	for( i=0; i<locnjob; i++ ) strcpy( bseq[i], aseq[i] );

	writePre( locnjob, name, nlen, aseq, 0 );

	if( utree )
	{
		if( weight == 2 ) 
			countnode_int( locnjob, topol, node );
		else if( weight == 4 )
		{
			treeCnv( stopol, locnjob, topol, len, branchWeight );
			calcBranchWeight( branchWeight, locnjob, stopol, topol, len );
		}
	}

	converged = 0;
	if( cooling ) cut *= 2.0;
	for( iterate = 0; iterate<100; iterate++ ) 
	{
		if( cooling ) cut *= 0.5; /* ... */

		fprintf( trap_g, "\n" );
		Writeoption2( trap_g, iterate, cut );
		fprintf( trap_g, "\n" );

		if( utree == 0 )
		{
			if( devide )
			{
				static char *buff1 = NULL;
				static char *buff2 = NULL;
				if( !buff1 )
				{
					buff1 = AllocateCharVec( alloclen );
					buff2 = AllocateCharVec( alloclen );
				}	

				for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) 	
				{
					buff1[0] = buff2[0] = 0;
					strcat( buff1, res_g[i] ); strcat( buff2, res_g[j] );
					strcat( buff1,  bseq[i] ); strcat( buff2,  bseq[j] );
					strcat( buff1, seq_g[i] ); strcat( buff2, seq_g[j] );

					mtx[i][j] = (double)substitution_hosei( buff1, buff2 );
				}
			}
			else
			{
				for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) 	
					mtx[i][j] = (double)substitution_hosei( bseq[i], bseq[j] );
			}

			if     ( treemethod == 'n' )
				nj( locnjob, mtx, topol, len );
			else if( treemethod == 's' )
				spg( locnjob, mtx, topol, len );
			else if( treemethod == 'x' )
				supg( locnjob, mtx, topol, len );
			else if( treemethod == 'p' )
				upg2( locnjob, mtx, topol, len );

			if( weight == 2 )
				countnode_int( locnjob, topol, node );
			else if( weight == 4 )
			{
				treeCnv( stopol, locnjob, topol, len, branchWeight );
				calcBranchWeight( branchWeight, locnjob, stopol, topol, len );
			}
			trap = fopen( "hat2", "w" );
			if( !trap ) ErrorExit( "Cannot open hat2." );
			WriteHat2( trap, locnjob, name, mtx );
			fclose( trap );
		}

		if( iterate % 2 == 0 ) 
		{
			lin = 0; ldf = +1;
		}
		else
		{
			lin = locnjob - 2; ldf = -1;
		}	
		for( l=lin; l < locnjob-1 && l >= 0 ; l+=ldf )
		{
			for( k=0; k<2; k++ ) 
			{
				if( l == locnjob-2 ) k = 1;
#if 1
				fprintf( stderr, "STEP %03d-%03d-%d %s", iterate+1, l+1, k, use_fft?"\n":"" );
#else
				fprintf( stderr, "STEP %03d-%03d-%d %s", iterate+1, l+1, k, use_fft?"\n":"\n" );
#endif
				for( i=0; i<locnjob; i++ ) for( j=0; j<locnjob; j++ ) pair[i][j] = 0;

				OneClusterAndTheOther( locnjob, pair, &s1, &s2, topol, l, k );
#if 0
				fprintf( stderr, "STEP%d-%d\n", l, k );
				for( i=0; i<locnjob; i++ ) 
				{
					for( j=0; j<locnjob; j++ ) 
					{
						fprintf( stderr, "%#3d", pair[i][j] );
					}
					fprintf( stderr, "\n" );
				}
#endif
				if( !weight )
					for( i=0; i<locnjob; i++ ) effarr[i] = 1.0;
				else if( weight == 2 ) 
				{
					nodeFromABranch( locnjob, branchnode, node, topol, len, l, k );
					node_eff( locnjob, effarr, branchnode );
				}
				else if( weight == 4 )
				{
					weightFromABranch( locnjob, effarr, stopol, topol, l, k );
				}

				for( i=0; i<locnjob; i++ ) strcpy( aseq[i], bseq[i] );

				clus1 = conjuction( pair, s1, aseq, mseq1, effarr1, effarr, name, name1, indication1 );
				clus2 = conjuction( pair, s2, aseq, mseq2, effarr2, effarr, name, name2, indication2 );

				mscore = intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2 );

#if DEBUG
				if( !weight )
				{
					score_for_check0 = score_calc0( aseq, locnjob, effmtx, 0 );
				}
#endif

				commongappick( clus1, mseq1 );
				commongappick( clus2, mseq2 );
#if DEBUG
				fprintf( stderr, "mscore = %f\n", mscore );
#endif
	
#if DEBUG
				if( !devide )
				{
		       		fprintf( trap_g, "\nSTEP%d-%d-%d\n", iterate+1, l+1, k );
       		    	fprintf( trap_g, "group1 = %s\n", indication1 );
   	   		    	fprintf( trap_g, "group2 = %s\n", indication2 );
					fflush( trap_g );
				}
	
#endif
#if 0
				printf( "STEP %d-%d-%d\n", iterate, l, k );
				for( i=0; i<clus2; i++ ) printf( "%f ", effarr2[i] );
				printf( "\n" );
#endif
				if( use_fft )
				{
					Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
				}
				else
				{
					if( alg == 'A' )
					{
						A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					}
					else if( alg == 'a' ) 
					{
						Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					}

					else if( alg == 'C' ) 
					{
						if( clus1 == 1 && clus2 != 1 || clus1 != 1 && clus2 == 1 )
							translate_and_Calign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
						else
							A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					}
					else ErrorExit( "Sorry!" );
				}
							
					if( checkC )
					{
						extern double DSPscore();
						extern double SSPscore();
						static double cur;
						static double pre;
	
/*
						pre = SSPscore( locnjob, bseq );
						cur = SSPscore( locnjob, aseq );
*/
						pre = DSPscore( locnjob, bseq );
						cur = DSPscore( locnjob, aseq );
	
						fprintf( stderr, "Previous Sscore = %f\n", pre );
						fprintf( stderr, "Currnet  Sscore = %f\n\n", cur );
					}
					
				identity = !strcmp( aseq[s1], bseq[s1] );
				identity *= !strcmp( aseq[s2], bseq[s2] );

/* Bug?  : idnetitcal but score change when scoreing mtx != JTT  */
	
				if( identity )
				{
					tscore = mscore;
					if( !devide ) fprintf( trap_g, "tscore =  %f   identical.\n", tscore );
					fprintf( stderr, " identical." );
					converged++;
				}
				else
				{
					tscore = intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2 );

					if( tscore > mscore - cut/100.0*mscore ) 
					{
						writePre( locnjob, name, nlen, aseq, 0 );
						for( i=0; i<locnjob; i++ ) strcpy( bseq[i], aseq[i] );
	
#if DEBUG
						fprintf( trap_g, "tscore =  %f   mscore = %f  accepted.\n", tscore, mscore );
#endif
						fprintf( stderr, " accepted. " );
						converged = 0;
	
					}
					else 
					{
#if DEBUG
						fprintf( trap_g, "tscore =  %f   mscore = %f  rejected.\n", tscore, mscore );
#endif
						fprintf( stderr, " rejected. " );
						tscore = mscore;
						converged++;
					}
				}
				if( use_fft ) fprintf( stderr, "\n\n" );
				else          fprintf( stderr, "\r" );

#if DEBUG
				if( !weight )
				{
					score_for_check1 = score_calc0( aseq, locnjob, effmtx, 0 );
					fprintf( trap_g, "\n" );
					fprintf( trap_g, "all-pair-score : %f -> %f\n", score_for_check0, score_for_check1 );
					if     ( score_for_check0 < score_for_check1 ) fprintf( trap, "acceptable\n" );
					else if(score_for_check0 == score_for_check1 ) fprintf( trap, "equivalent\n" );
					else                                           fprintf( trap, "should be rejected\n" );
				}
#endif

#if RECORD
				score_for_check1 = score_calc0( aseq, locnjob, effmtx, 0 );
				fprintf( trap_g, "all-pair-score : %f\n", score_for_check1 );
#endif

				history[iterate][l][k] = (float)tscore;
	
				if( converged >= locnjob * 2 )
				{
					fprintf( trap_g, "Converged.\n\n" );
					fprintf( stderr, "\nConverged.\n\n" );
					return( 0 );
				}
				if( iterate >= 1 )
				{
	/*   oscillation check    */
					oscillating = 0;
					for( ii=iterate-2; ii>=0; ii-=2 ) 
						oscillating += ( tscore == history[ii][l][k] );
					if( ( oscillating && !cooling ) || ( oscillating && cut < 0.001 && cooling ) )
					{
						fprintf( trap_g, "Oscillating.\n" );
						fprintf( stderr, "\nOscillating.\n\n" );
#if 1 /* hujuubun */
						return( -1 );
#endif
					}
				}      /* if( iterate ) */
			}          /* for( k ) */
		}              /* for( l ) */
	}                  /* for( iterate ) */
}                  	   /* int Tree... */
