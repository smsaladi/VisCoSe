#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0

void arguments( int argc, char *argv[] )
{
    int c;
 
	calledByXced = 0;
	devide = 0;
    weight = 3;
    utree = 1;
	tbutree = 0;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'C';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'x';

    while( --argc > 0 && (*++argv)[0] == '-' )
        while ( c = *++argv[0] )
            switch( c )
            {
				case 'Q':
					calledByXced = 1;
					break;
				case 's':
					treemethod = 's';
					break;
				case 'x':
					treemethod = 'x';
					break;
				case 'p':
					treemethod = 'p';
					break;
				case 'a':
					alg = 'a';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'C':
					alg = 'C';
					break;
				case 'd':
					disp = 1;
					break;
				case 'o':
					outgap = 0;
					break;
				case 'j':
					tbutree = 1;
					break;
				case 'f':
					tbrweight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'u':
					tbitr = 1;
					tbrweight = 0;
					break;
				case 'm':
					tbitr = 1;
					tbrweight = 3;
					break;
				case 'w':
					tbweight = 2;
					break;
				case 'Z':
					checkC = 1;
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
	if( alg == 'C' && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : C, o\n" );
		exit( 1 );
	}
	readOtherOptions( &ppid, &fftThreshold, &fftWinSize );
}


void treebase( char name[M][B], int nlen[M], char **seq, char **aseq, char **mseq1, char **mseq2, double **mtx, int ***topol, double **len, double **eff, int alloclen )
{
	int i, j, l;
	int clus1, clus2;
	int s1, s2, r1, r2;
	float pscore;
	char indication1[njob*3+100], indication2[njob*3+100];
	FILE *trap;
	static char name1[M][B], name2[M][B];
	static double **partialmtx = NULL;
	static int ***partialtopol = NULL;
	static double **partiallen = NULL;
	static double **partialeff = NULL;
	static double *effarr = NULL;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
#if 1
	char pair[njob][njob];
#else
	char **pair;
#endif
	if( partialtopol == NULL ) 
	{
		partialmtx = AllocateDoubleMtx( njob, njob );
		partialtopol = AllocateIntCub( njob, 2, njob );
		partialeff = AllocateDoubleMtx( njob, njob );
		partiallen = AllocateDoubleMtx( njob, 2 );
		effarr = AllocateDoubleVec( njob );
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
#if 1
#else
		pair = AllocateCharMtx( njob, njob );
#endif
	}

	if( checkC )
		for( i=0; i<njob; i++ ) fprintf( stderr, "eff in tb-%d %f\n", i, eff[i][i] );

	for( i=0; i<njob; i++ ) effarr[i] = eff[i][i];
	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );

	if( checkC )
		for( i=0; i<njob; i++ ) fprintf( stderr, "effarr for aseq-%d %f\n", i, effarr[i] );

	writePre( njob, name, nlen, aseq, 0 );

	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) pair[i][j] = 0;
	for( i=0; i<njob; i++ ) pair[i][i] = 1;
	for( l=0; l<njob-1; l++ )
	{
		s1 = topol[l][0][0];
		for( i=0; (r1=topol[l][0][i])>-1; i++ ) 
			if( pair[s1][r1] != 1 ) exit( 1 );
		s2 = topol[l][1][0];
		for( i=0; (r2=topol[l][1][i])>-1; i++ ) 
			if( pair[s2][r2] != 1 ) exit( 1 );

		clus1 = conjuction( pair, s1, aseq, mseq1, effarr1, effarr, name, name1, indication1 );
		clus2 = conjuction( pair, s2, aseq, mseq2, effarr2, effarr, name, name2, indication2 );
		fprintf( trap_g, "\nSTEP-%d\n", l );
		fprintf( trap_g, "group1 = %s\n", indication1 );
		fprintf( trap_g, "group2 = %s\n", indication2 );

		if( checkC )
			for( i=0; i<clus1; i++ ) fprintf( stderr, "STEP%d-eff for mseq1-%d %f\n", l+1, i, effarr1[i] );
/*
		fprintf( stderr, "before align all\n" );
		display( aseq, njob );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

		if( alg == 'a' ) 
			pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
		else if( alg == 'A' )
			pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
		else if( alg == 'C' ) 
		{
			if( outgap && ( clus1 == 1 && clus2 != 1 || clus1 != 1 && clus2 == 1 ) )
				pscore = translate_and_Calign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
			else
				pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
		}
		else ErrorExit( "ERROR IN SOURCE FILE" );
/*
		fprintf( stderr, "after align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "after align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/
		for( i=0; (r2=topol[l][1][i])>-1; i++ ) 
		{
			pair[s1][r2] = 1;
			pair[s2][r2] = 0;
		}

		writePre( njob, name, nlen, aseq, 0 );

		if( clus1+clus2 > 2 )
		{
			mdfymtx( pair, s1, partialmtx, mtx );

			if( tbitr )
			{	
				mseqcat( mseq1, mseq2, partialeff, effarr1, effarr2, name1, name2, clus1, clus2 );
				fprintf( stderr, "treebase before itr %s /%s     %f\n", indication1, indication2, pscore );
				display( mseq1, clus1+clus2 );
				fprintf( stderr, "\n" );

				utree = 1; weight = tbrweight;
				iteration( clus1+clus2, name1, nlen, mseq1, seq, partialtopol, partiallen, partialeff );

				writePre( njob, name, nlen, aseq, 0 );
			}

			if( tbweight )
			{	
				mseqcat( mseq1, mseq2, partialeff, effarr1, effarr2, name1, name2, clus1, clus2 );
				fprintf( stderr, "treebase before itr %s /%s     %f\n", indication1, indication2, pscore );
				display( mseq1, clus1+clus2 );
				fprintf( stderr, "\n" );
	
				weight = 2;

				if( tbutree == 0 )
				{	
					utree = 0; 
					treeconstruction( mseq1, clus1+clus2, partialtopol, partiallen, partialeff );
					utree = 1;
				}
				else 
				{
					spg( clus1+clus2, partialmtx, partialtopol, partiallen );
					counteff( clus1+clus2, partialtopol, partiallen, partialeff );
				}
#if DEBUG

#endif

				iteration( clus1+clus2, name1, nlen, mseq1, seq, partialtopol, partiallen, partialeff );
				weight = 0;

				writePre( njob, name, nlen, aseq, 0 );
			}
		}

		fprintf( stderr, "treebase %d %s /%s\n", l+1, indication1, indication2 );
		if( disp ) display( aseq, njob );
		fprintf( stderr, "\n" );

	}
}

static void WriteOptions( FILE *fp )
{
	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "reversely weighted\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}
	if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
	else if( scoremtx ==  1 ) fprintf( fp, "Dayhoff( machigai ga aru )\n" );
	else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	else if( scoremtx == -1 ) fprintf( fp, "DNA\n" );

    if( scoremtx == 0 )
        fprintf( fp, "Gap Penalty = %+d, %+d\n", penalty, offset );
    else
        fprintf( fp, "Gap Penalty = %+d\n", penalty );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Apgorithm A+\n" );
	else if( alg == 'S' ) 
		fprintf( fp, "Apgorithm S\n" );
	else if( alg == 'C' ) 
		fprintf( fp, "Apgorithm A+/C\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'x' )
		fprintf( fp, "simplified UPGMA and UPGMA method.\n" );
	else if( treemethod == 's' )
		fprintf( fp, "simplified UPGMA method.\n" );
	else if( treemethod == 'p' )
		fprintf( fp, "UPGMA method.\n" );
	else
		fprintf( fp, "Unknown tree.\n" );
}
	 

int main1( int argc, char *argv[] )
{
	int  nlen[M];	
	static char name[M][B], **seq;
	static char **mseq1, **mseq2;
	static char **aseq;
	static char **bseq;
	static double **pscore;
	static double **eff;
	static double **node0, **node1;
	int i, j;
	int identity;
	static int ***topol;
	static double **len;
	FILE *prep;
	char c;
	int alloclen;

	arguments( argc, argv );


	seq = AllocateCharMtx( njob, nlenmax*3 );
	aseq = AllocateCharMtx( njob, nlenmax*3 );
	bseq = AllocateCharMtx( njob, nlenmax*3 );
	mseq1 = AllocateCharMtx( njob, 1 );
	mseq2 = AllocateCharMtx( njob, 1 );
	alloclen = nlenmax*3;

	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	pscore = AllocateDoubleMtx( njob, njob );
	eff = AllocateDoubleMtx( njob, njob );
	node0 = AllocateDoubleMtx( njob, njob );
	node1 = AllocateDoubleMtx( njob, njob );

	Read( name, nlen, seq );

	constants();

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illeagal character %c\n", c );
		exit( 1 );
	}

	writePre( njob, name, nlen, seq, 0 );

	if( tbutree == 0 )
	{
		for( i=1; i<njob; i++ ) 
		{
			if( nlen[i] != nlen[0] ) 
			{
				fprintf( stderr, "Input pre-aligned seqences or make hat2.\n" );
				exit( 1 );
			}
		}
		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
		{
		/*
			pscore[i][j] = (double)score_calc1( seq[i], seq[j] );
		*/
			pscore[i][j] = (double)substitution_hosei( seq[i], seq[j] );
		}
	}
	else
	{
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2( prep, njob, name, pscore );
		fclose( prep );
	}

	if( treemethod == 'x' )
		supg( njob, pscore, topol, len );
	else if( treemethod == 's' )
		spg( njob, pscore, topol, len );
	else if( treemethod == 'p' )
		upg2( njob, pscore, topol, len );
	else 
		ErrorExit( "Incorrect tree\n" );

	countnode( njob, topol, node0 );
	if( tbrweight )
	{
		weight = 3; 
		utree = 0; counteff( njob, topol, len, eff ); utree = 1;
	}
	else
	{
		for( i=0; i<njob; i++ ) eff[i][i] = 1.0;
	}


	for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );

	treebase( name, nlen, bseq, aseq, mseq1, mseq2, pscore, topol, len, eff, alloclen );

	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
	{
		pscore[i][j] = (double)substitution_hosei( aseq[i], aseq[j] );
	}
		spg( njob, pscore, topol, len );

	if( treemethod == 'x' )
		supg( njob, pscore, topol, len );
	else if( treemethod = 's' )
		spg( njob, pscore, topol, len );
	else if( treemethod = 'p' )
		upg2( njob, pscore, topol, len );
	else
		ErrorExit( "Incorrect tree\n" );

	countnode( njob, topol, node1 );

	identity = 1;
	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
		identity *= ( node1[i][j] == node0[i][j] );

	fprintf( trap_g, "done\n" );

	writePre( njob, name, nlen, aseq, 1 );
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif

	if( identity )
	{
		fprintf( stderr, "topology not changed\n" );
		return( 0 );
	}
	else 
	{
		fprintf( stderr, "topology changed\n" );
		return( -1 );
	}


}

int main( int argc, char *argv[] )
{
    int i, nlen;
    char b[B];
    char a[] = "=";
	int value;

    gets( b ); njob = atoi( b );
/*
	if( strstr( b, "ayhoff" ) ) scoremtx = 1;
	else if( strstr( b, "dna" ) || strstr( b, "DNA" ) ) scoremtx = -1; 
	else if( strstr( b, "M-Y" ) || strstr( b, "iyata" ) ) scoremtx = 2; 
	else scoremtx = 0;
*/
    nlenmax = 0;
    i = 0;
    while( i<njob )
    {
        gets( b );
        if( strncmp( b, a, 1 ) == 0 ) 
        {
            i++;
            gets( b ); nlen = atoi( b );
            if( nlen > nlenmax ) nlenmax = nlen;
        }
    }
    if( nlen > N || njob > M ) ErrorExit( "ERROR in main\n" );
    rewind( stdin );
    value = main1( argc, argv );

	SHOWVERSION;
	return( 0 );
}

