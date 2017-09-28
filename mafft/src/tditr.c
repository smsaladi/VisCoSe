 /* Tree-dependent-iteration */

#include "mltaln.h"

void arguments( int argc, char *argv[] )
{
	int c;

	calledByXced = 0;
	devide = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    scmtd = 5;
	cooling = 0;
    weight = 4;
    utree = 0;
    refine = 1;
    check = 1;
    cut = 0.0;
	disp = 0;
	outgap = 1;
	use_fft = 0;
	alg = 'C';  /* chuui */
	mix = 0;
	checkC = 0;
	tbitr = 0;
	treemethod = 'x';

	while( --argc > 0 && (*++argv)[0] == '-' )
		while ( c = *++argv[0] )
			switch( c )
			{
				case 'Q':
					calledByXced = 1;
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
				case 'c':
					cooling = 1;
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
				case 'F':
					use_fft = 1;
					break;
				case 'w':
					weight = 2;
					break;
				case 't':
					weight = 4;
					break;
				case 'u':
					weight = 0;
					break;
				case 'j':
					utree = 1;
					break;
				case 'd':
					disp = 1;
					break;
				case 'Z':
					checkC = 1;
					break;
				case 'n' :
					treemethod = 'n';
					break;
				case 's' :
					treemethod = 's';
					break;
				case 'x' :
					treemethod = 'x';
					break;
				case 'p' :
					treemethod = 'p';
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
		fprintf( stderr, "options : Check source file!\n" );
		exit( 1 );
	}
#if 0
	if( alg == 'A' && weight == 0 ) 
		ErrorExit( "ERROR : Algorithm A+ and un-weighted\n" ); 
#endif
	readOtherOptions( &ppid, &fftThreshold, &fftWinSize );
}


static int main1( int nlen[M], int argc, char *argv[] )
{
    int identity;
	static char name[M][B], **seq, **aseq, **bseq;
	int i, j;
	int ***topol;
	double **len;
	double **eff;
	FILE *prep;
	int alloclen;
	int returnvalue;
	char c;

	identity = 1;
	for( i=0; i<njob; i++ ) 
	{
		identity *= ( nlen[i] == nlen[0] );
	}
	if( !identity ) 
	{
		fprintf( stderr, "Input pre-aligned data\n" );
		exit( 1 );
	}

	arguments( argc, argv );

	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	eff = AllocateDoubleMtx( njob, njob );
	seq = AllocateCharMtx( njob, nlenmax+1 );
	aseq = AllocateCharMtx( njob, nlenmax*5+1 );
	bseq = AllocateCharMtx( njob, nlenmax*5+1 );
	alloclen = nlenmax * 5;

	Read( name, nlen, seq );
	constants();

	initSignalSM();

	initFiles();


	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illeagal character %c\n", c );
		exit( 1 );
	}
	commongappick( njob, seq );

	if( utree )
	{
		prep = fopen( "hat2", "r" );
		if( !prep ) ErrorExit( "Make hat2." );
		readhat2( prep, njob, name, eff );
		fclose( prep );
#if DEBUG
		for( i=0; i<njob-1; i++ ) 
		{
			for( j=i+1; j<njob; j++ ) 
			{
				printf( " %f", eff[i][j] );
			}
			printf( "\n" );
		}
#endif
		if     ( treemethod == 'x' ) 
			supg( njob, eff, topol, len );
		else if( treemethod == 'n' ) 
			nj( njob, eff, topol, len );
		else if( treemethod == 's' )
			spg( njob, eff, topol, len );
		else if( treemethod == 'p' )
			upg2( njob, eff, topol, len );
		else ErrorExit( "Incorrect treemethod.\n" );
	}
#if DEBUG
	printf( "utree = %d\n", utree );
#endif

	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );
	
	returnvalue = TreeDependentIteration( njob, name, nlen, aseq, bseq, topol, len, alloclen );

#if 0
	Write( stdout, njob, name, nlen, bseq );
#endif

	if( returnvalue == 0 ) fprintf( stderr, "converged\n" );
	else if( returnvalue == -1 ) fprintf( stderr, "oscillating\n" );
	fprintf( trap_g, "done\n" );
	fclose( trap_g );

	writePre( njob, name, nlen, bseq, 1 );

	SHOWVERSION;
	return( returnvalue );
}

signed int main( int argc, char *argv[] )
{
	int i, nlen[M];
	char b[B];
	char a[] = "=";
	int value;

	gets( b ); njob = atoi( b );

/*
	scoremtx = 0;
	if( strstr( b, "ayhoff" ) ) scoremtx = 1;
	else if( strstr( b, "dna" ) || strstr( b, "DNA" ) ) scoremtx = -1;
	else if( strstr( b, "M-Y" ) || strstr( b, "iyata" ) ) scoremtx = 2;
	else scoremtx = 0;
*/
	if( strstr( b, "constraint" ) ) cnst = 1;

	nlenmax = 0;
	i = 0;
	while( i<njob )
	{
		gets( b );
		if( !strncmp( b, a, 1 ) ) 
		{
			gets( b ); nlen[i] = atoi( b );
			if( nlen[i] > nlenmax ) nlenmax = nlen[i];
			i++;
		}
	}
	if( nlenmax > N || njob > M ) 
	{
		fprintf( stderr, "ERROR in main\n" );
		exit( 1 );
	}
	/*
	nlenmax = Na;
	*/
	rewind( stdin );
	value = main1( nlen, argc, argv );
	exit( 0 );
}
