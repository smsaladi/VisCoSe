 /* iteration  ( algorithm C ) */
#include "mltaln.h"

double score_calc0( char **, int, double **, int );

void arguments( int argc, char *argv[] )
{
	int c;

	calledByXced = 0;
	devide = 0;
    scmtd = 5;
    weight = 0;
    utree = 0;
    refine = 1;
    check = 1;
    cut = 2.0;
	disp = 0;
	outgap = 1;
	alg = 'C';
	mix = 0;
	checkC = 0;
	tbitr = 0;

	while( --argc > 0 && (*++argv)[0] == '-' )
		while ( c = *++argv[0] )
			switch( c )
			{
				case 'Q':
					calledByXced = 1;
					break;
				case 'x':
					mix = 1;
					weight = 2;
					break;
				case 'C':
					alg = 'C';
					break;
				case 'l':
					scmtd = 5;
					break;
				case 'm':
					weight = 3;
					break;
				case 'w':
					weight = 2;
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
#if 0
				case 'o':
					outgap = 0;
					break;
				case 'B':
					alg = 'B';
					break;
				case 'p':
					check = -1;
					break;
				case 'n':
					scmtd = 4;
					break;
				case 'r':
					scmtd = 3;
					break;
#endif
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
	readOtherOptions( &ppid, &fftThreshold, &fftWinSize );
}


int main1( int nlen[M], int argc, char *argv[] )
{
    int identity;
	static char name[M][B], **seq, **aseq, **bseq;
	int i, j;
	int ***topol;
	double **len;
	double **eff;
	FILE *prep;
	char sai[M];
	char sai1[M];
	char sai2[M];
#if 0
    double mscore;
	char ssco1[N], ssco2[N], ssco3[N];
#endif
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

	for( i=0; i<njob; i++ ) 
	{
		sai1[i] = 1;
		sai2[i] = 2;
	}
	sai[njob] = sai1[njob] = sai2[njob] = 0;

	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	eff = AllocateDoubleMtx( njob, njob );
	seq = AllocateCharMtx( njob, nlenmax+1 );
	aseq = AllocateCharMtx( njob, nlenmax*3+1 );
	bseq = AllocateCharMtx( njob+3, nlenmax*3+1 );

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
		for( i=0; i<njob-1; i++ ) 
		{
			for( j=i+1; j<njob; j++ ) 
			{
				printf( " %f", eff[i][j] );
			}
			printf( "\n" );
		}
		fclose( prep );

/*
		for( i=0; i<njob-1; i++ ) for( j=i+1; i<njob; i++ )
		{
			eff[i][j] = 1.0 - exp( - eff[i][j] );
			eff[i][j] = - log( 1 - eff[i][j] )
		}
*/
			
		/*
		upg2( njob, eff, topol, len );
		*/
		spg( njob, eff, topol, len );
		counteff( njob, topol, len, eff );
	}
#if DEBUG
	printf( "utree = %d\n", utree );
#endif


	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );
	
	returnvalue = iteration( njob, name, nlen, aseq, bseq, topol, len, eff );

	if( weight == 0 ){ weight = 2; treeconstruction( bseq, njob, topol, len, eff ); weight = 0; }

#if 0
	sitescore( bseq, eff, ssco1, ssco2, ssco3 ); 
	bseq[njob+0] = ssco1;
	bseq[njob+1] = ssco2;
	bseq[njob+2] = ssco3;
	strcpy( name[njob+0], "=reliability ( unweighted )" );
	strcpy( name[njob+1], "=            ( weighted )" ); 
	strcpy( name[njob+2], "=            ( weighted, smoothed )" );
	Write( stdout, njob+3, name, nlen, bseq );
	printf( "score = %f\n", mscore );
#endif

	if( returnvalue == 0 ) fprintf( stderr, "converged\n" );
	else if( returnvalue == -1 ) fprintf( stderr, "oscillating\n" );
	else if( returnvalue == -2 ) fprintf( stderr, "iteration over \n" );
	fprintf( trap_g, "done\n" );
	fclose( trap_g );

	writePre( njob, name, nlen, bseq, 1 );

	SHOWVERSION;
	return( returnvalue );
}

int main( int argc, char *argv[] )
{
	int i, nlen[M];
	char b[B];
	char a[] = "=";

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
	return( main1( nlen, argc, argv ) );
}
