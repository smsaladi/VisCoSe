#include "mltaln.h"

#define DEBUG 0

int ***AllocateTopol( int );
double Cscore_m_1_o( char **, int, double ** );
char **Calignm1();
char **Calignm1_o();

void arguments( int argc, char *argv[] )
{
	int c;

	calledByXced = 0;
	devide = 0;
    weight = 0;
    utree = 0;
    refine = 0;
    check = 1;
    cut = 2.0;
	disp = 0;
	outgap = 1;
	alg = 'C';
	mix = 0;

	while( --argc > 0 && (*++argv)[0] == '-' )
		while ( c = *++argv[0] )
			switch( c )
			{
				case 'Q':
					calledByXced = 1;
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
#if 0
				case 'A':
					alg = 'A';
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
		fprintf( stderr, "options: Check source file !\n" );
		exit( 1 );
	}
	(void)readOtherOptions( &ppid, &fftThreshold, &fftWinSize );
}

char **align0( float *wm, char **aseq, char *seq, double effarr[M], int icyc, int ex )
{
	char **result;

	if( alg == 'A' )
	{
		ErrorExit( "Sorry!" );
		
/*
		result = translate_and_Aalign( wm, aseq, seq, effarr, icyc, ex, alloclen );
*/
	}
	else if( alg == 'C' )
	{
		result = Calignm1( wm, aseq, seq, effarr, icyc, ex );
	}
	return( result );
}
	
void Writeoptions( FILE *fp )
{
	if( scoremtx )
    	fprintf( fp, "Gap Penalty = %d\n", penalty );
	else
    	fprintf( fp, "Gap Penalty = %d, %d\n", penalty, offset );
    fprintf( fp, "unweighted\n" );
	fprintf( fp, "Algorithm %c\n", alg );
	if( scoremtx == 1 ) fprintf( fp, "Dayhoff\n" );
	else if( scoremtx == -1 ) fprintf( fp, "DNA\n" );
	else if( scoremtx == 2 ) fprintf( fp, "M-Y\n" );
	else                     fprintf( fp, "JTT\n" );
}



void main1( int nlen[M], int argc, char *argv[] )
{
	static char name[M][B], **seq, **aseq, **bseq;
    static char *mseq1, **mseq2;
	char **result;
	int icyc;
	int i;
	int alloclen;
	int nlenmax0 = nlenmax;
	double effarr[M];
	float wm;
	char c;
	int pid;

	pid = getpid();

	arguments( argc, argv );


/*
	initDataSM( njob, nlenmax );
*/

	alloclen = nlenmax * 5;
	seq = AllocateCharMtx( njob, nlenmax+1 );
	aseq = AllocateCharMtx( njob, nlenmax*5+1 );
	bseq = AllocateCharMtx( njob+3, nlenmax*5+1 );
	AllocateTmpSeqs( &mseq2, &mseq1, alloclen );

	for( i=0; i<njob; i++ ) effarr[i] = 1.0;

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

	for( i=0; i<njob; i++ ) gappick0( aseq[i], seq[i] );
#if 1
	writePre( njob, name, nlen, aseq, 0 );
#else
	writeDataSM( njob, name, nlen, aseq, 0 );
#endif	

	for( icyc=0; icyc<njob-1; icyc++ ) 
	{
		result = align0( &wm, aseq, aseq[icyc+1], effarr, icyc, icyc );
#if DEBUG
		fprintf( stderr, "result = %d alloclen = %d\n", strlen( result[0] ), alloclen );
#endif
		if( strlen( result[0] ) > alloclen )   /* ???? */
		{
			if( strlen( result[0] ) > nlenmax0*3+1 )
			{
				fprintf(stderr, "Error in main1\n");
				exit( 1 );
			}
			FreeTmpSeqs( mseq2, mseq1 );
			alloclen = strlen( result[0] ) * 1.5;
			fprintf( stderr, "\n\ntrying to allocate TmpSeqs\n\n" );
			AllocateTmpSeqs( &mseq2, &mseq1, alloclen );
		}
		for( i=0; i<=icyc+1; i++ ) strcpy( aseq[i], result[i] );

		fprintf( stderr, "%3d / %3d \n", icyc+1, njob-1 );
		if( disp )
		{
			display( aseq );
			fprintf( stderr, "%3d / %3d \n", icyc+1, njob-1 );
		}
		if( strlen( aseq[0] ) > nlenmax )
			nlenmax = strlen( aseq[0] );
#if 1
		writePre( njob, name, nlen, aseq, 0 );
#else
		writeDataSM( njob, name, nlen, aseq, 0 );
#endif
	}

	for( i=0; i<njob; i++ ) strcpy( bseq[i], aseq[i] );
	Writeoptions( trap_g );
	fprintf( trap_g, "done" );

#if 1
	writePre( njob, name, nlen, aseq, 1 );
#else
	writeDataSM( njob, name, nlen, aseq, 0 );
#endif

	SHOWVERSION;
}

int main( int argc, char *argv[] )
{
	int i, nlen[M];
	char b[B];
	char a[] = "=";

	fgets( b, B-1, stdin ); njob = atoi( b );

/*
	scoremtx = 0;
	if( strstr( b, "dna" ) || strstr( b, "DNA" ) ) scoremtx = -1;
	else if( strstr( b, "ayhoff" ) ) scoremtx = 1;
	else if( strstr( b, "M-Y" ) ) scoremtx = 2;
	else scoremtx = 0;
*/

	if( strstr( b, "constraint" ) ) cnst = 1;
	nlenmax = 0;
	i = 0;
	while( i<njob )
	{
		fgets( b, B-1, stdin );
		if( strncmp( b, a, 1 ) == 0 ) 
		{
			i++;
			fgets( b, B-1, stdin ); nlen[i] = atoi( b );
			if( nlen[i] > nlenmax ) nlenmax = nlen[i];
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
#if DEBUG 
	printf( "njob = %d, nalenmax = %d\n", njob, nlenmax );
#endif
	main1( nlen, argc, argv );
	return( 0 );
}
