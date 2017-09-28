#include "mltaln.h"

#define DEBUG 0

void main1( int nlen[M], int argc, char *argv[] )
{
	static char name[M][B], **seq;
	int i, j;
	FILE *prep;
	int nlenmax0 = nlenmax;
	int pid, signalsmid;

	char amino_clw[] = "CSTPAGNDEQHRKMILVFYW";

	seq = AllocateCharMtx( njob, nlenmax+1 );

	Read( name, nlen, seq );

	constants();

	for( i=0; i<20; i++ )
	{
		for( j=0; j<=i; j++ )
		{
			fprintf( stdout, "%3.0f", (double)amino_dis[amino_clw[i]][amino_clw[j]]/50.0 );
		}
		fprintf( stdout, "\n" );
	}

	SHOWVERSION;
}

void main( int argc, char *argv[] )
{
	int i, nlen[M];
	char b[B];
	char a[] = "=";

	gets( b ); njob = atoi( b );

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
		gets( b );
		if( strncmp( b, a, 1 ) == 0 ) 
		{
			i++;
			gets( b ); nlen[i] = atoi( b );
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
}
