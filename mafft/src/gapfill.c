#include "mltaln.h"


void main1( int nlen[M], int argc, char *argv[] )
{
	int i;
	char **seq;
	static char lseq[N];
	static char name[M][B];
	static char gap[N];

	for( i=0; i<nlenmax; i++ ) gap[i] = '-';
	gap[nlenmax] = 0;
	seq = AllocateCharMtx( njob, nlenmax+10 );

	Read( name, nlen, seq );
	constants(); /* parameters no tsugou */

	for( i=0; i<njob; i++ ) 
	{
		strncpy( lseq, gap, 10 ); lseq[10] = 0;
		strcat( lseq, seq[i] );
		strncat( lseq, gap, nlenmax-nlen[i] ); lseq[nlenmax+10] = 0;
		strcpy( seq[i], lseq );
		nlen[i] = nlenmax+10;
	}

	Write( stdout, njob, name, nlen, seq );
}

void main( int argc, char *argv[] )
{
    int i, nlen[M];
    char b[B];
    char a[] = "=";

    gets( b ); njob = atoi( b );

    scoremtx = 0;
    for( i=0; i<B-5; i++ ) 
        if( !strncmp( b+10+i, "dayhoff", 8 ) )
        {
            scoremtx = 1;
            break;
        }

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
    main1( nlen, argc, argv );
	exit( 0 );
}

