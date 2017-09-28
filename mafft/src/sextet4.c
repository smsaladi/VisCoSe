#include "mltaln.h"
#include "mtxutl.h"

#define DEBUG 0
#define TEST  0

#define END_OF_VEC -1

static int maxl;

void arguments( int argc, char *argv[] )
{
	int c;

	disopt = 0;

    while( --argc > 0 && (*++argv)[0] == '-' )
        while ( c = *++argv[0] )
            switch( c )
            {
				case 'i':
					disopt = 1;
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
    if( argc != 0 )
    {
        fprintf( stderr, "options: -i\n" );
        exit( 1 );
    }
}

void seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	while( *seq )
	{
		tmp = amino_grp[*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			fprintf( stderr, "WARNING : Unknown character %c\n", *(seq-1) );
	}
	*grp = END_OF_VEC;
}

void seq_grp( int *grp, char *seq )
{
	int tmp;
	while( *seq )
	{
		tmp = amino_grp[*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			fprintf( stderr, "WARNING : Unknown character %c\n", *(seq-1) );
	}
	*grp = END_OF_VEC;
}

int commonsextet_nuc_clear( short *table, int *n )
{
	register int *p;
	int value = 0;
	int i, j, k;
	short tmp;
	int point;
	static short *memo = NULL;
	static int *ct = NULL;
	static int *cp;

	if( !memo )
	{
		memo = calloc( 4096, sizeof( short ) );
		if( !memo ) ErrorExit( "Cannot allocate memo\n" );
		ct = calloc( MIN( maxl, 4096 ), sizeof( int ) );
		if( !ct ) ErrorExit( "Cannot allocate ct\n" );
	}

	cp = ct;
	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;

	tmp = memo[point]++;
	if( tmp < table[point] ) value++;
	*cp++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;

		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
	}
	*cp = END_OF_VEC;

	cp =  ct;
	while( *cp != END_OF_VEC )
		memo[*cp++] = 0;

	return( value );
}

int commonsextet_clear( short *table, int *n )
{
	register int *p;
	int value = 0;
	int i, j, k;
	short tmp;
	int point;
	static short *memo = NULL;
	static int *ct = NULL;
	static int *cp;

	if( !memo )
	{
		memo = calloc( 46656, sizeof( short ) );
		if( !memo ) ErrorExit( "Cannot allocate memo\n" );
		ct = calloc( MIN( 46656, maxl ), sizeof( int ) );
		if( !ct ) ErrorExit( "Cannot allocate memo\n" );
	}

	cp = ct;
	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;

	tmp = memo[point]++;
	if( tmp < table[point] ) value++;
	*cp++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;

#if 1
		if( point >= 46656 || point < 0 )
		{
			fprintf( stderr, "ERROR\n" );
			exit( 1 );
		}
#endif

		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
	}
	*cp = END_OF_VEC;
	
	cp =  ct;
	while( *cp != END_OF_VEC )
		memo[*cp++] = 0;

	return( value );
}
	
void makecompositiontable_nuc( short *table, int *n )
{
	int i, point;
	register int *p;

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	table[point]++;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		table[point]++;
	}
}

void makecompositiontable( short *table, int *n )
{
	int i, point;
	register int *p;

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	table[point]++;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		table[point]++;
	}
}

void main( int argc, char **argv )
{
	int i, j;
	FILE *fp;
	char **seq;
	int **grpseq;
	char *tmpseq;
	char name[M][B];
	int nlen[M];
	float **mtx;
	float **mtx2;
	float score, score0;
	static short *table1;
	char b[B];

#if 0
	PreRead( stdin, &njob, &nlenmax );
#else
	getnumlen( stdin );
#endif
	rewind( stdin );

	tmpseq = AllocateCharVec( nlenmax );
	seq = AllocateCharMtx( njob, nlenmax );
	grpseq = AllocateIntMtx( njob, nlenmax );
	mtx = AllocateFloatMtx( njob, njob );
	mtx2 = AllocateFloatMtx( njob, njob );

#if 0
	FRead( stdin, name, nlen, seq );
#else
	readData( stdin, name, nlen, seq );
#endif
	constants();


	if( scoremtx == -1 ) /* nuc */
	{
		maxl = 0;
		for( i=0; i<njob; i++ ) 
		{
			gappick0( tmpseq, seq[i] );
			nlen[i] = strlen( tmpseq );
			if( nlen[i] > maxl ) maxl = nlen[i];
			seq_grp_nuc( grpseq[i], tmpseq );
		}
		fprintf( stderr, "query :    " );
		for( i=0; i<njob; i++ ) 
		{
			table1 = calloc( 4096, sizeof( short ) );
            if( i % 10 == 0 )
            {
				fprintf( stderr, "\b\b\b\b\b\b\b\b\b\b\b" );
				fprintf( stderr, "%#4d / %#4d", i+1, njob ); 
            }
			makecompositiontable_nuc( table1, grpseq[i] );
	
			for( j=i; j<njob; j++ ) 
			{
				score = (float)commonsextet_nuc_clear( table1, grpseq[j] );
				mtx[i][j] = score;
			} 
			free( table1 );
		}
	}
	else
	{
		maxl = 0;
		for( i=0; i<njob; i++ ) 
		{
			gappick0( tmpseq, seq[i] );
			nlen[i] = strlen( tmpseq );
			if( nlen[i] > maxl ) maxl = nlen[i];
			seq_grp( grpseq[i], tmpseq );
		}
		fprintf( stderr, "query :    " );
		for( i=0; i<njob; i++ ) 
		{
			table1 = calloc( 46656, sizeof( short ) );
            if( i % 10 == 0 )
            {
                fprintf( stderr, "\b\b\b\b\b\b\b\b\b\b\b" );
				fprintf( stderr, "%#4d / %#4d", i+1, njob ); 
            }
			makecompositiontable( table1, grpseq[i] );
	
			for( j=i; j<njob; j++ ) 
			{
				score = (float)commonsextet_clear( table1, grpseq[j] );
				mtx[i][j] = score;
			} 
			free( table1 );
		}
	}
	for( i=0; i<njob; i++ )
	{
		score0 = mtx[i][i];
		for( j=0; j<njob; j++ ) 
			mtx2[i][j] = ( score0 - mtx[MIN(i,j)][MAX(i,j)] ) / score0 * 3.0;
	}
	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
	{
#if TEST
                double jscore;
                jscore = mtx[i][j] / ( MIN( strlen( seq[i] ), strlen( seq[j] ) ) - 2 );
                fprintf( stdout, "jscore = %f\n", jscore );

		fprintf( stdout, "mtx2[%d][%d] = %f, mtx2[%d][%d] = %f\n", i, j, mtx2[i][j], j, i, mtx2[j][i] );
#endif
		mtx2[i][j] = MIN( mtx2[i][j], mtx2[j][i] );
#if TEST
		fprintf( stdout, "sonokekka mtx2[%d][%d] %f\n", i, j, mtx2[i][j] );
#endif
	}

	if( disopt )
	{
		for( i=0; i<njob; i++ ) 
		{
			sprintf( b, "=lgth = %#04d", nlen[i] );
			strins( b, name[i] );
		}
	}
		
	fp = fopen( "hat2", "w" );
	if( !fp ) ErrorExit( "Cannot open hat2." );
#if 0
	WriteHat2( fp, njob, name, mtx2 );
#else
	WriteFloatHat2( fp, njob, name, mtx2 );
#endif
	fclose( fp );

	fprintf( stderr, "\b\b\b\b\b\b\b\b\b\b\b" );
	SHOWVERSION;
	exit( 0 );
}
