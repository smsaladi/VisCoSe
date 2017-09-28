#include "mltaln.h"
#include "mtxutl.h"

#define DEBUG 0
#define TEST  0

#define END_OF_VEC -1

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
	}
	*grp = END_OF_VEC;
}

int commonoctet_nuc( short *table, int *n )
{
	register int *p;
	int value = 0;
	int i, j, k;
	short tmp;
	int point;
	static short *memo = NULL;

	memo = calloc( 65536, sizeof( short ) );

	p = n;
	n += 8;
	point = p[0] * 16384 +
			p[1] *  4096 +
			p[2] *  1024 +
			p[3] *   256 +
			p[4] *    64 +
			p[5] *    16 +
			p[6] *     4 +
			p[7];
#if 0
		if( point >= 65536 || point < 0 )
		{
			fprintf( stderr, "ERROR\n" );
			exit( 1 );
		}
#endif

	tmp = memo[point]++;
	if( tmp < table[point] ) value++;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 16384;
		point *= 4;
		point += *n++;

		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
	}
	free( memo );
	return( value );
}

int commonoctet( short *table, int *n )
{
	register int *p;
	int value = 0;
	int i, j, k;
	short tmp;
	int point;
	static short *memo = NULL;

	memo = calloc( 1679616, sizeof( short ) );

	p = n;
	n += 8;
	point = p[0] * 279936 +
			p[1] *  46656 +
			p[2] *   7776 +
			p[3] *   1296 +
			p[4] *    216 +
			p[5] *     36 +
			p[6] *      6 +
			p[7];

	tmp = memo[point]++;
	if( tmp < table[point] ) value++;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 279936;
		point *= 6;
		point += *n++;

#if 0
		if( point >= 1679616 || point < 0 )
		{
			fprintf( stderr, "ERROR\n" );
			exit( 1 );
		}
#endif

		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
	}
	free( memo );
	return( value );
}

int commonoctet_clear( short *table, int *n )
{
	register int *p;
	int value = 0;
	int i, j, k;
	short tmp;
	int point;
	static short *memo = NULL;
	static int *ct = NULL;
	static int *cp;

	if( memo == NULL )
	{
		memo = calloc( 1679616, sizeof( short ) );
		ct = calloc( 1679616, sizeof( int ) );
	}

	cp = ct;
	p = n;
	n += 8;
	point = p[0] * 279936 +
			p[1] *  46656 +
			p[2] *   7776 +
			p[3] *   1296 +
			p[4] *    216 +
			p[5] *     36 +
			p[6] *      6 +
			p[7];

	tmp = memo[point]++;
	if( tmp < table[point] ) value++;

	*cp++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 279936;
		point *= 6;
		point += *n++;

#if 0
		if( point >= 1679616 || point < 0 )
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

	n += 7;
	while( *(++n) != END_OF_VEC )
	{
		point = n[-7] * 16384 +
				n[-6] *  4096 +
				n[-5] *  1024 +
				n[-4] *   256 +
				n[-3] *    64 +
				n[-2] *    16 +
				n[-1] *     4 +
				n[ 0];
#if 0
		if( point >= 65536 || point < 0 )
		{
			fprintf( stderr, "ERROR point = %d\n", point );
			exit( 1 );
		}
#endif
		table[point]++;
	}
}

void makecompositiontable( short *table, int *n )
{
	int i, point;

	n += 7;
	while( *(++n) != END_OF_VEC )
	{
		point = n[-7] *279936 +
				n[-6] * 46656 +
				n[-5] *  7776 +
				n[-4] *  1296 +
				n[-3] *   216 +
				n[-2] *    36 +
				n[-1] *     6 +
				n[ 0];
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
	double **mtx;
	double **mtx2;
	double score, score0;
	static short *table1;
	char b[B];

	getnumlen( stdin );
	rewind( stdin );

	tmpseq = AllocateCharVec( nlenmax );
	seq = AllocateCharMtx( njob, nlenmax );
	grpseq = AllocateIntMtx( njob, nlenmax );
	mtx = AllocateDoubleMtx( njob, njob );
	mtx2 = AllocateDoubleMtx( njob, njob );

	readData( stdin, name, nlen, seq );

	constants();

	for( i=0; i<njob; i++ )
	{
		gappick0( tmpseq, seq[i] );
#if 0
		fprintf( stderr, "%s\n", tmpseq );
#endif
		nlen[i] = strlen( tmpseq );
	}


	fprintf( stderr, "query :	" );
	if( scoremtx == -1 ) /* nuc */
	{
		for( i=0; i<njob; i++ ) 
			seq_grp_nuc( grpseq[i], tmpseq );
		for( i=0; i<njob; i++ ) 
		{
			table1 = calloc( 65536, sizeof( short ) );
			fprintf( stderr, "\b\b\b\b" );
			fprintf( stderr, "%#4d", i+1 );
			makecompositiontable_nuc( table1, grpseq[i] );
	
			for( j=i; j<njob; j++ ) 
			{
				score = (double)commonoctet_nuc( table1, grpseq[j] );
				mtx[i][j] = score;
			} 
			free( table1 );
		}
	}
	else
	{
		for( i=0; i<njob; i++ ) 
			seq_grp( grpseq[i], tmpseq );
		for( i=0; i<njob; i++ ) 
		{
			table1 = calloc( 1679616, sizeof( short ) );
			fprintf( stderr, "\b\b\b\b" );
			fprintf( stderr, "%#4d", i+1 );
			makecompositiontable( table1, grpseq[i] );
	
			for( j=i; j<njob; j++ ) 
			{
				score = (double)commonoctet_clear( table1, grpseq[j] );
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
	WriteHat2( fp, njob, name, mtx2 );
	fclose( fp );

	fprintf( stderr, "\b\b\b\b\b\b\b\b\b\b\b" );
	SHOWVERSION;
	exit( 0 );
}
