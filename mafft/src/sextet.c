#include "mltaln.h"
#include "mtxutl.h"

#define DEBUG 0
#define TEST  0

void seq_grp( int *grp, char *seq )
{
	while( *seq )
		*grp++ = amino_grp[*seq++];
	*grp = -1;
}

void makeGroupCompositionTable( short table[7][7][7][7][7][7], int *seq )
{
	int i, j, k, l, m, n;
	int len = intlen( seq );

#if DEBUG
	for( i=0; i<len; i++ ) fprintf( stderr, "seq[%d] = %d\n", i, seq[i] );
#endif

	for( i=0; i<7; i++ ) for( j=0; j<7; j++ ) for( k=0; k<7; k++ ) for( l=0; l<7; l++ ) for( m=0; m<7; m++ ) for( n=0; n<7; n++ )   /* ichiou */
		table[i][j][k][l][m][n] = 0;
	for( i=0; i<len-5; i++ )
		table[seq[i]][seq[i+1]][seq[i+2]][seq[i+3]][seq[i+4]][seq[i+5]]++;
}

int commonSextet( short table[7][7][7][7][7][7], short memo[7][7][7][7][7][7], int *seq )
{
	int i, j, k, l, m, n, tmp;
	int value = 0;
	int len = intlen( seq );

#if 1
	for( i=0; i<6; i++ ) for( j=0; j<6; j++ ) for( k=0; k<6; k++ ) for( l=0; l<6; l++ ) for( m=0; m<6; m++ ) for( n=0; n<6; n++ ) 
		memo[i][j][k][l][m][n] = 0;
	for( i=0; i<len-5; i++ )
	{
		tmp = memo[*seq][*(seq+1)][*(seq+2)][*(seq+3)][*(seq+4)][*(seq+5)]++;
		if( tmp < table[*seq][*(seq+1)][*(seq+2)][*(seq+3)][*(seq+4)][*(seq+5)] )
			value++;
		seq++;
	}
#else
	makeGroupCompositionTable( memo, seq );
	for( i=0; i<6; i++ ) for( j=0; j<6; j++ ) for( k=0; k<6; k++ ) for( l=0; l<6; l++ ) for( m=0; m<6; m++ ) for( n=0; n<6; n++ ) 
		value += MIN( table[i][j][k][l][m][n], memo[i][j][k][l][m][n] );
#endif
	return( value );
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
	double score;
	static short table1[7][7][7][7][7][7];
	static short table2[7][7][7][7][7][7];
	double byChance;

#if 0
	PreRead( stdin, &njob, &nlenmax );
#else
	getnumlen( stdin );
#endif
	rewind( stdin );

	tmpseq = AllocateCharVec( nlenmax );
	seq = AllocateCharMtx( njob, nlenmax );
	grpseq = AllocateIntMtx( njob, nlenmax );
	mtx = AllocateDoubleMtx( njob, njob );

#if 0
	FRead( stdin, name, nlen, seq );
#else
	readData( stdin, name, nlen, seq );
#endif

	constants();

	if( scoremtx == -1 ) 
		byChance = pow( 0.25, 6 );
	else
		byChance = 0.000043083;
#if TEST
	fprintf( stdout, "byChance = %f\n", byChance );
#endif

	for( i=0; i<njob; i++ )
	{
		gappick0( tmpseq, seq[i] );
#if DEBUG
		fprintf( stderr, "%s\n", tmpseq );
#endif
		seq_grp( grpseq[i], tmpseq );
		nlen[i] = strlen( tmpseq );
	}

	fprintf( stderr, "query :    " );
	for( i=0; i<njob-1; i++ ) 
	{
		fprintf( stderr, "\b\b\b" );
		fprintf( stderr, "%#3d", i+1 );
		makeGroupCompositionTable( table1, grpseq[i] );
		for( j=i+1; j<njob; j++ ) 
		{
			score = (double)commonSextet( table1, table2, grpseq[j] );
			score -= byChance * ( nlen[i] - 5 ) * ( nlen[j] - 5 );
			score /= (double)( MIN( nlen[i], nlen[j] ) - 5 );
#if TEST
			fprintf( stdout, "sore = %f\n", log( score ) );
#endif
			if( score > 0.0 ) 
				score = -0.293 * log( score );
			else score = 1.0;
			if     ( score < 0.0  ) score = 0.0;
			else if( score < 0.95 ) score = -log( 1.0 - score );
			else                    score = 3.0;
			mtx[i][j] = score;
		} 
	}
	fp = fopen( "hat2", "w" );
	if( !fp ) ErrorExit( "Cannot open hat2." );
	WriteHat2( fp, njob, name, mtx );
	fclose( fp );

	fprintf( stderr, "\b\b\b\b\b\b\b\b\b\b\b" );
	SHOWVERSION;
	exit( 0 );
}
