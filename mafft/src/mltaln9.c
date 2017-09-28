#include "mltaln.h"

#define DEBUG 0

int intlen( int *num )
{
	int value = 0;
	while( *num++ != -1 ) value++;
	return( value );
}

char seqcheck( char **seq )
{
	int i, len;
	while( *seq )	
	{
		len = strlen( *seq );
		for( i=0; i<len; i++ ) 
			if( amino_n[(int)(*seq)[i]] == -1 ) return( (int)(*seq)[i] );
		seq++;
	}
	return( 0 );
}
	
void scmx_calc( int icyc, char **aseq, double *effarr, float **scmx )
{
	int  i, j, lgth;
	 
	lgth = strlen( aseq[0] );
	for( j=0; j<lgth; j++ )
	{
		for( i=0; i<26; i++ )
		{
			scmx[i][j] = 0;
		}
	}
	for( i=0; i<icyc+1; i++ )
	{
		int id;
		id = amino_n[(int)aseq[i][0]];
		scmx[id][0] += (float)effarr[i];
	}
	for( j=1; j<lgth-1; j++ )
	{
		for( i=0; i<icyc+1; i++ )
		{
			int id;
			id = amino_n[(int)aseq[i][j]];
			scmx[id][j] += (float)effarr[i];
		}
	}
	for( i=0; i<icyc+1; i++ )
	{
		int id;
		id = amino_n[(int)aseq[i][lgth-1]];
		scmx[id][lgth-1] += (float)effarr[i];
	}
}

void exitall( char arr[] )
{
	fprintf( stderr, "%s\n", arr );
	exit( 1 );
}

void display( char **seq, int nseq )
{
	int i, imax;
	char b[121];

	if( !disp ) return;
		if( nseq > DISPSEQF ) imax = DISPSEQF;
		else                  imax = nseq;
		fprintf( stderr, "    ....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+\n" );
		for( i=0; i<+imax; i++ )
		{
			strncpy( b, seq[i]+DISPSITEI, 120 );
			b[120] = 0;
			fprintf( stderr, "%3d %s\n", i+1, b );
		}
}

double intergroup_score( char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2 )
{
	int i, j, k;
	int len = strlen( seq1[0] );
	double score;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

	score = 0.0;
	for( i=0; i<clus1; i++ ) for( j=0; j<clus2; j++ ) 
	{
		efficient = eff1[i] * eff2[j];
		mseq1 = seq1[i];
		mseq2 = seq2[j];
		tmpscore = 0.0;
		for( k=0; k<len; k++ ) 
		{
			if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
			tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];

			if( mseq1[k] == '-' ) 
			{
				tmpscore += penalty;
				tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
				while( mseq1[++k] == '-' )
					tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
				k--;
				if( k >len-2 ) break;
				continue;
			}
			if( mseq2[k] == '-' )
			{
				tmpscore += penalty;
				tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
				while( mseq2[++k] == '-' )
					tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
				k--;
				if( k > len-2 ) break;
				continue;
			}
		}
		score += (double)tmpscore * efficient;
	}
#if DEBUG
	fprintf( stderr, "score in intergroup_score = %f\n", score );
#endif
	return( score );
}

double score_calc3( char **seq, int s, double **eff, int ex )  /* method 3 */
{
    int i, j, k, c;
    int len = strlen( seq[0] );
    double score;
    long tmpscore;
    static char mseq1[N*2], mseq2[N*2];
	double totaleff;

	switch( weight )
	{
		case 0:
			totaleff = ( (double)s * ((double)s-1.0) ) / 2.0;
			break;
		case 2:
			totaleff = s / 2; 
			break;
		case 3:
			totaleff = 0.0; for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ ) totaleff += eff[i][j];
			break;
		default:
			fprintf( stderr, "error\n" );
			exit( 1 );
	}

    score = 0.0;
    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
            strcpy( mseq1, seq[i] );
            strcpy( mseq2, seq[j] );
            tmpscore = 0;
            c = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]] + 400 * !scoremtx;
                c++;
                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty - n_dis[0][24];
                    while( mseq1[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty - n_dis[0][24];
                    while( mseq2[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            /*
            if( mseq1[0] == '-' || mseq2[0] == '-' )
            {
                for( k=0; k<len; k++ )
                {
                    if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                    if( !( mseq1[k] != '-' && mseq2[k] != '-' ) )
                    {
                        c--;
                        tmpscore -= SGAPP;
                        break;
                    }
                    else break;
                }
            }
            if( mseq1[len-1] == '-' || mseq2[len-1] == '-' )
            {
                for( k=len-1; k>=0; k-- )
                {
                    if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                    if( !( mseq1[k] != '-' && mseq2[k] != '-' ) )
                    {
                        c--;
                        tmpscore -= SGAPP;
                        break;
                    }
                    else break;
                }
            }
            */
            /*
            if( x == 65 ) printf( "i=%d j=%d tmpscore=%d l=%d\n", i, j, tmpscore, len )
;
            */
            score += (double)tmpscore / (double)c * eff[i][j];
        }
    }
	if( weight == 0 )
    	score /= totaleff; 
    return( (double)score );
}
double score_calc5( char **seq, int s, double **eff, int ex )  /* method 3 deha nai */
{
    int i, j, k;
    double c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;
    double efficient;
#if DEBUG
	FILE *fp;
#endif

    score = 0.0;
    c = 0.0;

	for( i=0; i<s; i++ ) 
	{
		
			if( i == ex ) continue;
            efficient = eff[i][ex];
            mseq1 = seq[i];
            mseq2 = seq[ex];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore * efficient;
/*
			fprintf( stdout, "%d-%d tmpscore = %f, eff = %f, tmpscore*eff = %f\n", i, ex, tmpscore, efficient, tmpscore*efficient );
*/
	}
	/*
	fprintf( stdout, "total score = %f\n", score );
	*/

    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
			if( i == ex || j == ex ) continue;

            efficient = eff[i][j];
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore * efficient;
        }
    }
/*
	fprintf( stderr, "score in score_calc5 = %f\n", score );
*/
    return( (double)score );
/*

fprintf( trap_g, "score by fast = %f\n", (float)score );

tmpscore = score = 0.0;
	for( i=0; i<s; i++ ) 
	{
		if( i == ex ) continue;
		tmpscore = Cscore_m_1( seq, i, eff );
		fprintf( stdout, "%d %f\n", i, tmpscore );

		score += tmpscore;
	}
	tmpscore = Cscore_m_1( seq, ex, eff );
	fprintf( stdout, "ex%d %f\n", i, tmpscore );
	score += tmpscore;

	return( score );
*/
}


	
double score_calc4( char **seq, int s, double **eff, int ex )  /* method 3 deha nai */
{
    int i, j, k;
	double c;
    int len = strlen( seq[0] );
    double score;
    long tmpscore;
	char *mseq1, *mseq2;
	double efficient;

    score = 0.0;
	c = 0.0;
/*
	printf( "in score_calc4\n" );
	for( i=0; i<s; i++ )
	{
		for( j=0; j<s; j++ ) 
		{
			printf( "% 5.3f", eff[i][j] ); 
		}
		printf( "\n" );
		
	}
*/
    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
			efficient = eff[i][j];
			if( mix == 1 ) efficient = 1.0;
			/*
			printf( "weight for %d v.s. %d = %f\n", i, j, efficient );
			*/
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]] + 400 * !scoremtx ;

				c += efficient;

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty - n_dis[24][0];
                    while( mseq1[++k] == '-' )
						;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty - n_dis[24][0];
                    while( mseq2[++k] == '-' )
						;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
			/*
			if( x == 65 ) printf( "i=%d j=%d tmpscore=%d l=%d\n", i, j, tmpscore, len );
			*/
            score += (double)tmpscore * efficient;
        }
    }
    score /= c;
    return( (double)score );
}



void upg2( int nseq, double **eff, int ***topol, double **len )
{
    int i, j, k;
	double tmplen[M];

    static char **pair = NULL;

	if( !pair )
	{
		pair = AllocateCharMtx( njob, njob );
	}

	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

    for( k=0; k<nseq-1; k++ )
    {
        float minscore = 9999.0;
        int im = -1, jm = -1;
        int count;

        for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
        {
            if( eff[i][j] < minscore )
            {
                minscore = eff[i][j];
                im = i; jm = j;
            }
        }
        for( i=0, count=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
            {
                topol[k][0][count] = i;
                count++;
            }
        topol[k][0][count] = -1;
        for( i=0, count=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
            {
                topol[k][1][count] = i;
                count++;
            }
        topol[k][1][count] = -1;

		len[k][0] = minscore / 2.0 - tmplen[im];
		len[k][1] = minscore / 2.0 - tmplen[jm];

		tmplen[im] = minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
                eff[MIN(i,im)][MAX(i,im)] =
                ( eff[MIN(i,im)][MAX(i,im)] + eff[MIN(i,jm)][MAX(i,jm)] ) / 2.0;
                eff[MIN(i,jm)][MAX(i,jm)] = 9999.0;
            }
            eff[im][jm] = 9999.0;
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

void supg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k;
#if 0
	double eff[nseq][nseq];
    char pair[njob][njob];
#else
	static double *tmplen;
	static double **eff = NULL;
    static char **pair = NULL;
	if( !eff )
	{
		eff = AllocateDoubleMtx( njob, njob );
		pair = AllocateCharMtx( njob, njob );
		tmplen = AllocateDoubleVec( njob );
	}
#endif

	
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			eff[i][j] = oeff[i][j];
		}
	}
	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;


    for( k=0; k<nseq-1; k++ )
    {
        float minscore = 9999.0;
        int im = -1, jm = -1;
        int count;


        for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
        {
            if( eff[i][j] < minscore )
            {
                minscore = eff[i][j];
                im = i; jm = j;
            }
        }
        for( i=0, count=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
            {
                topol[k][0][count] = i;
                count++;
            }
        topol[k][0][count] = -1;
        for( i=0, count=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
            {
                topol[k][1][count] = i;
                count++;
            }
        topol[k][1][count] = -1;

		len[k][0] = minscore / 2.0 - tmplen[im];
		len[k][1] = minscore / 2.0 - tmplen[jm];

		tmplen[im] = minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
                eff[MIN(i,im)][MAX(i,im)] =
                MIN( eff[MIN(i,im)][MAX(i,im)], eff[MIN(i,jm)][MAX(i,jm)] ) * ( 1.0 - SUEFF ) +
				( eff[MIN(i,im)][MAX(i,im)] + eff[MIN(i,jm)][MAX(i,jm)] ) * 0.5 * SUEFF;
                eff[MIN(i,jm)][MAX(i,jm)] = 9999.0;
            }
            eff[im][jm] = 9999.0;
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

void spg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k;
	double tmplen[M];
#if 0
	double eff[nseq][nseq];
    char pair[njob][njob];
#else
	double **eff = NULL;
    char **pair = NULL;
	if( !eff )
	{
		eff = AllocateDoubleMtx( njob, njob );
		pair = AllocateCharMtx( njob, njob );
	}
#endif
	
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) eff[i][j] = oeff[i][j];
	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

    for( k=0; k<nseq-1; k++ )
    {
        float minscore = 9999.0;
        int im = -1, jm = -1;
        int count;

        for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
        {
            if( eff[i][j] < minscore )
            {
                minscore = eff[i][j];
                im = i; jm = j;
            }
        }
        for( i=0, count=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
            {
                topol[k][0][count] = i;
                count++;
            }
        topol[k][0][count] = -1;
        for( i=0, count=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
            {
                topol[k][1][count] = i;
                count++;
            }
        topol[k][1][count] = -1;

		len[k][0] = minscore / 2.0 - tmplen[im];
		len[k][1] = minscore / 2.0 - tmplen[jm];

		tmplen[im] = minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
                eff[MIN(i,im)][MAX(i,im)] =
                MIN( eff[MIN(i,im)][MAX(i,im)], eff[MIN(i,jm)][MAX(i,jm)] );
                eff[MIN(i,jm)][MAX(i,jm)] = 9999.0;
            }
            eff[im][jm] = 9999.0;
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

double ipower( double x, int n )    /* n > 0  */
{
    double r;

    r = 1;
    while( n != 0 )
    {
        if( n & 1 ) r *= x;
        x *= x; n >>= 1;
    }
    return( r );
}

void countnode( int nseq, int ***topol, double **node ) /* node[j][i] != node[i][j] */
{
    int i, j, k, s1, s2;
    double rootnode[M];

    if( nseq-2 < 0 )
	{
		fprintf( stderr, "Too few sequence for countnode: nseq = %d\n", nseq );
		exit( 1 );
    }

    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    for( i=0; i<nseq-2; i++ )
    {
        for( j=0; topol[i][0][j]>-1; j++ )
            rootnode[topol[i][0][j]]++;
        for( j=0; topol[i][1][j]>-1; j++ )
            rootnode[topol[i][1][j]]++;
        for( j=0; topol[i][0][j]>-1; j++ )
        {
            s1 = topol[i][0][j];
            for( k=0; topol[i][1][k]>-1; k++ )
            {
                s2 = topol[i][1][k];
                node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            }
        }
    }
    for( j=0; topol[nseq-2][0][j]>-1; j++ )
    {
        s1 = topol[nseq-2][0][j];
        for( k=0; topol[nseq-2][1][k]>-1; k++ )
        {
            s2 = topol[nseq-2][1][k];
            node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        }
    }
}

void countnode_int( int nseq, int ***topol, int **node )  /* node[i][j] == node[j][i] */
{
    int i, j, k, s1, s2;
    int rootnode[M];

    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    for( i=0; i<nseq-2; i++ )
    {
        for( j=0; topol[i][0][j]>-1; j++ )
            rootnode[topol[i][0][j]]++;
        for( j=0; topol[i][1][j]>-1; j++ )
            rootnode[topol[i][1][j]]++;
        for( j=0; topol[i][0][j]>-1; j++ )
        {
            s1 = topol[i][0][j];
            for( k=0; topol[i][1][k]>-1; k++ )
            {
                s2 = topol[i][1][k];
                node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            }
        }
    }
    for( j=0; topol[nseq-2][0][j]>-1; j++ )
    {
        s1 = topol[nseq-2][0][j];
        for( k=0; topol[nseq-2][1][k]>-1; k++ )
        {
            s2 = topol[nseq-2][1][k];
            node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        }
    }
	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ ) 
		node[j][i] = node[i][j];
#if DEBUG
	fprintf( stderr, "node[][] in countnode_int" );
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			fprintf( stderr, "%#3d", node[i][j] );
		}
		fprintf( stderr, "\n" );
	}
#endif
}

void counteff( int nseq, int ***topol, double **len, double **node )
{
    int i, j, k, s1, s2;
	double rootnode[M];
	double eff[M];

	if( mix ) 
	{
		switch( weight )
		{
			case( 2 ): 
				weight = 3;
				break;
			case( 3 ): 
				weight = 2;
				break;
			default: 
				ErrorExit( "mix error" );
				break;
		}
	}

	if( weight == 2 )
	{
	    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    	for( i=0; i<nseq-2; i++ )
    	{
        	for( j=0; topol[i][0][j]>-1; j++ )
            	rootnode[topol[i][0][j]]++;
        	for( j=0; topol[i][1][j]>-1; j++ )
            	rootnode[topol[i][1][j]]++;
        	for( j=0; topol[i][0][j]>-1; j++ )
        	{
            	s1 = topol[i][0][j];
            	for( k=0; topol[i][1][k]>-1; k++ )
            	{
                	s2 = topol[i][1][k];
                	node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            	}
        	}
    	}
    	for( j=0; topol[nseq-2][0][j]>-1; j++ )
    	{
        	s1 = topol[nseq-2][0][j];
        	for( k=0; topol[nseq-2][1][k]>-1; k++ )
        	{
            	s2 = topol[nseq-2][1][k];
            	node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        	}
    	}
   		for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
   	   		node[i][j] = ipower( 0.5, (int)node[i][j] ) + geta2;
		for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ ) 
			node[j][i] = node[i][j];
	}

	if( weight == 3 )
	{
#if DEBUG
		for( i=0; i<nseq; i++ ){
			fprintf( stderr, "len0 = %f\n", len[i][0] );
			fprintf( stderr, "len1 = %f\n", len[i][1] );
		}
#endif
	    for( i=0; i<nseq; i++ )
		{
			rootnode[i] = 0.0;
			eff[i] = 1.0;
/*
			rootnode[i] = 1.0;
*/
		}
    	for( i=0; i<nseq-1; i++ )
    	{
        	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
			{
   	        	rootnode[s1] += len[i][0] * eff[s1];
				eff[s1] *= 0.5;
/*
   	        	rootnode[s1] *= 0.5;
*/
				
			}
   	    	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
			{
   	        	rootnode[s2] +=  len[i][1] * eff[s2];
				eff[s2] *= 0.5;
/*
   	        	rootnode[s2] *= 0.5;
*/
				
			}
		}
		for( i=0; i<nseq; i++ ) 
		{
#if 1 /* 97.9.29 */
			rootnode[i] += GETA3;
#endif
#if DEBUG
			fprintf( stderr, "rootnode for %d = %f\n", i, rootnode[i] );
#endif
		}
		for( i=0; i<nseq; i++ ) 
		{
			for( j=0; j<nseq; j++ ) 
				if( j != i )
					node[i][j] = (double)rootnode[i] * rootnode[j];
				else node[i][i] = rootnode[i];
		}
	}

#if 0
	printf( "weight matrix in counteff\n" );
	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<nseq; j++ ) 
		{
			printf( "%f ", node[i][j] );
		}
		printf( "\n" );
	}
#endif
}

float score_calc1( char *seq1, char *seq2 )   /* method 1 */
{
	int k;
	float score = 0.0;
	int count = 0;
	int len = strlen( seq1 );

	for( k=0; k<len; k++ )
	{	
		if( seq1[k] != '-' && seq2[k] != '-' )
		{
			score += (float)amino_dis[(int)seq1[k]][(int)seq2[k]];
			count++;
		}
	}
	if( count ) score /= (float)count;
	else score = 1.0;
	return( score );
}

float substitution_hosei( char *seq1, char *seq2 )   /* method 1 */
{
	int k;
	float score = 0.0;
	int count = 0;
	int len = strlen( seq1 );

	for( k=0; k<len; k++ )
	{	
		if( seq1[k] != '-' && seq2[k] != '-' )
		{
			score += (float)( seq1[k] != seq2[k] );
			count++;
		}
	}
	if( count ) score /= (float)count;
	else score = 1.0;
	if( score < 0.95 ) score = - log( 1.0 - score );
	else score = 3.0;
	return( score );
}

float substitution( char *seq1, char *seq2 )   /* method 1 */
{
	int k;
	float score = 0.0;
	int count = 0;
	int len = strlen( seq1 );

	for( k=0; k<len; k++ )
	{	
		if( seq1[k] != '-' && seq2[k] != '-' )
		{
			score += (float)( seq1[k] != seq2[k] );
			count++;
		}
	}
	if( count ) score /= (float)count;
	else score = 1.0;
	return( score );
}


void treeconstruction( char **seq, int nseq, int ***topol, double **len, double **eff )
{
    int i, j;

	if( weight > 1 )
	{
		if( utree == 0 )
		{
	    	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
   		 	{
/*
		       	 eff[i][j] = (double)score_calc1( seq[i], seq[j] );
*/
		       	 eff[i][j] = (double)substitution_hosei( seq[i], seq[j] );
 /*
				 fprintf( stderr, "%f\n", eff[i][j] );
 */
   		 	}
/*
			fprintf( stderr, "distance matrix\n" );
			for( i=0; i<nseq; i++ )
			{
				for( j=0; j<nseq; j++ ) 
				{
					fprintf( stderr, "%f ", eff[i][j] );
				}
				fprintf( stderr, "\n" );
			}
*/
/*
   			upg( nseq, eff, topol, len );
   			upg2( nseq, eff, topol, len );
*/
   			spg( nseq, eff, topol, len );
   			counteff( nseq, topol, len, eff );
		}
	}
	else
	{
		for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) 
			eff[i][j] = 1.0;
	}
/*
fprintf( stderr, "weight matrix\n" );
for( i=0; i<nseq; i++ )
{
	for( j=0; j<nseq; j++ ) 
	{
		fprintf( stderr, "%f ", eff[i][j] );
	}
	fprintf( stderr, "\n" );
}
*/
}

float bscore_calc( char **seq, int s, double **eff )  /* algorithm B */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    long score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 + !gb1  * !gc1 
                 * !gb2  *  gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
				 + gb1   * !gc1
				 * gb2   *  gc2      *BEFF

				 + gb1   *  gc1
				 * gb2   * !gc2      *BEFF
                 ;
			score += (long)cob * penalty * efficient;
			score += (long)amino_dis[(int)seq[i][k]][(int)seq[j][k]] * efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (float)score / nglen + 400.0 * !scoremtx );
}

void AllocateTmpSeqs( char ***mseq2pt, char **mseq1pt, int locnlenmax )
{
	*mseq2pt = AllocateCharMtx( njob, locnlenmax+1 );
	*mseq1pt = AllocateCharVec( locnlenmax+1 );
}

void FreeTmpSeqs( char **mseq2, char *mseq1 )
{
	FreeCharMtx( mseq2 );
	free( (char *)mseq1 );
}

void gappick0( char *aseq, char *seq )
{


	for( ; *seq != 0; seq++ )
	{
		if( *seq != '-' )
			*aseq++ = *seq;
	}
	*aseq = 0;

}

void gappick( int nseq, int s, char **aseq, char **mseq2, 
			  double **eff, double *effarr )
{
	int i, j, count, countjob, len, allgap;
	len = strlen( aseq[0] );
	for( i=0, count=0; i<len; i++ ) 
	{
		allgap = 1;
		for( j=0; j<nseq; j++ ) if( j != s ) allgap *= ( aseq[j][i] == '-' );
        if( allgap == 0 )
		{
			for( j=0, countjob=0; j<nseq; j++ ) 
			{
				if( j != s )
				{
					mseq2[countjob][count] = aseq[j][i];
					countjob++;
				}
			}
			count++;
		}
	}
	for( i=0; i<nseq-1; i++ ) mseq2[i][count] = 0;

	for( i=0, countjob=0; i<nseq; i++ ) 
	{
		if( i != s )
		{
			effarr[countjob] = eff[s][i];
			countjob++;
		}
	}
/*
fprintf( stdout, "effarr in gappick s = %d\n", s+1 );
for( i=0; i<countjob; i++ ) 
	fprintf( stdout, " %f", effarr[i] );
printf( "\n" );
*/
}

void commongappick( int nseq, char **seq )
{
	int i, j, count;
	int len = strlen( seq[0] );

	for( i=0, count=0; i<=len; i++ ) 
	{
	/*
		allgap = 1;
		for( j=0; j<nseq; j++ ) 
			allgap *= ( seq[j][i] == '-' );
		if( !allgap )
	*/
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			for( j=0; j<nseq; j++ )
			{
				seq[j][count] = seq[j][i];
			}
			count++;
	 	}
	}
}
		
double score_calc0( char **seq, int s, double **eff, int ex )
{
	double tmp;

	if( scmtd == 3 ) tmp = score_calc3( seq, s, eff, ex );
	if( scmtd == 4 ) tmp = score_calc4( seq, s, eff, ex );
	if( scmtd == 5 ) tmp = score_calc5( seq, s, eff, ex );
	else             tmp = score_calc5( seq, s, eff, ex );

	return( tmp );

}

/*
float score_m_1( char **seq, int ex, double **eff )
{
	int i, j, k;
	int len = strlen( seq[0] );
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
	double score;

	score = 0.0;
	nglen = 0;
	for( i=0; i<njob; i++ ) 
	{
		double efficient = eff[MIN(i,ex)][MAX(i,ex)];
		if( i == ex ) continue;

		gc1 = 0; 
		gc2 = 0;
		for( k=0; k<len; k++ ) 
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[ex][k] == '-' );
      
            cob = 
                   !gb1  *  gc1
                 * !gb2  * !gc2

                 + !gb1  * !gc1
                 * !gb2  *  gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
                 +  gb1  * !gc1
                 *  gb2  *  gc2      *BEFF

                 +  gb1  *  gc1
                 *  gb2  * !gc2      *BEFF
                 ;
			score += (double)cob * penalty * efficient;
			score += (double)amino_dis[seq[i][k]][seq[ex][k]] * efficient;
			*
			nglen += ( !gc1 * !gc2 );
			*
			if( !gc1 && !gc2 ) fprintf( stdout, "%f\n", score );
		}
	}
	return( (float)score / nglen + 400.0 * !scoremtx );
}
*/

#if 0
void sitescore( char **seq, double **eff, char sco1[], char sco2[], char sco3[] )
{
	int i, j, k;
	int len = strlen( seq[0] );
	double tmp;
	double count;
	int ch;
	double sco[N];

	for( i=0; i<len; i++ ) 
	{
		tmp = 0.0; count = 0;
		for( j=0; j<njob-1; j++ ) for( k=j+1; k<njob; k++ ) 
		{
		/*
			if( seq[j][i] != '-' && seq[k][i] != '-' )
		*/
			{
				tmp += amino_dis[seq[j][i]][seq[k][i]] + 400 * !scoremtx;
				count++; 
			}
		}
		if( count > 0.0 ) tmp /= count;
		else( tmp = 0.0 );
		ch = (int)( tmp/100.0 - 0.000001 );
		sprintf( sco1+i, "%c", ch+0x61 );
	}
	sco1[len] = 0;

    for( i=0; i<len; i++ ) 
    {
        tmp = 0.0; count = 0;
        for( j=0; j<njob-1; j++ ) for( k=j+1; k<njob; k++ ) 
        {
		/*
            if( seq[j][i] != '-' && seq[k][i] != '-' )
		*/
            {
                tmp += eff[j][k] * ( amino_dis[seq[j][i]][seq[k][i]] + 400 * !scoremtx );
                count += eff[j][k]; 
            }
        }
		if( count > 0.0 ) tmp /= count;
		else( tmp = 0.0 );
		tmp = ( tmp - 400 * !scoremtx ) * 2;
		if( tmp < 0 ) tmp = 0;
        ch = (int)( tmp/100.0 - 0.000001 );
        sprintf( sco2+i, "%c", ch+0x61 );
		sco[i] = tmp;
    }
    sco2[len] = 0;

	for( i=WIN; i<len-WIN; i++ )
	{
		tmp = 0.0;
		for( j=i-WIN; j<=i+WIN; j++ )
		{
			tmp += sco[j];
		}
		for( j=0; j<njob; j++ ) 
		{
			if( seq[j][i] == '-' )
			{
				tmp = 0.0;
				break;
			}
		}
		tmp /= WIN * 2 + 1;
		ch = (int)( tmp/100.0 - 0.0000001 );
		sprintf( sco3+i, "%c", ch+0x61 );
	}
	for( i=0; i<WIN; i++ ) sco3[i] = '-';
	for( i=len-WIN; i<len; i++ ) sco3[i] = '-';
	sco3[len] = 0;
}
#endif

void strins( char *str1, char *str2 )
{
	char *bk;
	int len1 = strlen( str1 );
	int len2 = strlen( str2 );

	bk = str2;
	str2 += len1+len2;
	str1 += len1-1;

	while( str2 >= bk+len1 ) *str2-- = *(str2-len1);
	while( str2 >= bk ) *str2-- = *str1--;
}

int isaligned( int nseq, char **seq )
{
	int i;
	int len = strlen( seq[0] );
	for( i=1; i<nseq; i++ ) 
	{
		if( strlen( seq[i] ) != len ) return( 0 );
	}
	return( 1 );
}

double score_calc_for_score( int nseq, char **seq )
{
    int i, j, k, c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;

    score = 0.0;
    for( i=0; i<nseq-1; i++ )
    {
        for( j=i+1; j<nseq; j++ )
        {
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            c = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                c++;
                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty - n_dis[0][24];
                    while( mseq1[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty - n_dis[0][24];
                    while( mseq2[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore / (double)c;
#if DEBUG
			printf( "tmpscore in mltaln9.c = %f\n", tmpscore );
			printf( "tmpscore / c          = %f\n", tmpscore/(double)c );
#endif
        }
    }
	fprintf( stderr, "raw score = %f\n", score );
	score /= (double)nseq * ( nseq-1.0 ) / 2.0;
	score += 400.0;
#if DEBUG
	printf( "score in mltaln9.c = %f\n", score );
#endif
    return( (double)score );
}

void floatncpy( float *vec1, float *vec2, int len )
{
	while( len-- )
		*vec1++ = *vec2++;
}

float score_calc_a( char **seq, int s, double **eff )  /* algorithm A+ */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
                 ;
			score += 0.5 * (float)cob * penalty * efficient;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]] * (float)efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (float)score / nglen + 400.0 * !scoremtx );
}


float score_calc_s( char **seq, int s, double **eff )  /* algorithm S, not used */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
#if 0
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
#endif
                 ;
			score += 0.5 * (float)cob * penalty * efficient;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]] * (float)efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (float)score / nglen + 400.0 );
}

double score_calc_for_score_s( int s, char **seq )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
#if 0
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
#endif
                 ;
			score += 0.5 * (float)cob * penalty;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]];
			nglen += ( !gc1 * !gc2 );
		}
#if 0
		fprintf( stderr, "i = %d, j=%d\n", i+1, j+1 );
		fprintf( stderr, "score = %f\n", score );
#endif
	}
	return( (double)score / nglen + 400.0 );
}

double SSPscore___( int s, char **seq, int ex )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	i=ex; for( j=0; j<s; j++ )
	{

		if( j == ex ) continue;

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2 * 2.0 

                 +  gb1  * !gc1
                 * !gb2  *  gc2 * 2.0 
      
#if 0
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
#endif
                 ;
			score += 0.5 * (float)cob * penalty;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]];
			nglen += ( !gc1 * !gc2 ); /* tsukawanai */
		}
#if 0
		fprintf( stderr, "i = %d, j=%d\n", i+1, j+1 );
		fprintf( stderr, "score = %f\n", score );
#endif
	}
	return( (double)score );
}

double SSPscore( int s, char **seq )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
#if 0
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
#endif
                 ;
			score += 0.5 * (float)cob * penalty;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]];
			nglen += ( !gc1 * !gc2 ); /* tsukawanai */
		}
#if 0
		fprintf( stderr, "i = %d, j=%d\n", i+1, j+1 );
		fprintf( stderr, "score = %f\n", score );
#endif
	}
	return( (double)score );
}



double DSPscore( int s, char **seq )  /* method 3 deha nai */
{
    int i, j, k;
    double c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;
#if DEBUG
	FILE *fp;
#endif

    score = 0.0;
    c = 0.0;

    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore;
        }
    }

	return( score );
}


#define SEGMENTSIZE 150

int searchAnchors( int nseq, char **seq, Segment *seg )
{
	int i, j, k;
	int status;
	double score;
	int value = 0;
	int len;
	int length;
	static double *stra = NULL;
	static int alloclen = 0;
	double cumscore;
	static double threshold;

	len = strlen( seq[0] );
	if( alloclen < len )
	{
		if( alloclen )
		{
			FreeDoubleVec( stra );
		}
		else
		{
			threshold = (int)divThreshold / 100.0 * 600.0 * divWinSize;
		}
		stra = AllocateDoubleVec( len );
		alloclen = len;
	}

	for( i=0; i<len; i++ )
	{
#if 0
		/* make prf */
		for( j=0; j<26; j++ )
		{
			prf[j] = 0.0;
		}
		for( j=0; j<nseq; j++ ) prf[amino_n[seq[j][i]]] += 1.0;

		/* make hat */
		pre = 26;
		for( j=25; j>=0; j-- )
		{
			if( prf[j] )
			{
				hat[pre] = j;
				pre = j;
			}
		}
		hat[pre] = -1;

		/* make site score */
		stra[i] = 0.0;
		for( k=hat[26]; k!=-1; k=hat[k] ) 
			for( j=hat[26]; j!=-1; j=hat[j] ) 
				stra[i] += n_dis[k][j] * prf[k] * prf[j];
#else
		stra[i] = 0.0;
		for( k=0; k<nseq-1; k++ ) for( j=k+1; j<nseq; j++ )
			stra[i] += n_dis[(int)amino_n[(int)seq[k][i]]][(int)amino_n[(int)seq[j][i]]];
		stra[i] /= (double)nseq * ( nseq-1 ) / 2;
#endif
	}

	(seg+0)->skipForeward = 0;
	(seg+1)->skipBackward = 0;
	status = 0;
	cumscore = 0.0;
	score = 0.0;
	length = 0; /* modified at 01/09/11 */
	for( j=0; j<divWinSize; j++ ) score += stra[j];
	for( i=1; i<len-divWinSize; i++ )
	{
		score = score - stra[i-1] + stra[i+divWinSize-1];
#if DEBUG
		fprintf( stderr, "%d %f   ? %f", i, score, threshold );
		if( score > threshold ) fprintf( stderr, "YES\n" );
		else                    fprintf( stderr, "NO\n" );
#endif

		if( score > threshold )
		{
			if( !status )
			{
				status = 1;
				seg->start = i;
				length = 0;
				cumscore = 0.0;
			}
			length++;
			cumscore += score;
		}
		if( score <= threshold || length > SEGMENTSIZE )
		{
			if( status )
			{
				seg->end = i;
				seg->center = ( seg->start + seg->end + divWinSize ) / 2 ;
				seg->score = cumscore;
#if DEBUG
				fprintf( stderr, "%d-%d length = %d\n", seg->start, seg->end, length );
#endif
				if( length > SEGMENTSIZE )
				{
					(seg+0)->skipForeward = 1;
					(seg+1)->skipBackward = 1;
				}
				else
				{
					(seg+0)->skipForeward = 0;
					(seg+1)->skipBackward = 0;
				}
				length = 0;
				cumscore = 0.0;
				status = 0;
				value++;
				seg++;
				if( value > MAXSEG - 3 ) ErrorExit( "TOO MANY SEGMENTS!");
			}
		}
	}
	if( status )
	{
		seg->end = i;
		seg->center = ( seg->start + seg->end + divWinSize ) / 2 ;
		seg->score = cumscore;
#if DEBUG
fprintf( stderr, "%d-%d length = %d\n", seg->start, seg->end, length );
#endif
		value++;
	}
	return( value );
}

char *progName( char *str )
{
	char *value;
	if( ( value = strrchr( str, '/' ) ) != NULL )
		return( value+1 );
	else
		return( str );
}
