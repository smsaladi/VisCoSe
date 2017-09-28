#include "mltaln.h"

#define DEBUG 0

#if 0
void mdfymtx( char pair[njob][njob], int s1, double **partialmtx, double **mtx )
#else
void mdfymtx( char **pair, int s1, double **partialmtx, double **mtx )
#endif
{
	int i, j;
	int icount, jcount;
#if DEBUG
	FILE *fp;
	static char name[M][B];

	for( i=0; i<M; i++ ) name[i][0] = 0;
	fprintf( stdout, "s1 = %d\n", s1 );
	for( i=0; i<njob; i++ ) 
	{
		for( j=0; j<njob; j++ ) 
		{
			printf( "%#2d", pair[i][j] );
		}
		printf( "\n" );
	}
#endif

	for( i=0, icount=0; i<njob-1; i++ )
	{
		if( !pair[s1][i] ) continue;
		for( j=i+1, jcount=icount+1; j<njob; j++ ) 
		{
			if( !pair[s1][j] ) continue;
			partialmtx[icount][jcount] = mtx[i][j];
			jcount++;
		}
		icount++;
	}
#if DEBUG
	fp = fopen( "hat2.org", "w" );
	WriteHat2( fp, njob, name, mtx );
	fclose( fp );
	fp = fopen( "hat2.mdf", "w" );
	WriteHat2( fp, icount, name, partialmtx );
	fclose( fp );
#endif
		
}

		
float score_calc( char **seq, int s )  /* method 3  */
{
    int i, j, k, c;
    int len = strlen( seq[0] );
    float score;
    int tmpscore;
    char *mseq1, *mseq2;

    score = 0.0;
    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0;
            c = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                c++;
                tmpscore += amino_dis[mseq1[k]][mseq2[k]];
                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
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
                        tmpscore -= penalty;
                        break;
                    }
                    else break;
                }
            }
            if( mseq1[len-1] == '-' || mseq2[len-1] == '-' )
            {
                for( k=0; k<len; k++ )
                {
                    if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                    if( !( mseq1[k] != '-' && mseq2[k] != '-' ) ) 
                    {
                        c--;
                        tmpscore -= penalty;
                        break;
                    }
                    else break;
                }
            }
            */
            score += (double)tmpscore / (double)c;
        }
    }
    score = (float)score / ( ( (double)s * ((double)s-1.0) ) / 2.0 );
	fprintf( stderr, "score in score_calc = %f\n", score );
    return( score );
}

void cpmx_calc( char **seq, float **cpmx, double *eff, int lgth, int clus )
{
	int  i, j, k;
	double totaleff = 0.0;

	for( i=0; i<clus; i++ ) totaleff += eff[i]; 
	for( i=0; i<26; i++ ) for( j=0; j<lgth; j++ ) cpmx[i][j] = 0.0;
	for( j=0; j<lgth; j++ ) for( k=0; k<clus; k++ )
			cpmx[amino_n[seq[k][j]]][j] += (float)eff[k] / totaleff;
}

void mseqcat( char **seq1, char **seq2, double **eff, double *effarr1, double *effarr2, char name1[M][B], char name2[M][B], int clus1, int clus2 )
{
	int i, j;
    for( i=0; i<clus2; i++ )
        seq1[clus1+i] = seq2[i];
    for( i=0; i<clus2; i++ ) strcpy( name1[clus1+i], name2[i] );

	for( i=0; i<clus1; i++ ) for( j=0; j<clus1; j++ ) eff[i][j] = effarr1[i]* effarr1[j]; 
	for( i=0; i<clus1; i++ ) for( j=clus1; j<clus1+clus2; j++ ) eff[i][j] = effarr1[i]* effarr2[j-clus1]; 
	for( i=clus1; i<clus1+clus2; i++ ) for( j=0; j<clus1; j++ ) eff[i][j] = effarr2[i-clus1] * effarr1[j]; 
	for( i=clus1; i<clus1+clus2; i++ ) for( j=clus1; j<clus1+clus2; j++ ) eff[i][j] = effarr2[i-clus1] * effarr2[j-clus1]; 
}

	

void strnbcat( char *s1, char *s2, int m )   /* s1 + s2 ---> s2  */
{
	static char b[N];
	strncpy( b, s1, m ); b[m] = 0;
	strcat( b, s2 );
	strcpy( s2, b );
}

#if 0
int conjuction( char pair[njob][njob], int s, char **seq, char **aseq, double *peff, double *eff, char name[M][B], char aname[M][B], char *d )
#else
int conjuctionforgaln( int s0, int s1, char **seq, char **aseq, double *peff, double *eff, char **name, char **aname, char *d )
#endif
{
	int m, k;
	char b[B];

#if 1
	fprintf( stderr, "s0 = %d, s1 = %d\n", s0, s1 );
#endif

	d[0] = 0;
	for( m=s0, k=0; m<s1; m++ )
	{
		{
			sprintf( b, " %d", m+1 ); 
#if 1
			if( strlen( d ) < 100 ) strcat( d, b );
#else
			strcat( d, b );
#endif
			aseq[k] = seq[m];
			peff[k] = eff[m];
#if 0
			strcpy( aname[k], name[m] );
#endif
			k++;
		}
	}
	return( k );
}

int conjuction( char **pair, int s, char **seq, char **aseq, double *peff, double *eff, char name[M][B], char aname[M][B], char *d )
{
	int m, k;
	char b[B];

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	d[0] = 0;
	for( m=s, k=0; m<njob; m++ )
	{
		if( pair[s][m] != 0 )
		{
			sprintf( b, " %d", m+1 ); 
#if 1
			if( strlen( d ) < 100 ) strcat( d, b );
#else
			strcat( d, b );
#endif
			aseq[k] = seq[m];
			peff[k] = eff[m];
#if 0
			strcpy( aname[k], name[m] );
#endif
			k++;
		}
	}
	return( k );
}
			
void floatdelete( float **cpmx, int d, int len )
{
	int i, j;

	for( i=d; i<len-1; i++ )
	{
		for( j=0; j<26; j++ )
		{
			cpmx[j][i] = cpmx[j][i+1]; 
		}
	}
}
		
void chardelete( char *seq, int d )
{
	char b[N]; 

		strcpy( b, seq+d+1 );
		strcpy( seq+d, b );
}

int RootBranchNode( int nseq, int ***topol, int step, int branch )
{
	int i, j, s1, s2, value;

	value = 1;
	for( i=step+1; i<nseq-2; i++ ) 
	{
		for( j=0; (s1=topol[i][0][j])>-1; j++ ) 
			if( s1 == topol[step][branch][0] ) value++;
		for( j=0; (s2=topol[i][1][j])>-1; j++ ) 
			if( s2 == topol[step][branch][0] ) value++;
	}
	return( value );
}

void BranchLeafNode( int nseq, int ***topol, int *node, int step, int branch )
{
	int i, j, k, s;

	for( i=0; i<nseq; i++ ) node[i] = 0;
	for( i=0; i<step-1; i++ )
		for( k=0; k<2; k++ ) 
			for( j=0; (s=topol[i][k][j])>-1; j++ ) 
				node[s]++;
	for( k=0; k<branch+1; k++ ) 
		for( j=0; (s=topol[step][k][j])>-1; j++ )
			node[s]++;
}

void RootLeafNode( int nseq, int ***topol, int *node )
{
	int i, j, k, s;
	for( i=0; i<nseq; i++ ) node[i] = 0;
	for( i=0; i<nseq-2; i++ )
		for( k=0; k<2; k++ ) 
			for( j=0; (s=topol[i][k][j])>-1; j++ ) 
				node[s]++;
}
		
void nodeFromABranch( int nseq, int *result, int **pairwisenode, int ***topol, double **len, int step, int num )
{
	int i, s, count;
	int *innergroup;
	int *outergroup1;
#if 0
	int outergroup2[nseq];
	int table[nseq];
#else
	static int *outergroup2 = NULL;
	static int *table = NULL;
	if( outergroup2 == NULL )
	{
		outergroup2 = AllocateIntVec( nseq );
		table = AllocateIntVec( nseq );
	}
#endif
	innergroup = topol[step][num];
	outergroup1 = topol[step][!num];
	for( i=0; i<nseq; i++ ) table[i] = 1;
	for( i=0; (s=innergroup[i])>-1; i++ ) table[s] = 0;
	for( i=0; (s=outergroup1[i])>-1; i++ ) table[s] = 0;
	for( i=0, count=0; i<nseq; i++ ) 
	{
		if( table[i] )
		{
			outergroup2[count] = i;
			count++;
		}
	}
	outergroup2[count] = -1;

	for( i=0; (s=innergroup[i])>-1; i++ )
	{
		result[s] = pairwisenode[s][outergroup1[0]]
				  + pairwisenode[s][outergroup2[0]]
				  - pairwisenode[outergroup1[0]][outergroup2[0]] - 1;
		result[s] /= 2;
	}
	for( i=0; (s=outergroup1[i])>-1; i++ ) 
	{
		result[s] = pairwisenode[s][outergroup2[0]]
				  + pairwisenode[s][innergroup[0]]
				  - pairwisenode[innergroup[0]][outergroup2[0]] + 1;
		result[s] /= 2;
	}
	for( i=0; (s=outergroup2[i])>-1; i++ ) 
	{
		result[s] = pairwisenode[s][outergroup1[0]]
				  + pairwisenode[s][innergroup[0]]
				  - pairwisenode[innergroup[0]][outergroup1[0]] + 1;
		result[s] /= 2;
	}

#if 0
	for( i=0; i<nseq; i++ ) 
		fprintf( stderr, "result[%d] = %d\n", i+1, result[i] );
#endif
}
		




	

		
		
					
					

				
		

#if 0
void OneClusterAndTheOther( int locnjob, char pair[njob][njob], int *s1, int *s2, int ***topol, int step, int branch )
#else
void OneClusterAndTheOther( int locnjob, char **pair, int *s1, int *s2, int ***topol, int step, int branch )
#endif
{
	int i;
	int r1;
	
    *s1 = topol[step][branch][0];
    for( i=0; (r1=topol[step][branch][i])>-1; i++ ) 
        pair[*s1][r1] = 1;
    for( i=0; i<locnjob; i++ ) 
    {
        if( !pair[*s1][i] ) 
        {
            *s2 = i;
            break;
        }
    }
    for( i=*s2; i<locnjob; i++ ) 
    {
        if( !pair[*s1][i] )
            pair[*s2][i] = 1;
    }
}

void makeEffMtx( int nseq, double **mtx, double *vec )
{
	int i, j;
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) 
		mtx[i][j] = vec[i] * vec[j];
}
	
void node_eff( int nseq, double *eff, int *node )
{
	int i;
	extern double ipower( double, int );
	for( i=0; i<nseq; i++ ) 
		eff[i] = ipower( 0.5, node[i] ) + geta2;
	/*
		eff[i] = ipower( 0.5, node[i] ) + 0;
	*/
#if DEBUG
	for( i=0; i<nseq; i++ ) 
#endif
}
