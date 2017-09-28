#include "fft.h"
#include "mltaln.h"

#define SEGMENTSIZE 150
#define TMPTMPTMP 0

#define DEBUG 0

void keika( char *str, int current, int all )
{
	if( current == 0 )
		fprintf( stderr, "%s :         ", str );

		fprintf( stderr, "\b\b\b\b\b\b\b\b" );
		fprintf( stderr, "%#3d /%#3d", current+1, all+1 );

	if( current+1 == all )
		fprintf( stderr, "\b\b\b\b\b\b\b\bdone.     \n" );
}

double maxItch( double *soukan, int size )
{
	int i;
	double value = 0.0;
	double cand;
	for( i=0; i<size; i++ ) 
		if( ( cand = soukan[i] ) > value ) value = cand;
	return( value );
}

void calcNaiseki( Fukusosuu *value, Fukusosuu *x, Fukusosuu *y )
{
	value->R =  x->R * y->R + x->I * y->I;
	value->I = -x->R * y->I + x->I * y->R;
}

Fukusosuu *AllocateFukusosuuVec( int l1 )
{
	Fukusosuu *value;
	value = (Fukusosuu *)calloc( l1, sizeof( Fukusosuu ) );
	if( !value ) ErrorExit( "Cannot allocate FukusosuuVec" );
	return( value );
}
	
Fukusosuu **AllocateFukusosuuMtx( int l1, int l2 )
{
	Fukusosuu **value;
	int j;
	value = (Fukusosuu **)calloc( l1+1, sizeof( Fukusosuu * ) );
	if( !value ) ErrorExit( "Cannot allocate FukusosuuVecMtx" );
	for( j=0; j<l1; j++ ) 
		value[j] = AllocateFukusosuuVec( l2 );
	value[l1] = NULL;
	return( value );
}

Fukusosuu ***AllocateFukusosuuCub( int l1, int l2, int l3 )
{
	Fukusosuu ***value;
	int i;
	value = calloc( l1+1, sizeof( Fukusosuu ** ) );
	if( !value ) ErrorExit( "Cannot allocate Fukusosuu" );
	for( i=0; i<l1; i++ ) value[i] = AllocateFukusosuuMtx( l2, l3 );
	value[l1] = NULL;
	return( value );
}

void FreeFukusosuuVec( Fukusosuu *vec )
{
	free( (void *)vec );
}

void FreeFukusosuuMtx( Fukusosuu **mtx )
{
	int i;

	for( i=0; mtx[i]; i++ ) 
		free( (void *)mtx[i] );
	free( (void *)mtx );
}

void getKouho( int *kouho, int nkouho, double *soukan, int nlen2 )
{
	int i, j;
	int nlen4 = nlen2 / 2;
	double max;
	double tmp;
	int ikouho;
	for( j=0; j<nkouho; j++ ) /* $B>/$7L5BL(B */
	{
		max = -9999.9;
		for( i=0; i<nlen2; i++ ) 
		{
			if( ( tmp = soukan[i] ) > max )
			{
				ikouho = i;
				max = tmp;
#if DEBUG
				fprintf( stderr, "Kouho No.%d, pos=%d, score=%f\n", j, i, tmp );
#endif
			}
		}
		soukan[ikouho] = -9999.9;
		kouho[j] = ( ikouho - nlen4 );
	}
}

static void gapfill( char *seq, int gaplen )
{
	while( gaplen-- )
		*seq++ = '-';
	*seq = 0;
}
	
void zurasu2( int lag, int    clus1, int    clus2, 
                       char  **seq1, char  **seq2, 
		 			   char **aseq1, char **aseq2 )
{
	int i;
#if DEBUG
	fprintf( stderr, "lag = %d\n", lag );
#endif
	if( lag > 0 )
	{
		for( i=0; i<clus1; i++ ) aseq1[i] = seq1[i];
		for( i=0; i<clus2; i++ ) aseq2[i] = seq2[i]+lag;
	}
	else
	{
		for( i=0; i<clus1; i++ ) aseq1[i] = seq1[i]-lag;
		for( i=0; i<clus2; i++ ) aseq2[i] = seq2[i];
	}
}

void zurasu( int lag, int    clus1, int    clus2, 
                      char  **seq1, char  **seq2, 
					  char **aseq1, char **aseq2 )
{
	int i;
#if DEBUG
	fprintf( stderr, "lag = %d\n", lag );
#endif
	if( lag > 0 )
	{
		for( i=0; i<clus1; i++ ) strcpy( aseq1[i], seq1[i] );
		for( i=0; i<clus2; i++ ) strcpy( aseq2[i], seq2[i]+lag );
	}
	else
	{
		for( i=0; i<clus1; i++ ) strcpy( aseq1[i], seq1[i]-lag );
		for( i=0; i<clus2; i++ ) strcpy( aseq2[i], seq2[i] );
	}
}


int alignableReagion( int    clus1, int    clus2, 
					  char  **seq1, char  **seq2,
					  double *eff1, double *eff2,
					  Segment *seg )
{
	int i, j, k, l;
	int status;
	double score;
	int value = 0;
	int len, maxlen;
	int length;
	static double *stra = NULL;
	static int alloclen = 0;
	double totaleff;
	double cumscore;
	static double threshold;
	static double prf1[26], prf2[26];
	static int hat1[27], hat2[27];
	int pre1, pre2;

	len = MIN( strlen( seq1[0] ), strlen( seq2[0] ) );
	maxlen = MAX( strlen( seq1[0] ), strlen( seq2[0] ) ) + fftWinSize;
	if( alloclen < maxlen )
	{
		if( alloclen )
		{
			FreeDoubleVec( stra );
		}
		else
		{
			threshold = (int)fftThreshold / 100.0 * 600.0 * fftWinSize;
		}
		stra = AllocateDoubleVec( maxlen );
		alloclen = maxlen;
	}

	totaleff = 0.0;
	for( i=0; i<clus1; i++ ) for( j=0; j<clus2; j++ ) totaleff += eff1[i] * eff2[j];
	for( i=0; i<len; i++ )
	{
		/* make prfs */
		for( j=0; j<26; j++ )
		{
			prf1[j] = 0.0;
			prf2[j] = 0.0;
		}
		for( j=0; j<clus1; j++ ) prf1[amino_n[seq1[j][i]]] += eff1[j];
		for( j=0; j<clus2; j++ ) prf2[amino_n[seq2[j][i]]] += eff2[j];

		/* make hats */
		pre1 = pre2 = 26;
		for( j=25; j>=0; j-- )
		{
			if( prf1[j] )
			{
				hat1[pre1] = j;
				pre1 = j;
			}
			if( prf2[j] )
			{
				hat2[pre2] = j;
				pre2 = j;
			}
		}
		hat1[pre1] = -1;
		hat2[pre2] = -1;

		/* make site score */
		stra[i] = 0.0;
		for( k=hat1[26]; k!=-1; k=hat1[k] ) 
			for( j=hat2[26]; j!=-1; j=hat2[j] ) 
				stra[i] += n_dis[k][j] * prf1[k] * prf2[j];
		stra[i] /= totaleff;
	}

	(seg+0)->skipForeward = 0;
	(seg+1)->skipBackward = 0;
	status = 0;
	cumscore = 0.0;
	score = 0.0;
	for( j=0; j<fftWinSize; j++ ) score += stra[j];

	for( i=1; i<len-fftWinSize; i++ )
	{
		score = score - stra[i-1] + stra[i+fftWinSize-1];
#if TMPTMPTMP
		fprintf( stdout, "%d %10.0f   ? %10.0f\n", i, score, threshold );
#endif

		if( score > threshold )
		{
#if 0
				seg->start = i;
				seg->end = i;
				seg->center = ( seg->start + seg->end + fftWinSize ) / 2 ;
				seg->score = score;
				status = 0;
				value++;
#else
			if( !status )
			{
				status = 1;
				seg->start = i;
				length = 0;
				cumscore = 0.0;
			}
			length++;
			cumscore += score;
#endif
		}
		if( score <= threshold || length > SEGMENTSIZE )
		{
			if( status )
			{
				seg->end = i;
				seg->center = ( seg->start + seg->end + fftWinSize ) / 2 ;
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
		seg->center = ( seg->start + seg->end + fftWinSize ) / 2 ;
		seg->score = cumscore;
#if DEBUG
fprintf( stderr, "%d-%d length = %d\n", seg->start, seg->end, length );
#endif
		value++;
	}
#if TMPTMPTMP
	exit( 0 );
#endif
	return( value );
}

void blockAlign( int *cut1, int *cut2, double **ocrossscore, int *ncut )
{
	int i, j, shift, cur1, cur2, count;
	static int result1[MAXSEG], result2[MAXSEG];
	static int ocut1[MAXSEG], ocut2[MAXSEG];
	double maximum;
	static double max[MAXSEG], maxi;
	static double point[MAXSEG], pointi;
	static double **crossscore = NULL;
	static int **track = NULL;

	if( crossscore == NULL )
	{
		crossscore = AllocateDoubleMtx( MAXSEG, MAXSEG );
		track = AllocateIntMtx( MAXSEG, MAXSEG );
	}

#if DEBUG
	for( i=0; i<*ncut; i++ )
		fprintf( stderr, "i=%d, cut1 = %d, cut2 = %d\n", i, cut1[i], cut2[i] );
	for( i=0; i<*ncut; i++ ) 
	{
		for( j=0; j<*ncut; j++ )
			fprintf( stderr, "%#4.0f ", ocrossscore[i][j]/100 );
		fprintf( stderr, "\n" );
	}
#endif

	for( i=0; i<*ncut; i++ ) for( j=0; j<*ncut; j++ )  /* mudadanaa */
		crossscore[i][j] = ocrossscore[i][j];
	for( i=0; i<*ncut; i++ ) 
	{
		ocut1[i] = cut1[i];
		ocut2[i] = cut2[i];
	}

	for( j=0; j<*ncut; j++ )
	{
		max[j] = 0.0;
		point[j] = 0;
	}

	for( i=1; i<*ncut; i++ )
	{
		maxi = 0.0;
		pointi = 0;
		for( j=1; j<*ncut; j++ )
		{
			maximum = crossscore[i-1][j-1];
			track[i][j] = 0;

			if( maximum < max[j] + penalty )
			{
				maximum = max[j] + penalty;
				track[i][j] = point[j] - i;
			}
			if( maximum < maxi + penalty )
			{
				maximum = maxi + penalty;
				track[i][j] = j - pointi;
			}
			crossscore[i][j] += maximum;

			if( maxi < crossscore[i-1][j-1] )
			{
				maxi = crossscore[i-1][j-1];
				pointi = j-1;
			}
			if( max[j] < crossscore[i-1][j-1] )
			{
				max[j] = crossscore[i-1][j-1];
				point[j] = i-1;
			}
		}
	}
#if 0
	for( i=0; i<*ncut; i++ ) 
	{
		for( j=0; j<*ncut; j++ )
			fprintf( stderr, "%3d ", track[i][j] );
		fprintf( stderr, "\n" );
	}
#endif


	result1[MAXSEG-1] = *ncut-1;
	result2[MAXSEG-1] = *ncut-1;

	for( i=MAXSEG-1; i>=1; i-- )
	{
		cur1 = result1[i];
		cur2 = result2[i];
		if( cur1 == 0 || cur2 == 0 ) break;
		shift = track[cur1][cur2];
		if( shift == 0 )
		{
			result1[i-1] = cur1 - 1;
			result2[i-1] = cur2 - 1;
			continue;
		}
		else if( shift > 0 )
		{
			result1[i-1] = cur1 - 1;
			result2[i-1] = cur2 - shift;
		}
		else if( shift < 0 )
		{
			result1[i-1] = cur1 + shift;
			result2[i-1] = cur2 - 1;
		}
	}

	count = 0;
	for( j=i; j<MAXSEG; j++ )
	{
		if( ocrossscore[result1[j]][result2[j]] == 0.0 ) continue;

		if( result1[j] == result1[j-1] || result2[j] == result2[j-1] )
			if( ocrossscore[result1[j]][result2[j]] > ocrossscore[result1[j-1]][result2[j-1]] )
				count--;
				
		cut1[count] = ocut1[result1[j]];
		cut2[count] = ocut2[result2[j]];
		count++;
	}

	*ncut = count;
#if 0
	for( i=0; i<*ncut; i++ )
		fprintf( stderr, "i=%d, cut1 = %d, cut2 = %d\n", i, cut1[i], cut2[i] );
#endif
}

static int permit( Segment *seg1, Segment *seg2 )
{
	return( 0 );
	if( seg1->end >= seg2->start ) return( 0 );
	if( seg1->pair->end >= seg2->pair->start ) return( 0 );
	else return( 1 );
}

void blockAlign2( int *cut1, int *cut2, Segment **seg1, Segment **seg2, double **ocrossscore, int *ncut )
{
	int i, j, k, l, shift, cur1, cur2, count;
	static int crossscoresize = 0;
	static int result1[MAXSEG], result2[MAXSEG];
	static int ocut1[MAXSEG], ocut2[MAXSEG];
	double maximum;
	static double **crossscore = NULL;
	static int **track = NULL;
	static double maxj, maxi;
	static int pointj, pointi;

    if( crossscoresize < *ncut+2 )
    {
        crossscoresize = *ncut+2;
		fprintf( stderr, "allocating crossscore and track, size = %d\n", crossscoresize );
		if( track ) FreeIntMtx( track );
        if( crossscore ) FreeDoubleMtx( crossscore );
		track = AllocateIntMtx( crossscoresize, crossscoresize );
        crossscore = AllocateDoubleMtx( crossscoresize, crossscoresize );
    }

#if DEBUG
	for( i=0; i<*ncut-2; i++ )
		fprintf( stderr, "%d.start = %d, score = %f\n", i, seg1[i]->start, seg1[i]->score );

	for( i=0; i<*ncut; i++ )
		fprintf( stderr, "i=%d, cut1 = %d, cut2 = %d\n", i, cut1[i], cut2[i] );
	for( i=0; i<*ncut; i++ ) 
	{
		for( j=0; j<*ncut; j++ )
			fprintf( stderr, "%#4.0f ", ocrossscore[i][j]/100 );
		fprintf( stderr, "\n" );
	}
#endif

	for( i=0; i<*ncut; i++ ) for( j=0; j<*ncut; j++ )  /* mudadanaa */
		crossscore[i][j] = ocrossscore[i][j];
	for( i=0; i<*ncut; i++ ) 
	{
		ocut1[i] = cut1[i];
		ocut2[i] = cut2[i];
	}

	for( i=1; i<*ncut; i++ )
	{
		for( j=1; j<*ncut; j++ )
		{
#if 0
			fprintf( stderr, "i=%d, j=%d\n", i, j );
#endif
			pointi = 0; maxi = 0.0;
			for( k=0; k<j-2; k++ )
			{
/*
				fprintf( stderr, "k=%d, i=%d\n", k, i );
*/
				if( k && k<*ncut-1 && j<*ncut-1 && !permit( seg1[k-1], seg1[j-1] ) ) continue;
				if( crossscore[i-1][k] > maxj )
				{
					pointi = k;
					maxi = crossscore[i-1][k];
				}
			}

			pointj = 0; maxj = 0.0;
			for( k=0; k<i-2; k++ )
			{
				if( k && k<*ncut-1 && i<*ncut-1 && !permit( seg2[k-1], seg2[i-1] ) ) continue;
				if( crossscore[k][j-1] > maxj )
				{
					pointj = k;
					maxj = crossscore[k][j-1];
				}
			}	

			maxi += penalty;
			maxj += penalty;

			maximum = crossscore[i-1][j-1];
			track[i][j] = 0;

			if( maximum < maxi )
			{
				maximum = maxi ;
				track[i][j] = j - pointi;
			}

			if( maximum < maxj )
			{
				maximum = maxj ;
				track[i][j] = pointj - i;
			}

			crossscore[i][j] += maximum;
		}
	}
#if 0
	for( i=0; i<*ncut; i++ ) 
	{
		for( j=0; j<*ncut; j++ )
			fprintf( stderr, "%3d ", track[i][j] );
		fprintf( stderr, "\n" );
	}
#endif


	result1[MAXSEG-1] = *ncut-1;
	result2[MAXSEG-1] = *ncut-1;

	for( i=MAXSEG-1; i>=1; i-- )
	{
		cur1 = result1[i];
		cur2 = result2[i];
		if( cur1 == 0 || cur2 == 0 ) break;
		shift = track[cur1][cur2];
		if( shift == 0 )
		{
			result1[i-1] = cur1 - 1;
			result2[i-1] = cur2 - 1;
			continue;
		}
		else if( shift > 0 )
		{
			result1[i-1] = cur1 - 1;
			result2[i-1] = cur2 - shift;
		}
		else if( shift < 0 )
		{
			result1[i-1] = cur1 + shift;
			result2[i-1] = cur2 - 1;
		}
	}

	count = 0;
	for( j=i; j<MAXSEG; j++ )
	{
		if( ocrossscore[result1[j]][result2[j]] == 0.0 ) continue;

		if( result1[j] == result1[j-1] || result2[j] == result2[j-1] )
			if( ocrossscore[result1[j]][result2[j]] > ocrossscore[result1[j-1]][result2[j-1]] )
				count--;
				
		cut1[count] = ocut1[result1[j]];
		cut2[count] = ocut2[result2[j]];
		count++;
	}

	*ncut = count;
#if 0
	for( i=0; i<*ncut; i++ )
		fprintf( stderr, "i=%d, cut1 = %d, cut2 = %d\n", i, cut1[i], cut2[i] );
#endif
}

