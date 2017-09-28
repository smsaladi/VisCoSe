#include "mltaln.h"
#include "fft.h"

static FILE *fftfp;
static int n20or4or2;

#define KEIKA 0
#define RND   0
#define DEBUG 0

extern int fft( int, Fukusosuu *, int );


static void generateRndSeq( char *seq, int len )
{
	while( len-- )
#if 1
		*seq++ = (int)( rnd() * n20or4or2 );
#else
		*seq++ = (int)1;
#endif
}

static void vec_init( Fukusosuu *result, int nlen )
{
	while( nlen-- ) result->R = result++->I = 0.0;
}

static void vec_init2( Fukusosuu **result, char *seq, double eff, int st, int ed )
{
	int i;
	int n;
	for( i=st; i<ed; i++ )
		result[*seq++][i].R += eff;
}

static void seq_vec_2( Fukusosuu *result, double *score, double incr, char *seq )
{
	static int n;
	for( ; *seq; result++ )
	{
		n = amino_n[(int)*seq++];
		if( n < 20 && n >= 0 ) result->R += incr * score[n];
#if 0
		fprintf( stderr, "n=%d, score=%f, inc=%f R=%f\n",n,  score[n], incr * score[n], result->R );
#endif
	}
}

static void seq_vec_3( Fukusosuu **result, double incr, char *seq )
{
	int i;
	int n;
	for( i=0; *seq; i++ )
	{
		n = amino_n[(int)*seq++];
		if( n < n20or4or2 && n >= 0 ) result[n][i].R += incr;
	}
}

	
static void seq_vec( Fukusosuu *result, char query, double incr, char *seq )
{
#if 0
	int bk = nlen;
#endif
	while( *seq )
	{
		if( *seq++ == query ) result->R += incr;
		result++;
#if 0
fprintf( stderr, "i = %d result->R = %f\n", bk-nlen, (result-1)->R );
#endif
	}
}

int checkRepeat( int num, int *cutpos )
{
	int tmp, buf;

	buf = *cutpos;
	while( num-- )
	{
		if( ( tmp = *cutpos++ ) < buf ) return( 1 );
		buf = tmp;
	}
	return( 0 );
}

int segcmp( void *ptr1, void *ptr2 )
{
	int diff;
	Segment **seg1 = (Segment **)ptr1;
	Segment **seg2 = (Segment **)ptr2;
#if 0
	return( (*seg1)->center - (*seg2)->center );
#else
	diff = (*seg1)->center - (*seg2)->center;
	if( diff ) return( diff );

	diff = (*seg1)->start - (*seg2)->start;
	if( diff ) return( diff );

	diff = (*seg1)->end - (*seg2)->end;
	if( diff ) return( diff );

	fprintf( stderr, "USE STABLE SORT !!\n" );
	exit( 1 );
	return( 0 );
#endif
}


void mymergesort( int first, int last, Segment **seg )
{
	int middle;
	static int i, j, k, p;
	static int allo = 0;
	static Segment **work = NULL;
	if( last > allo )
	{
		allo = last;
		work = (Segment **)calloc( allo / 2 + 1, sizeof( Segment *) );
	}

	if( first < last )
	{
		middle = ( first + last ) / 2;
		mymergesort( first, middle, seg );
		mymergesort( middle+1, last, seg );
		p = 0;
		for( i=first; i<=middle; i++ ) work[p++] = seg[i];
		i = middle + 1; j = 0; k = first;
		while( i <= last && j < p )
		{
			if( work[j]->center <= seg[i]->center ) 
				seg[k++] = work[j++];
			else
				seg[k++] = seg[i++];
		}
		while( j < p ) seg[k++] = work[j++];
	}
}

void tmpsort( int num, Segment **seg )
{
	Segment *buf;
	int i, j;
	for( i=0; i<num-1; i++ )
	{
		for( j=i+1; j<num; j++ )
		{
			if( seg[j]->center < seg[i]->center ) 
			{
				buf = seg[i];
				seg[i] = seg[j];
				seg[j] = buf;
			}
		}
#if D0EBUG
		fprintf( stderr, "[%d], %d, %d\n", i, seg[i]->center, seg[i]->center );
#endif
	}
}

double Fgetlag( char  **seq1, char  **seq2, 
			    double *eff1, double *eff2, 
			    int    clus1, int    clus2,
			    int alloclen )
{
	int i, j, k, l, m;
	int nlen, nlen2, nlen4;
	static int crossscoresize = 0;
	static char **tmpseq1 = NULL;
	static char **tmpseq2 = NULL;
	static char **tmpptr1 = NULL;
	static char **tmpptr2 = NULL;
	static char **tmpres1 = NULL;
	static char **tmpres2 = NULL;
	static char **result1 = NULL;
	static char **result2 = NULL;
#if RND
	static char **rndseq1 = NULL;
	static char **rndseq2 = NULL;
#endif
	static Fukusosuu **seqVector1 = NULL;
	static Fukusosuu **seqVector2 = NULL;
	static Fukusosuu **naiseki = NULL;   
	static Fukusosuu *naisekiNoWa = NULL; 
	static double *soukan = NULL;
	static double **crossscore = NULL;
	int nlentmp;
	static int *kouho = NULL;
	static Segment *segment = NULL;
	static Segment *segment1 = NULL;
	static Segment *segment2 = NULL;
	static Segment **sortedseg1 = NULL;
	static Segment **sortedseg2 = NULL;
	static int *cut1 = NULL;
	static int *cut2 = NULL;
	static int localalloclen = 0;
	int lag;
	int tmpint;
	int count, count0;
	int len1, len2;
	int totallen;
	double totaleff1, totaleff2;

	extern Fukusosuu   *AllocateFukusosuuVec();
	extern Fukusosuu  **AllocateFukusosuuMtx();

	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

	len1 = strlen( seq1[0] );
	len2 = strlen( seq2[0] );
	nlentmp = MAX( len1, len2 );

	nlen = 1;
	while( nlentmp >= nlen ) nlen <<= 1;
#if 0
	fprintf( stderr, "###   nlen    = %d\n", nlen );
#endif

	nlen2 = nlen/2; nlen4 = nlen2 / 2;

#if DEBUG
	fprintf( stderr, "len1 = %d, len2 = %d\n", len1, len2 );
	fprintf( stderr, "nlentmp = %d, nlen = %d\n", nlentmp, nlen );
#endif

	if( !localalloclen )
	{
		kouho = AllocateIntVec( NKOUHO );
		cut1 = AllocateIntVec( MAXSEG );
		cut2 = AllocateIntVec( MAXSEG );
		tmpptr1 = AllocateCharMtx( njob, 1 );
		tmpptr2 = AllocateCharMtx( njob, 1 );
		result1 = AllocateCharMtx( njob, alloclen );
		result2 = AllocateCharMtx( njob, alloclen );
		tmpres1 = AllocateCharMtx( njob, alloclen );
		tmpres2 = AllocateCharMtx( njob, alloclen );
//		crossscore = AllocateDoubleMtx( MAXSEG, MAXSEG );
		segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment1 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment2 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		sortedseg1 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		sortedseg2 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		if( !( segment && segment1 && segment2 && sortedseg1 && sortedseg2 ) )
			ErrorExit( "Allocation error\n" );

		if     ( scoremtx == -1 ) n20or4or2 = 4;
		else if( fftscore == 1  ) n20or4or2 = 2;
		else                      n20or4or2 = 20;
	}
	if( localalloclen < nlen )
	{
		if( localalloclen )
		{
#if 1
			FreeFukusosuuMtx ( seqVector1 );
			FreeFukusosuuMtx ( seqVector2 );
			FreeFukusosuuVec( naisekiNoWa );
			FreeFukusosuuVec( naiseki );
			FreeDoubleVec( soukan );
			FreeCharMtx( tmpseq1 );
			FreeCharMtx( tmpseq2 );
#endif
#if RND
			FreeCharMtx( rndseq1 );
			FreeCharMtx( rndseq2 );
#endif
		}

		tmpseq1 = AllocateCharMtx( njob, nlen );
		tmpseq2 = AllocateCharMtx( njob, nlen );
		naisekiNoWa = AllocateFukusosuuVec( nlen );
		naiseki = AllocateFukusosuuMtx( n20or4or2, nlen );
		seqVector1 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
		seqVector2 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
		soukan = AllocateDoubleVec( nlen+1 );
#if RND
		rndseq1 = AllocateCharMtx( njob, nlen );
		rndseq2 = AllocateCharMtx( njob, nlen );
		for( i=0; i<njob; i++ )
		{
			generateRndSeq( rndseq1[i], nlen );
			generateRndSeq( rndseq2[i], nlen );
		}
#endif
		localalloclen = nlen;
	}
	
	for( j=0; j<clus1; j++ ) strcpy( tmpseq1[j], seq1[j] );
	for( j=0; j<clus2; j++ ) strcpy( tmpseq2[j], seq2[j] );

#if 0
fftfp = fopen( "input_of_Falign", "w" );
fprintf( fftfp, "nlen = %d\n", nlen );
fprintf( fftfp, "seq1: ( %d sequences ) \n", clus1 );
for( i=0; i<clus1; i++ )
	fprintf( fftfp, "%s\n", seq1[i] );
fprintf( fftfp, "seq2: ( %d sequences ) \n", clus2 );
for( i=0; i<clus2; i++ )
	fprintf( fftfp, "%s\n", seq2[i] );
fclose( fftfp );
system( "less input_of_Falign < /dev/tty > /dev/tty" );
#endif

	fprintf( stderr,  "FFT ... " );

	for( j=0; j<n20or4or2; j++ ) vec_init( seqVector1[j], nlen );
	if( fftscore && scoremtx != -1 )
	{
		for( i=0; i<clus1; i++ )
		{
			seq_vec_2( seqVector1[0], polarity, eff1[i], tmpseq1[i] );
			seq_vec_2( seqVector1[1], volume,   eff1[i], tmpseq1[i] );
		}
	}
	else
	{
#if 0
		for( i=0; i<clus1; i++ ) for( j=0; j<n20or4or2; j++ ) 
			seq_vec( seqVector1[j], amino[j], eff1[i], tmpseq1[i] );
#else
		for( i=0; i<clus1; i++ )
			seq_vec_3( seqVector1, eff1[i], tmpseq1[i] );
#endif
	}
#if RND
	for( i=0; i<clus1; i++ )
	{
		vec_init2( seqVector1, rndseq1[i], eff1[i], len1, nlen );
	}
#endif
#if 0
fftfp = fopen( "seqVec", "w" );
fprintf( fftfp, "before transform\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "nlen=%d\n", nlen );
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector1[k][l].R, seqVector1[k][l].I );
}
fclose( fftfp );
system( "less seqVec < /dev/tty > /dev/tty" );
#endif

	for( j=0; j<n20or4or2; j++ ) vec_init( seqVector2[j], nlen );
	if( fftscore && scoremtx != -1 )
	{
		for( i=0; i<clus2; i++ )
		{
			seq_vec_2( seqVector2[0], polarity, eff2[i], tmpseq2[i] );
			seq_vec_2( seqVector2[1], volume,   eff2[i], tmpseq2[i] );
		}
	}
	else
	{
#if 0
		for( i=0; i<clus2; i++ ) for( j=0; j<n20or4or2; j++ ) 
			seq_vec( seqVector2[j], amino[j], eff2[i], tmpseq2[i] );
#else
		for( i=0; i<clus2; i++ )
			seq_vec_3( seqVector2, eff2[i], tmpseq2[i] );
#endif
	}
#if RND
	for( i=0; i<clus2; i++ )
	{
		vec_init2( seqVector2, rndseq2[i], eff2[i], len2, nlen );
	}
#endif

#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "before fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

	for( j=0; j<n20or4or2; j++ )
	{
		fft( nlen, seqVector2[j], (j==0) );
		fft( nlen, seqVector1[j], 0 );
	}
#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "#after fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "#%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
	   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

	for( k=0; k<n20or4or2; k++ ) 
	{
		for( l=0; l<nlen; l++ ) 
			calcNaiseki( naiseki[k]+l, seqVector1[k]+l, seqVector2[k]+l );
	}
	for( l=0; l<nlen; l++ ) 
	{
		naisekiNoWa[l].R = 0.0;
		naisekiNoWa[l].I = 0.0;
		for( k=0; k<n20or4or2; k++ ) 
		{
			naisekiNoWa[l].R += naiseki[k][l].R;
			naisekiNoWa[l].I += naiseki[k][l].I;
		}
	}

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#Before fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f %f\n", l, naisekiNoWa[l].R, naisekiNoWa[l].I ); 
	fclose( fftfp );
	system( "less naisekiNoWa < /dev/tty > /dev/tty " );
#endif

	fft( -nlen, naisekiNoWa, 0 );

	for( m=0; m<=nlen2; m++ ) 
		soukan[m] = naisekiNoWa[nlen2-m].R;
	for( m=nlen2+1; m<nlen; m++ ) 
		soukan[m] = naisekiNoWa[nlen+nlen2-m].R;

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#After fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l, naisekiNoWa[l].R ); 
	fclose( fftfp );
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'naisekiNoWa'\npause -1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot &" );
#endif
#if 0
	fprintf( stderr, "frt write start\n" );
	fftfp = fopen( "frt", "w" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l-nlen2, soukan[l] ); 
	fclose( fftfp );
	system( "less frt < /dev/tty > /dev/tty" );
#if 0
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'frt'\n pause +1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot" );
#endif
#endif


	getKouho( kouho, NKOUHO, soukan, nlen );

#if 1
	for( i=0; i<NKOUHO; i++ )
	{
		fprintf( stdout, "kouho[%d] = %d\n", i, kouho[i] );
	}
#endif

#if KEIKA
	fprintf( stderr, "Searching anchors ... " );
#endif
	count = 0;



#define CAND 0
#if CAND
	fftfp = fopen( "cand", "w" );
	fclose( fftfp );
#endif

	for( k=0; k<NKOUHO; k++ ) 
	{
		lag = kouho[k];
		zurasu2( lag, clus1, clus2, seq1, seq2, tmpptr1, tmpptr2 );
#if CAND
		fftfp = fopen( "cand", "a" );
		fprintf( fftfp, "Candidate No.%d lag = %d\n", k+1, lag );
		fprintf( fftfp, "%s\n", tmpptr1[0] );
		fprintf( fftfp, "%s\n", tmpptr2[0] );
		fclose( fftfp );
#endif
		tmpint = alignableReagion( clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment+count );
		
		if( count+tmpint > MAXSEG -3 ) ErrorExit( "TOO MANY SEGMENTS.\n" );


		while( tmpint-- > 0 )
		{
			if( lag > 0 )
			{
				segment1[count].start  = segment[count].start ;
				segment1[count].end    = segment[count].end   ;
				segment1[count].center = segment[count].center;
				segment1[count].score  = segment[count].score;

				segment2[count].start  = segment[count].start  + lag;
				segment2[count].end    = segment[count].end    + lag;
				segment2[count].center = segment[count].center + lag;
				segment2[count].score  = segment[count].score       ;
			}
			else
			{
				segment1[count].start  = segment[count].start  - lag;
				segment1[count].end    = segment[count].end    - lag;
				segment1[count].center = segment[count].center - lag;
				segment1[count].score  = segment[count].score       ;

				segment2[count].start  = segment[count].start ;
				segment2[count].end    = segment[count].end   ;
				segment2[count].center = segment[count].center;
				segment2[count].score  = segment[count].score ;
			}
#if 0
			fftfp = fopen( "cand", "a" );
			fprintf( fftfp, "Goukaku=%dko\n", tmpint ); 
			fprintf( fftfp, "in 1 %d\n", segment1[count].center );
			fprintf( fftfp, "in 2 %d\n", segment2[count].center );
			fclose( fftfp );
#endif
			segment1[count].pair = &segment2[count];
			segment2[count].pair = &segment1[count];
			count++;
#if 0
			fprintf( stderr, "count=%d\n", count );
#endif
		}
	}
#if 1
	fprintf( stderr, "%d segments found\n", count );
#endif
	if( !count && fftNoAnchStop )
		ErrorExit( "Cannot detect anchor!" );
#if 1
	fprintf( stdout, "RESULT before sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( stdout, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( stdout, "%d score = %f\n", segment2[l].center, segment1[l].score );
	}
	exit( 1 );
#endif

#if KEIKA
	fprintf( stderr, "Aligning anchors ... " );
#endif
	for( i=0; i<count; i++ )
	{
		sortedseg1[i] = &segment1[i];
		sortedseg2[i] = &segment2[i];
	}
#if 0
	tmpsort( count, sortedseg1 ); 
	tmpsort( count, sortedseg2 ); 
	qsort( sortedseg1, count, sizeof( Segment * ), segcmp );
	qsort( sortedseg2, count, sizeof( Segment * ), segcmp );
#else
	mymergesort( 0, count-1, sortedseg1 ); 
	mymergesort( 0, count-1, sortedseg2 ); 
#endif
	for( i=0; i<count; i++ ) sortedseg1[i]->number = i;
	for( i=0; i<count; i++ ) sortedseg2[i]->number = i;

	if( crossscoresize < count+2 )
	{
		crossscoresize = count+2;
		fprintf( stderr, "####################################################################################################################################allocating crossscore, size = %d\n", crossscoresize );
		if( crossscore ) FreeDoubleMtx( crossscore );
		crossscore = AllocateDoubleMtx( crossscoresize, crossscoresize );
	}

	for( i=0; i<count+2; i++ ) for( j=0; j<count+2; j++ )
		crossscore[i][j] = 0.0;
	for( i=0; i<count; i++ )
	{
		crossscore[segment1[i].number+1][segment1[i].pair->number+1] = segment1[i].score;
		cut1[i+1] = sortedseg1[i]->center;
		cut2[i+1] = sortedseg2[i]->center;
	}

#if DEBUG
	fprintf( stderr, "AFTER SORT\n" );
	for( i=0; i<count; i++ ) fprintf( stderr, "%d, %d\n", segment1[i].start, segment2[i].start );
#endif

	crossscore[0][0] = 10000000.0;
	cut1[0] = 0; 
	cut2[0] = 0;
	crossscore[count+1][count+1] = 10000000.0;
	cut1[count+1] = len1;
	cut2[count+1] = len2;
	count += 2;
	count0 = count;

	blockAlign2( cut1, cut2, sortedseg1, sortedseg2, crossscore, &count );
	if( count0 > count )
	{
#if 0
		fprintf( stderr, "\7 REPEAT!? \n" ); 
#else
		fprintf( stderr, "REPEAT!? \n" ); 
#endif
		if( fftRepeatStop ) exit( 1 );
	}
#if KEIKA
	else fprintf( stderr, "done\n" );
#endif

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d\n", segment2[l].center );
	}
	fclose( fftfp );
#endif

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut : %d %d\n", cut1[l], cut2[l] );
	}
	fclose( fftfp );
#endif

#if KEIKA
	fprintf( trap_g, "Devided to %d segments\n", count-1 );
	fprintf( trap_g, "%d  %d forg\n", MIN( clus1, clus2 ), count-1 );
#endif

	totallen = 0;
	for( j=0; j<clus1; j++ ) result1[j][0] = 0;
	for( j=0; j<clus2; j++ ) result2[j][0] = 0;
	for( i=0; i<count-1; i++ )
	{
#if DEBUG
		fprintf( stderr, "DP %03d / %03d %4d to ", i+1, count-1, totallen );
#else
#if KEIKA
		fprintf( stderr, "DP %03d / %03d\r", i+1, count-1 );
#endif
#endif
		for( j=0; j<clus1; j++ )
		{
			strncpy( tmpres1[j], seq1[j]+cut1[i], cut1[i+1]-cut1[i] );
			tmpres1[j][cut1[i+1]-cut1[i]] = 0;
		}
		for( j=0; j<clus2; j++ )
		{
			strncpy( tmpres2[j], seq2[j]+cut2[i], cut2[i+1]-cut2[i] );
			tmpres2[j][cut2[i+1]-cut2[i]] = 0;
		}
		switch( alg )
		{
			case( 'a' ):
				Aalign( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			case( 'A' ):
				A__align( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			case ( 'C' ):
				if( clus1 == 1 && clus2 != 1 || clus1 != 1 && clus2 == 1 )
					translate_and_Calign( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				else
					A__align( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			default:
				fprintf( stderr, "alg = %c\n", alg );
				ErrorExit( "ERROR IN SOURCE FILE Falign.c" );
				break;
		}

		nlen = strlen( tmpres1[0] );
		if( totallen + nlen > alloclen ) ErrorExit( "LENGTH OVER in Falign\n " );
		for( j=0; j<clus1; j++ ) strcat( result1[j], tmpres1[j] );
		for( j=0; j<clus2; j++ ) strcat( result2[j], tmpres2[j] );
		totallen += nlen;
#if 0
		fprintf( stderr, "%4d\r", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
	}
#if KEIKA
	fprintf( stderr, "DP ... done   \n" );
#endif

	for( j=0; j<clus1; j++ ) strcpy( seq1[j], result1[j] );
	for( j=0; j<clus2; j++ ) strcpy( seq2[j], result2[j] );
#if 0
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, "%s\n", result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, "%s\n", result2[j] );
	}
#endif
	return( 0.0 );
}
float Falign( char  **seq1, char  **seq2, 
			  double *eff1, double *eff2, 
			  int    clus1, int    clus2,
			  int alloclen )
{
	int i, j, k, l, m;
	int nlen, nlen2, nlen4;
	static int crossscoresize = 0;
	static char **tmpseq1 = NULL;
	static char **tmpseq2 = NULL;
	static char **tmpptr1 = NULL;
	static char **tmpptr2 = NULL;
	static char **tmpres1 = NULL;
	static char **tmpres2 = NULL;
	static char **result1 = NULL;
	static char **result2 = NULL;
#if RND
	static char **rndseq1 = NULL;
	static char **rndseq2 = NULL;
#endif
	static Fukusosuu **seqVector1 = NULL;
	static Fukusosuu **seqVector2 = NULL;
	static Fukusosuu **naiseki = NULL;   
	static Fukusosuu *naisekiNoWa = NULL; 
	static double *soukan = NULL;
	static double **crossscore = NULL;
	int nlentmp;
	static int *kouho = NULL;
	static Segment *segment = NULL;
	static Segment *segment1 = NULL;
	static Segment *segment2 = NULL;
	static Segment **sortedseg1 = NULL;
	static Segment **sortedseg2 = NULL;
	static int *cut1 = NULL;
	static int *cut2 = NULL;
	static int localalloclen = 0;
	int lag;
	int tmpint;
	int count, count0;
	int len1, len2;
	int totallen;
	double totaleff1, totaleff2;
	float totalscore;

	extern Fukusosuu   *AllocateFukusosuuVec();
	extern Fukusosuu  **AllocateFukusosuuMtx();

	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

	len1 = strlen( seq1[0] );
	len2 = strlen( seq2[0] );
	nlentmp = MAX( len1, len2 );

	nlen = 1;
	while( nlentmp >= nlen ) nlen <<= 1;
#if 0
	fprintf( stderr, "###   nlen    = %d\n", nlen );
#endif

	nlen2 = nlen/2; nlen4 = nlen2 / 2;

#if DEBUG
	fprintf( stderr, "len1 = %d, len2 = %d\n", len1, len2 );
	fprintf( stderr, "nlentmp = %d, nlen = %d\n", nlentmp, nlen );
#endif

	if( !localalloclen )
	{
		kouho = AllocateIntVec( NKOUHO );
		cut1 = AllocateIntVec( MAXSEG );
		cut2 = AllocateIntVec( MAXSEG );
		tmpptr1 = AllocateCharMtx( njob, 1 );
		tmpptr2 = AllocateCharMtx( njob, 1 );
		result1 = AllocateCharMtx( njob, alloclen );
		result2 = AllocateCharMtx( njob, alloclen );
		tmpres1 = AllocateCharMtx( njob, alloclen );
		tmpres2 = AllocateCharMtx( njob, alloclen );
//		crossscore = AllocateDoubleMtx( MAXSEG, MAXSEG );
		segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment1 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment2 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		sortedseg1 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		sortedseg2 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		if( !( segment && segment1 && segment2 && sortedseg1 && sortedseg2 ) )
			ErrorExit( "Allocation error\n" );

		if     ( scoremtx == -1 ) n20or4or2 = 4;
		else if( fftscore == 1  ) n20or4or2 = 2;
		else                      n20or4or2 = 20;
	}
	if( localalloclen < nlen )
	{
		if( localalloclen )
		{
#if 1
			FreeFukusosuuMtx ( seqVector1 );
			FreeFukusosuuMtx ( seqVector2 );
			FreeFukusosuuVec( naisekiNoWa );
			FreeFukusosuuVec( naiseki );
			FreeDoubleVec( soukan );
			FreeCharMtx( tmpseq1 );
			FreeCharMtx( tmpseq2 );
#endif
#if RND
			FreeCharMtx( rndseq1 );
			FreeCharMtx( rndseq2 );
#endif
		}

		tmpseq1 = AllocateCharMtx( njob, nlen );
		tmpseq2 = AllocateCharMtx( njob, nlen );
		naisekiNoWa = AllocateFukusosuuVec( nlen );
		naiseki = AllocateFukusosuuMtx( n20or4or2, nlen );
		seqVector1 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
		seqVector2 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
		soukan = AllocateDoubleVec( nlen+1 );
#if RND
		rndseq1 = AllocateCharMtx( njob, nlen );
		rndseq2 = AllocateCharMtx( njob, nlen );
		for( i=0; i<njob; i++ )
		{
			generateRndSeq( rndseq1[i], nlen );
			generateRndSeq( rndseq2[i], nlen );
		}
#endif
		localalloclen = nlen;
	}
	
	for( j=0; j<clus1; j++ ) strcpy( tmpseq1[j], seq1[j] );
	for( j=0; j<clus2; j++ ) strcpy( tmpseq2[j], seq2[j] );

#if 0
fftfp = fopen( "input_of_Falign", "w" );
fprintf( fftfp, "nlen = %d\n", nlen );
fprintf( fftfp, "seq1: ( %d sequences ) \n", clus1 );
for( i=0; i<clus1; i++ )
	fprintf( fftfp, "%s\n", seq1[i] );
fprintf( fftfp, "seq2: ( %d sequences ) \n", clus2 );
for( i=0; i<clus2; i++ )
	fprintf( fftfp, "%s\n", seq2[i] );
fclose( fftfp );
system( "less input_of_Falign < /dev/tty > /dev/tty" );
#endif

	fprintf( stderr,  "FFT ... " );

	for( j=0; j<n20or4or2; j++ ) vec_init( seqVector1[j], nlen );
	if( fftscore && scoremtx != -1 )
	{
		for( i=0; i<clus1; i++ )
		{
			seq_vec_2( seqVector1[0], polarity, eff1[i], tmpseq1[i] );
			seq_vec_2( seqVector1[1], volume,   eff1[i], tmpseq1[i] );
		}
	}
	else
	{
#if 0
		for( i=0; i<clus1; i++ ) for( j=0; j<n20or4or2; j++ ) 
			seq_vec( seqVector1[j], amino[j], eff1[i], tmpseq1[i] );
#else
		for( i=0; i<clus1; i++ )
			seq_vec_3( seqVector1, eff1[i], tmpseq1[i] );
#endif
	}
#if RND
	for( i=0; i<clus1; i++ )
	{
		vec_init2( seqVector1, rndseq1[i], eff1[i], len1, nlen );
	}
#endif
#if 0
fftfp = fopen( "seqVec", "w" );
fprintf( fftfp, "before transform\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "nlen=%d\n", nlen );
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector1[k][l].R, seqVector1[k][l].I );
}
fclose( fftfp );
system( "less seqVec < /dev/tty > /dev/tty" );
#endif

	for( j=0; j<n20or4or2; j++ ) vec_init( seqVector2[j], nlen );
	if( fftscore && scoremtx != -1 )
	{
		for( i=0; i<clus2; i++ )
		{
			seq_vec_2( seqVector2[0], polarity, eff2[i], tmpseq2[i] );
			seq_vec_2( seqVector2[1], volume,   eff2[i], tmpseq2[i] );
		}
	}
	else
	{
#if 0
		for( i=0; i<clus2; i++ ) for( j=0; j<n20or4or2; j++ ) 
			seq_vec( seqVector2[j], amino[j], eff2[i], tmpseq2[i] );
#else
		for( i=0; i<clus2; i++ )
			seq_vec_3( seqVector2, eff2[i], tmpseq2[i] );
#endif
	}
#if RND
	for( i=0; i<clus2; i++ )
	{
		vec_init2( seqVector2, rndseq2[i], eff2[i], len2, nlen );
	}
#endif

#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "before fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

	for( j=0; j<n20or4or2; j++ )
	{
		fft( nlen, seqVector2[j], (j==0) );
		fft( nlen, seqVector1[j], 0 );
	}
#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "#after fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "#%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
	   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

	for( k=0; k<n20or4or2; k++ ) 
	{
		for( l=0; l<nlen; l++ ) 
			calcNaiseki( naiseki[k]+l, seqVector1[k]+l, seqVector2[k]+l );
	}
	for( l=0; l<nlen; l++ ) 
	{
		naisekiNoWa[l].R = 0.0;
		naisekiNoWa[l].I = 0.0;
		for( k=0; k<n20or4or2; k++ ) 
		{
			naisekiNoWa[l].R += naiseki[k][l].R;
			naisekiNoWa[l].I += naiseki[k][l].I;
		}
	}

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#Before fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f %f\n", l, naisekiNoWa[l].R, naisekiNoWa[l].I ); 
	fclose( fftfp );
	system( "less naisekiNoWa < /dev/tty > /dev/tty " );
#endif

	fft( -nlen, naisekiNoWa, 0 );

	for( m=0; m<=nlen2; m++ ) 
		soukan[m] = naisekiNoWa[nlen2-m].R;
	for( m=nlen2+1; m<nlen; m++ ) 
		soukan[m] = naisekiNoWa[nlen+nlen2-m].R;

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#After fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l, naisekiNoWa[l].R ); 
	fclose( fftfp );
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'naisekiNoWa'\npause -1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot &" );
#endif
#if 0
	fprintf( stderr, "frt write start\n" );
	fftfp = fopen( "frt", "w" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l-nlen2, soukan[l] ); 
	fclose( fftfp );
	system( "less frt < /dev/tty > /dev/tty" );
#if 0
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'frt'\n pause +1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot" );
#endif
#endif


	getKouho( kouho, NKOUHO, soukan, nlen );

#if 0
	for( i=0; i<NKOUHO; i++ )
	{
		fprintf( stderr, "kouho[%d] = %d\n", i, kouho[i] );
	}
#endif

#if KEIKA
	fprintf( stderr, "Searching anchors ... " );
#endif
	count = 0;



#define CAND 0
#if CAND
	fftfp = fopen( "cand", "w" );
	fclose( fftfp );
#endif

	for( k=0; k<NKOUHO; k++ ) 
	{
		lag = kouho[k];
		zurasu2( lag, clus1, clus2, seq1, seq2, tmpptr1, tmpptr2 );
#if CAND
		fftfp = fopen( "cand", "a" );
		fprintf( fftfp, "Candidate No.%d lag = %d\n", k+1, lag );
		fprintf( fftfp, "%s\n", tmpptr1[0] );
		fprintf( fftfp, "%s\n", tmpptr2[0] );
		fclose( fftfp );
#endif
		tmpint = alignableReagion( clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment+count );
		
		if( count+tmpint > MAXSEG -3 ) ErrorExit( "TOO MANY SEGMENTS.\n" );


		while( tmpint-- > 0 )
		{
			if( lag > 0 )
			{
				segment1[count].start  = segment[count].start ;
				segment1[count].end    = segment[count].end   ;
				segment1[count].center = segment[count].center;
				segment1[count].score  = segment[count].score;

				segment2[count].start  = segment[count].start  + lag;
				segment2[count].end    = segment[count].end    + lag;
				segment2[count].center = segment[count].center + lag;
				segment2[count].score  = segment[count].score       ;
			}
			else
			{
				segment1[count].start  = segment[count].start  - lag;
				segment1[count].end    = segment[count].end    - lag;
				segment1[count].center = segment[count].center - lag;
				segment1[count].score  = segment[count].score       ;

				segment2[count].start  = segment[count].start ;
				segment2[count].end    = segment[count].end   ;
				segment2[count].center = segment[count].center;
				segment2[count].score  = segment[count].score ;
			}
#if 0
			fftfp = fopen( "cand", "a" );
			fprintf( fftfp, "Goukaku=%dko\n", tmpint ); 
			fprintf( fftfp, "in 1 %d\n", segment1[count].center );
			fprintf( fftfp, "in 2 %d\n", segment2[count].center );
			fclose( fftfp );
#endif
			segment1[count].pair = &segment2[count];
			segment2[count].pair = &segment1[count];
			count++;
#if 0
			fprintf( stderr, "count=%d\n", count );
#endif
		}
	}
#if 1
	fprintf( stderr, "%d segments found\n", count );
#endif
	if( !count && fftNoAnchStop )
		ErrorExit( "Cannot detect anchor!" );
#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT before sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d score = %f\n", segment2[l].center, segment1[l].score );
	}
	fclose( fftfp );
#endif

#if KEIKA
	fprintf( stderr, "Aligning anchors ... " );
#endif
	for( i=0; i<count; i++ )
	{
		sortedseg1[i] = &segment1[i];
		sortedseg2[i] = &segment2[i];
	}
#if 0
	tmpsort( count, sortedseg1 ); 
	tmpsort( count, sortedseg2 ); 
	qsort( sortedseg1, count, sizeof( Segment * ), segcmp );
	qsort( sortedseg2, count, sizeof( Segment * ), segcmp );
#else
	mymergesort( 0, count-1, sortedseg1 ); 
	mymergesort( 0, count-1, sortedseg2 ); 
#endif
	for( i=0; i<count; i++ ) sortedseg1[i]->number = i;
	for( i=0; i<count; i++ ) sortedseg2[i]->number = i;

	if( crossscoresize < count+2 )
	{
		crossscoresize = count+2;
#if 1
		fprintf( stderr, "######allocating crossscore, size = %d\n", crossscoresize );
#endif
		if( crossscore ) FreeDoubleMtx( crossscore );
		crossscore = AllocateDoubleMtx( crossscoresize, crossscoresize );
	}
	for( i=0; i<count+2; i++ ) for( j=0; j<count+2; j++ )
		crossscore[i][j] = 0.0;
	for( i=0; i<count; i++ )
	{
		crossscore[segment1[i].number+1][segment1[i].pair->number+1] = segment1[i].score;
		cut1[i+1] = sortedseg1[i]->center;
		cut2[i+1] = sortedseg2[i]->center;
	}

#if DEBUG
	fprintf( stderr, "AFTER SORT\n" );
	for( i=0; i<count; i++ ) fprintf( stderr, "%d, %d\n", segment1[i].start, segment2[i].start );
#endif

	crossscore[0][0] = 10000000.0;
	cut1[0] = 0; 
	cut2[0] = 0;
	crossscore[count+1][count+1] = 10000000.0;
	cut1[count+1] = len1;
	cut2[count+1] = len2;
	count += 2;
	count0 = count;

	blockAlign2( cut1, cut2, sortedseg1, sortedseg2, crossscore, &count );
	if( count0 > count )
	{
#if 0
		fprintf( stderr, "\7 REPEAT!? \n" ); 
#else
		fprintf( stderr, "REPEAT!? \n" ); 
#endif
		if( fftRepeatStop ) exit( 1 );
	}
#if KEIKA
	else fprintf( stderr, "done\n" );
#endif

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d\n", segment2[l].center );
	}
	fclose( fftfp );
#endif

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut : %d %d\n", cut1[l], cut2[l] );
	}
	fclose( fftfp );
#endif

#if KEIKA
	fprintf( trap_g, "Devided to %d segments\n", count-1 );
	fprintf( trap_g, "%d  %d forg\n", MIN( clus1, clus2 ), count-1 );
#endif

	totallen = 0;
	for( j=0; j<clus1; j++ ) result1[j][0] = 0;
	for( j=0; j<clus2; j++ ) result2[j][0] = 0;
	totalscore = 0.0;
	for( i=0; i<count-1; i++ )
	{
#if DEBUG
		fprintf( stderr, "DP %03d / %03d %4d to ", i+1, count-1, totallen );
#else
#if KEIKA
		fprintf( stderr, "DP %03d / %03d\r", i+1, count-1 );
#endif
#endif
		for( j=0; j<clus1; j++ )
		{
			strncpy( tmpres1[j], seq1[j]+cut1[i], cut1[i+1]-cut1[i] );
			tmpres1[j][cut1[i+1]-cut1[i]] = 0;
		}
		for( j=0; j<clus2; j++ )
		{
			strncpy( tmpres2[j], seq2[j]+cut2[i], cut2[i+1]-cut2[i] );
			tmpres2[j][cut2[i+1]-cut2[i]] = 0;
		}
		switch( alg )
		{
			case( 'a' ):
				totalscore += Aalign( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			case( 'A' ):
				totalscore += A__align( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			case ( 'C' ):
				if( clus1 == 1 && clus2 != 1 || clus1 != 1 && clus2 == 1 )
					totalscore += translate_and_Calign( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				else
					totalscore += A__align( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			default:
				fprintf( stderr, "alg = %c\n", alg );
				ErrorExit( "ERROR IN SOURCE FILE Falign.c" );
				break;
		}

		nlen = strlen( tmpres1[0] );
		if( totallen + nlen > alloclen )
		{
			fprintf( stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", totallen, nlen, alloclen );
			ErrorExit( "LENGTH OVER in Falign\n " );
		}
		for( j=0; j<clus1; j++ ) strcat( result1[j], tmpres1[j] );
		for( j=0; j<clus2; j++ ) strcat( result2[j], tmpres2[j] );
		totallen += nlen;
#if 0
		fprintf( stderr, "%4d\r", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
	}
#if KEIKA
	fprintf( stderr, "DP ... done   \n" );
#endif

	for( j=0; j<clus1; j++ ) strcpy( seq1[j], result1[j] );
	for( j=0; j<clus2; j++ ) strcpy( seq2[j], result2[j] );
#if 0
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, "%s\n", result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, "%s\n", result2[j] );
	}
#endif
	return( totalscore );
}
