#include "mltaln.h"

static int upperCase = 0;

#define DEBUG   0
#define IODEBUG 0

void strncpy_caseC( char *str1, char *str2, int len )
{
	if( scoremtx == -1 && upperCase > 0 ) 
	{
		while( len-- )
			*str1++ = toupper( *str2++ );
	}
	else strncpy( str1, str2, len );
}
	
void seqUpper( int nseq, char **seq ) /* not used */
{
	int i, j, len;
	for( i=0; i<nseq; i++ ) 
	{
		len = strlen( seq[i] );
		for( j=0; j<len; j++ ) 
			seq[i][j] = toupper( seq[i][j] );
	}
}

void seqLower( int nseq, char **seq )
{
	int i, j, len;
	for( i=0; i<nseq; i++ ) 
	{
		len = strlen( seq[i] );
		for( j=0; j<len; j++ ) 
			seq[i][j] = tolower( seq[i][j] );
	}
}

int getaline_fp_eof( char *s, int l, FILE *fp )  /* end of file -> return 1 */
{
    int c, i = 0 ;
    int noteofflag;
    for( i=0; i<l && ( noteofflag = ( (c=getc(fp)) != EOF ) ) && c != '\n'; i++ ) 
    	*s++ = c;
    *s = '\0' ;
     return( !noteofflag );
}

int getaline_fp_eof_new(s, l, fp)  /* end of file -> return 1 */
char    s[] ; int l ; FILE *fp ;
{
        int     c, i = 0 ;
		int noteofflag;

		if( feof( fp ) ) return( 1 );

		for( i=0; i<l && ( noteofflag = ( (c=getc(fp)) != EOF ) ) && c != '\n'; i++ ) 
        { *s++ = c ; }
        *s = '\0' ;
		if( c != '\n' && c != EOF ) while( getc(fp) != '\n' )
			;
		return( !noteofflag );
}

int myfgets(s, l, fp)  /* l以上は、行末まで読み飛ばす */
char    s[] ; int l ; FILE *fp ;
{
        int     c, i = 0 ;

		if( feof( fp ) ) return( 1 );

		for( i=0; i<l && ( c=getc( fp ) ) != '\n'; i++ ) 
        	*s++ = c;
        *s = '\0' ;
		if( c != '\n' ) 
			while( getc(fp) != '\n' )
				;
}

#if 0
float input( FILE *fp, int lmax, int d )
{
    int i;
    char b[B];
    static char x[M*M/2*D+1];
    static int f = 0;
    static int l = 0;

    if( f == 0 )
    {
        f = 1;
        x[0] = 0;
        for( i=0; i<lmax*d; )
        {
            int j;
            fgets( b, B, fp );
            j = strlen( b ) - 1;
            b[j] = 0;
            i += j;
            strcat( x, b );
        }
    }

    strncpy( b, x+l*d, d ); b[d] = 0;
    l++;
    if( l == lmax )
    {
        f = 0; l = 0;
    }
    return( atof( b ) );
}
#endif

float input_new( FILE *fp, int d )
{
	char mojiretsu[10];
	int i, c;

	c = getc( fp );
	if( c != '\n' )
		ungetc( c, fp );

	for( i=0; i<d; i++ )
		mojiretsu[i] = getc( fp );

	return( atof( mojiretsu ) );
}


void PreRead( FILE *fp, int *locnjob, int *locnlenmax )
{
	int i, nleni;
	char b[B];

	fgets( b, B-1, fp ); *locnjob = atoi( b );
	*locnlenmax = 0;
	i=0; 
	while( i<*locnjob )
	{
		fgets( b, B-1, fp );
		if( !strncmp( b, "=", 1 ) )
		{
			i++;
			fgets( b, B-1, fp ); nleni = atoi( b );
			if( nleni > *locnlenmax ) *locnlenmax = nleni;
		}
	}
	if( *locnlenmax > N )
	{
		fprintf( stderr, "TOO LONG SEQUENCE!\n" );
		exit( 1 );
	}
	if( njob > M  ) 
	{
		fprintf( stderr, "TOO MANY SEQUENCE!\n" );
		fprintf( stderr, "%d > %d\n", njob, M );
		exit( 1 );
	}
}	

int allSpace( char *str )
{
	int value = 1;
	while( *str ) value *= ( !isdigit( *str++ ) );
	return( value );
}
	
void Read( char name[M][B], int nlen[M], char **seq )
{
	extern void FRead( FILE *x, char y[M][B], int z[M], char **w );
	FRead( stdin, name, nlen, seq );
}

void FRead( FILE *fp, char name[][B], int nlen[], char **seq )
{
	int i, j; 
	char b[B];

	fgets( b, B-1, fp );
#if DEBUG
	fprintf( stderr, "b = %s\n", b );
#endif

    if( strstr( b, "onnet" ) ) scoremtx = 1;
    else if( strstr( b, "DnA" ) ) 
	{
		scoremtx = -1; 
		upperCase = -1;
	}
    else if( strstr( b, "dna" ) ) 
	{
		scoremtx = -1; 
		upperCase = 0;
	}
	else if( strstr( b, "DNA" ) )
	{
		scoremtx = -1; 
		upperCase = 1;
	}
    else if( strstr( b, "M-Y" ) || strstr( b, "iyata" ) ) scoremtx = 2; 
    else scoremtx = 0;
#if DEBUG
	fprintf( stderr, " %s->scoremtx = %d\n", b, scoremtx );
#endif

	geta2 = GETA2;

#if 0
	if( strlen( b ) >=25 )
	{
		b[25] = 0;
	#if DEBUG
		fprintf( stderr, "kimuraR = %s\n", b+20 );
	#endif
		kimuraR = atoi( b+20 );

		if( kimuraR < 0 || 20 < kimuraR ) ErrorExit( "Illeagal kimuraR value.\n" );
		if( allSpace( b+20 ) ) kimuraR = NOTSPECIFIED;
	}
	else kimuraR = NOTSPECIFIED;
	#if DEBUG
	fprintf( stderr, "kimuraR = %d\n", kimuraR );
	#endif

	if( strlen( b ) >=20 )
	{
		b[20] = 0;
	#if DEBUG
		fprintf( stderr, "pamN = %s\n", b+15 );
	#endif
		pamN = atoi( b+15 );
		if( pamN < 0 || 400 < pamN ) ErrorExit( "Illeagal pam value.\n" );
		if( allSpace( b+15 ) ) pamN = NOTSPECIFIED;
	}
	else pamN = NOTSPECIFIED;

	if( strlen( b ) >= 15 )
	{
		b[15] = 0;
	#if DEBUG
		fprintf( stderr, "poffset = %s\n", b+10 );
	#endif
		poffset = atoi( b+10 );
		if( poffset > 500 ) ErrorExit( "Illegal extending gap ppenalty\n" );
		if( allSpace( b+10 ) ) poffset = NOTSPECIFIED;
	}
	else poffset = NOTSPECIFIED;

	if( strlen( b ) >= 10 )
	{
		b[10] = 0;
	#if DEBUG
		fprintf( stderr, "ppenalty = %s\n", b+5 );
	#endif
		ppenalty = atoi( b+5 );
		if( ppenalty > 0 ) ErrorExit( "Illegal opening gap ppenalty\n" );
		if( allSpace( b+5 ) ) ppenalty = NOTSPECIFIED;
	}
	else ppenalty = NOTSPECIFIED;
#endif

	for( i=0; i<njob; i++ )
	{
		getaline_fp_eof_new( b, B-1, fp );
		strcpy( name[i], b );
#if DEBUG
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		fgets( b, B-1, fp ); nlen[i] = atoi( b );      /* seq i no nagasa */
		seq[i][0] = 0;
		if( nlen[i] ) for( j=0; j <= (nlen[i]-1)/C; j++ )
		{
			getaline_fp_eof_new( b, B-1, fp );
			/*	b[C] = 0;  */
			strcat( seq[i], b );
		} 
		seq[i][nlen[i]] = 0;
	}
	if( scoremtx == -1 && upperCase != -1 ) seqLower( njob, seq );
}


static int countKUorWA( FILE *fp )
{
	int value;
	int c, b;

	value= 0;
	b = '\n';
	while( ( c = getc( fp ) ) != EOF )
	{
		if( b == '\n' && ( c == '=' || c == '>' ) )
			value++;
		b = c;
	}
	rewind( fp );
	return( value );
}

static void searchKUorWA( FILE *fp )
{
	int c, b;
	b = '\n';
	while( !( ( ( c = getc( fp ) ) == '>' || c == '=' || c == EOF ) && b == '\n' ) )
		b = c;
	ungetc( c, fp );
}

static int onlyAlpha_lower( char *str )
{
	char tmp;
	char *res = str;
	char *bk = str;

	while( (tmp=*str++) )
		if( isalpha( tmp ) || tmp == '-' || tmp == '*' || tmp == '.' )
			*res++ = tolower( tmp );
	*res = 0;
	return( res - bk );
}
static int onlyAlpha_upper( char *str )
{
	char tmp;
	char *res = str;
	char *bk = str;

	while( (tmp=*str++) )
		if( isalpha( tmp ) || tmp == '-' || tmp == '*' || tmp == '.' )
			*res++ = toupper( tmp );
	*res = 0;
	return( res - bk );
}

void kake2hiku( char *str )
{
	do
		if( *str == '*' ) *str = '-';
	while( *str++ );
}

static int load1SeqWithoutName_new( FILE *fpp, char *cbuf )
{
	int c, b;
	char *bk = cbuf;

	b = '\n';
	while( ( c = getc( fpp ) ) != EOF &&                    /* by T. Nishiyama */
          !( ( c == '>' || c == '=' || c == '(' || c == EOF ) && b == '\n' ) )
	{
		*cbuf++ = (char)c;  /* 長すぎてもしらない */
		b = c;
	}
	ungetc( c, fpp );
	*cbuf = 0;
	if( scoremtx == -1 )
		onlyAlpha_lower( bk );
	else
		onlyAlpha_upper( bk );
	kake2hiku( bk );
	return( 0 );
}


void readDataforgaln( FILE *fp, char **name, int *nlen, char **seq )
{
	int i, j; 
	static char *tmpseq = NULL;

	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		load1SeqWithoutName_new( fp, tmpseq );
		strcpy( seq[i], tmpseq );
		nlen[i] = strlen( seq[i] );
	}
	if( scoremtx == -1 && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
}

void readData( FILE *fp, char name[][B], int nlen[], char **seq )
{
	int i, j; 
	static char *tmpseq = NULL;

	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		load1SeqWithoutName_new( fp, tmpseq );
		strcpy( seq[i], tmpseq );
		nlen[i] = strlen( seq[i] );
	}
	if( scoremtx == -1 && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
}


double countATGC( char *s )
{
	int nATGC;
	int nChar;
	char c;
	nATGC = nChar = 0;

	do
	{
		c = tolower( *s );
		if( isalpha( c ) )
		{
			nChar++;
			if( c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'u' || c == 'n' )
				nATGC++;
		}
	}
	while( *++s );
	return( (double)nATGC / nChar );
}


void getnumlen( FILE *fp )
{
	int i, tmp;
	char *tmpseq;
	double atgcfreq;
	tmpseq = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	atgcfreq = 0.0;
	for( i=0; i<njob; i++ )
	{
		fgets( tmpseq, N-1, fp );
		load1SeqWithoutName_new( fp, tmpseq );
		tmp = strlen( tmpseq );
		if( tmp > nlenmax ) nlenmax  = tmp;
		atgcfreq += countATGC( tmpseq );
	}
	atgcfreq /= (double)njob;
	if( scoremtx == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) scoremtx =  upperCase = -1;
		else                  scoremtx = 0;
	}
	free( tmpseq );
}
	


void WriteGapFill( FILE *fp, int locnjob, char name[][B], int nlen[M], char **aseq )
{
	static char b[N];
	int i, j;
	int nalen[M];
	static char gap[N];
	static char buff[N];

#if IODEBUG
	fprintf( stderr, "IMAKARA KAKU\n" );
#endif
	nlenmax = 0;
	for( i=0; i<locnjob; i++ )
	{
		int len = strlen( aseq[i] );
		if( nlenmax < len ) nlenmax = len;
	}

	for( i=0; i<nlenmax; i++ ) gap[i] = '-';
	gap[nlenmax] = 0;

	fprintf( fp, "%5d", locnjob );
	fprintf( fp, "\n" );

	for( i=0; i<locnjob; i++ )
	{
		strcpy( buff, aseq[i] );
		strncat( buff, gap, nlenmax-strlen( aseq[i] ) );
		buff[nlenmax] = 0;
		nalen[i] = strlen( buff );
		fprintf( fp, "%s\n", name[i] );
		fprintf( fp, "%5d\n", nalen[i] );
		for( j=0; j<nalen[i]; j=j+C )
		{
			strncpy_caseC( b, buff+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
		}
	}
#if DEBUG
	fprintf( stderr, "nalen[0] = %d\n", nalen[0] );
#endif
#if IODEBUG
	fprintf( stderr, "KAKIOWATTA\n" );
#endif
}

void writeDataforgaln( FILE *fp, int locnjob, char **name, int *nlen, char **aseq )
{
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		nalen = strlen( aseq[i] );
		fprintf( fp, ">%s\n", name[i]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[i]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[i]+j );
#endif
		}
	}
}

void writeData( FILE *fp, int locnjob, char name[][B], int nlen[], char **aseq )
{
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		nalen = strlen( aseq[i] );
		fprintf( fp, ">%s\n", name[i]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[i]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[i]+j );
#endif
		}
	}
}





void readhat2( FILE *fp, int nseq, char name[M][B], double **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = (double)input_new( fp, D);
    }
}

void WriteFloatHat2( FILE *hat2p, int locnjob, char name[M][B], float **mtx )
{
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob-1; i++ )
	{
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, "%#6.3f", mtx[i][j] );
			if( (j-i) % 12 == 0 || j == locnjob-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2( FILE *hat2p, int locnjob, char name[M][B], double **mtx )
{
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob-1; i++ )
	{
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, "%#6.3f", mtx[i][j] );
			if( (j-i) % 12 == 0 || j == locnjob-1 ) fprintf( hat2p, "\n" );
		}
	}
}

int ReadFasta_sub( FILE *fp, double *dis, int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
            count++;
        }
    }

	for( i=0; i<nseq; i++ ) dis[i] = 0.0;
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
		if( fgets( b, B-1, fp ) ) break;
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            fgets( b, B-1, fp );
            dis[junban[count]] = atof( b );
            count++;
        }
    }
    return 0;
}


int ReadSsearch( FILE *fp, double *dis, int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];
	int opt;

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
			sscanf( b+75, "%d", &opt ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }

/*
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            dis[junban[count]] = atof( b+65 );
            count++;
        }
    }
*/
    return 0;
}

int ReadFasta3( FILE *fp, double *dis, int nseq, char name[M][B] )
{
    int count=0;
    char b[B];
	char *pt;
    int junban[M];
	int initn, init1, opt;
	double z;

    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %d %d %lf", &initn, &init1, &opt, &z ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }
    return 0;
}

int ReadFasta( FILE *fp, double *dis, int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];
	int initn, init1, opt;

    count = 0;
	for( i=0; i<nseq; i++ ) dis[i] = 0.0;
    for( i=0; !feof( fp ) && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			sscanf( b+50, "%d %d %d", &initn, &init1, &opt ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }

/*
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            dis[junban[count]] = atof( b+65 );
            count++;
        }
    }
*/
    return 0;
}


int ReadOpt( FILE *fp, int opt[M], int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];
	int optt, initn, init1;

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
			sscanf( b+50, "%d %d %d", &initn, &init1, &optt ); 
            opt[junban[count]] = (double)optt;
            count++;
        }
    }
    return 0;
}

int ReadOpt2( FILE *fp, int opt[M], int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
            opt[junban[count]] = atoi( b+65 );
            count++;
        }
    }
    return 0;
}


int ErrorExit( char *message )
{
	fprintf( stderr, "%s\n", message );
	exit( 1 );
}

int writePre( int nseq, char name[][B], int nlen[M], char **aseq, int force )
{
#if USE_XCED
	int i, value;
	if( !signalSM )
	{
		if( force ) 
		{
			rewind( prep_g );
			if( devide ) dvWrite( prep_g, nseq, name, nlen, aseq );
#if 0
			else    WriteGapFill( prep_g, nseq, name, nlen, aseq );
#else
			else    writeData( prep_g, nseq, name, nlen, aseq );
#endif
		}
		return( 0 );
	}
	for( i=0; i<10; i++ )
	{
#if IODEBUG
		fprintf( stderr, "SEMAPHORE = %d\n", signalSM[SEMAPHORE] );
#endif
		if( signalSM[SEMAPHORE]-- > 0 )
		{
#if 0 /* /tmp/pre の関係ではずした */
			if( ferror( prep_g ) ) prep_g = fopen( "pre", "w" );
			if( !prep_g ) ErrorExit( "Cannot re-open pre." ); 
#endif
			rewind( prep_g );
			signalSM[STATUS] = IMA_KAITERU;
#if IODEBUG
			if( force ) fprintf( stderr, "FINAL " );
#endif
			if( devide ) dvWrite( prep_g, nseq, name, nlen, aseq );
			else    WriteGapFill( prep_g, nseq, name, nlen, aseq );
			/*
			fprintf( prep_g, '\EOF' );
			*/
			fflush( prep_g );
			if( force ) signalSM[STATUS] = OSHIMAI;
			else        signalSM[STATUS] = KAKIOWATTA;
			value = 1;
			signalSM[SEMAPHORE]++;
#if IODEBUG
			fprintf( stderr, "signalSM[STATUS] = %c\n", signalSM[STATUS] );
#endif
			break;
		}
		else
		{
#if IODEBUG
			fprintf( stderr, "YONDERUKARA_AKIRAMERU\n" );
#endif
			value = 0;
			signalSM[SEMAPHORE]++;
			if( !force ) break;
#if IODEBUG
			fprintf( stderr, "MATSU\n" );
#endif
			sleep( 1 );
		}
	}
	if( force && !value ) ErrorExit( "xced ga pre wo hanasanai \n" );
	return( value );
#else
	if( force ) 
	{
		rewind( prep_g );
			writeData( prep_g, nseq, name, nlen, aseq );
	}
#endif
}

#include "fft.h"

void readOtherOptions( int *ppidptr, int *fftThresholdptr, int *fftWinSizeptr )
{
	if( calledByXced )
	{
		FILE *fp = fopen( "pre", "r" );
		char b[B];
		if( !fp ) ErrorExit( "Cannot open pre.\n" );
		fgets( b, B-1, fp );
		sscanf( b, "%d %d %d", ppidptr, fftThresholdptr, fftWinSizeptr );
		fclose( fp );
#if IODEBUG
	fprintf( stderr, "b = %s\n", b );
	fprintf( stderr, "ppid = %d\n", ppid );
	fprintf( stderr, "fftThreshold = %d\n", fftThreshold );
	fprintf( stderr, "fftWinSize = %d\n", fftWinSize );
#endif
	}
	else
	{
		*ppidptr = 0;
		*fftThresholdptr = FFT_THRESHOLD;
		if( scoremtx == -1 )
			*fftWinSizeptr = FFT_WINSIZE_D;
		else
			*fftWinSizeptr = FFT_WINSIZE_P;
	}
#if 0
	fprintf( stderr, "fftThresholdptr=%d\n", *fftThresholdptr );
	fprintf( stderr, "fftWinSizeptr=%d\n", *fftWinSizeptr );
#endif
}

void initSignalSM()
{
	int signalsmid;

#if IODEBUG
	if( ppid ) fprintf( stderr, "PID of xced = %d\n", ppid );
#endif
	if( !ppid )
	{
		signalSM = NULL;
		return;
	}

#if 0
	signalsmid = shmget( (key_t)ppid, 3, IPC_ALLOC | 0666 );
	if( signalsmid == -1 ) ErrorExit( "Cannot get Shared memory for signal.\n" );
	signalSM = shmat( signalsmid, 0, 0 );
	if( (int)signalSM == -1 ) ErrorExit( "Cannot attatch Shared Memory for signal!\n" );
	signalSM[STATUS] = IMA_KAITERU;
	signalSM[SEMAPHORE] = 1;
#endif
}

void initFiles( void )
{
	char pname[100];
	if( ppid )
		sprintf( pname, "/tmp/pre.%d", ppid );
	else
		sprintf( pname, "pre" );
	fprintf( stderr, "pre in align = %s\n", pname );
	prep_g = fopen( pname, "w" );
	if( !prep_g ) ErrorExit( "Cannot open pre" );

	trap_g = fopen( "trace", "w" );
	if( !trap_g ) ErrorExit( "cannot open trace" );
	fprintf( trap_g, "PID = %d\n", getpid() );
	fflush( trap_g );
}

/* not used */
#if 0

static char ***dataSM;
static char **data;

void initDataSM( int pid )
{
	int i;
	int dataSMID;

	dataSMID = shmget( pid, 2, IPC_EXCL | IPC_CREAT | 0666 );
	if( dataSMID == -1 ) ErrorExit( "data sonna..\n" );

	dataSM = (char ***)shmat( dataSMID, 0, 0 );
	if( (int)dataSM == -1 ) ErrorExit( "Cannot attatch Shared Memory for Sequences!" );

	data = AllocateCharMtx( njob, nlenmax * 3 );
	dataSM[STATUS] = data;
	dataSM[SEMAPHORE] = (char **)1;
}

void writeDataSM( int nseq, char name[][B], int nlen[M], char **aseq, int force )
{
	int i;

	for( i=0; i<10; i++ )
	{
		if( dataSM[SEMAPHORE]-- > 0 )
		{
			signalSM[STATUS] = IMA_KAITERU;
			for( i=0; i<nseq; i++ ) strcpy( data[i], aseq[i] );
			dataSM[STATUS] = 100;
			if( force ) signalSM[STATUS] = OSHIMAI;
			else        signalSM[STATUS] = KAKIOWATTA;
			dataSM[SEMAPHORE]++;
			break;
		}
		else
		{
			dataSM[SEMAPHORE]++;
			if( !force ) break;
			else         sleep( 1 );
		}
	}
}
#endif

void WriteForFasta( FILE *fp, int locnjob, char name[][B], int nlen[M], char **aseq )
{
    static char b[N];
    int i, j;
    int nalen[M];

    for( i=0; i<locnjob; i++ )
    {
        nalen[i] = strlen( aseq[i] );
        fprintf( fp, ">%s\n", name[i] );
        for( j=0; j<nalen[i]; j=j+C ) 
        {
            strncpy( b, aseq[i]+j, C ); b[C] = 0;
            fprintf( fp, "%s\n",b );
        }
    }
}
