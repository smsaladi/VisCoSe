#include "mltaln.h"
#include "dp.h"

#define OUTGAP0TRY 1
#define DEBUG 0
#define XXXXXXX    0
#define USE_PENALTY_EX  0

static void OpeningGapCount( double *ogcp, int clus, char **seq, double *eff )
{
	int i, j, gc, gb; 
	int len = strlen( seq[0] );
	double totaleff = 0.0;
	
	for( i=0; i<len; i++ ) ogcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		gc = 0;
		for( i=0; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			{
				if( !gb *  gc ) ogcp[i] += eff[j];
			}
		}
		totaleff+= eff[j];
	}
	for( i=0; i<len; i++ ) 
		ogcp[i] /= totaleff;
}

static void FinalGapCount( double *fgcp, int clus, char **seq, double *eff )
{
	int i, j, gc, gb; 
	int len = strlen( seq[0] );
	double totaleff = 0.0;
	
	for( i=0; i<len; i++ ) fgcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		gc = ( seq[j][0] == '-' );
		for( i=1; i<len+1; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			{
				if( gb * !gc ) fgcp[i-1] += eff[j];
			}
		}
		totaleff += eff[j];
	}
	for( i=0; i<len; i++ ) 
		fgcp[i] /= totaleff;
}


			
		
#if 0
		
static void v_calc11( char *seq1, char *seq2, float **v, int lgth1, int lgth2 )
{
	int i, j;

	for( i=0; i<lgth1; i++ )
		for( j=0; j<lgth2; j++ )
			v[i][j] = amino_dis[seq1[i]][seq2[j]];
}

static void v_calc1m( char *seq1, float **cpmx2, float **v, int lgth1, int lgth2 )
{
	int i, j, k, l;
	float scarr[26];

	for( j=0; j<lgth2; j++ )
	{
		for( l=0; l<26; l++ )
		{
			scarr[l] = 0.0;
			for( k=0; k<26; k++ )
				scarr[l] += n_dis[k][l] * cpmx2[k][j];
		}
		for( i=0; i<lgth1; i++ )
		{
			v[i][j] = scarr[amino_n[seq1[i]]];
		}
	}
}

static void v_calcm1( float **cpmx1, char *seq2, float **v, int lgth1, int lgth2 )
{
	int i, j, k, l;
	float scarr[26];

	for( i=0; i<lgth1; i++ )
	{
		for( l=0; l<26; l++ )
		{
			scarr[l] = 0.0;
			for( k=0; k<26; k++ )
				scarr[l] += n_dis[k][l] * cpmx1[k][i];
		}
		for( j=0; j<lgth2; j++ )
		{
			v[i][j] = scarr[amino_n[seq2[j]]];
		}
	}
}
			
static void v_calcmm( float **cpmx1, float **cpmx2, float **v, int lgth1, int lgth2 )
{
	int i, j, k, l;
	float scarr[26];
	float cpmxpd[26][N];
	int cpmxpdn[26][N];
	int count = 0;
	for( j=0; j<lgth2; j++ )
	{
		count = 0;
		for( l=0; l<26; l++ )
		{
			if( cpmx2[l][j] )
			{
				cpmxpd[count][j] = cpmx2[l][j];
				cpmxpdn[count][j] = l;
				count++;
			}
		}
		cpmxpdn[count][j] = -1;
	}
	for( i=0; i<lgth1; i++ )
	{
		for( l=0; l<26; l++ )
		{
			scarr[l] = 0.0;
			for( k=0; k<26; k++ )
				scarr[l] += n_dis[k][l] * cpmx1[k][i];
		}
		for( j=0; j<lgth2; j++ )
		{
			v[i][j] = 0;
			for( k=0; cpmxpdn[k][j] > -1;  k++ )
				v[i][j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
		} 
	}
}
#endif

static void match_calc( float *match, float **cpmx1, float **cpmx2, int i1, int lgth2, float **floatwork, int **intwork, int initialize )
{
	int j, k, l;
	float scarr[26];
	float **cpmxpd = floatwork;
	int **cpmxpdn = intwork;
	int count = 0;

	if( initialize )
	{
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<26; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[count][j] = cpmx2[l][j];
					cpmxpdn[count][j] = l;
					count++;
				}
			}
			cpmxpdn[count][j] = -1;
		}
	}

	for( l=0; l<26; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<26; k++ )
			scarr[l] += n_dis[k][l] * cpmx1[k][i1];
	}
#if 0 /* �����Ȥ��Ȥ���floatwork�Υ������Ȥ�դˤ��� */
	{
		float *fpt, **fptpt, *fpt2;
		int *ipt, **iptpt;
		fpt2 = match;
		iptpt = cpmxpdn;
		fptpt = cpmxpd;
		while( lgth2-- )
		{
			*fpt2 = 0.0;
			ipt=*iptpt,fpt=*fptpt;
			while( *ipt > -1 )
				*fpt2 += scarr[*ipt++] * *fpt++;
			fpt2++,iptpt++,fptpt++;
		} 
	}
#else
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[k][j]>-1; k++ )
			match[j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
	} 
#endif
}

static float Atracking( float *lasthorizontalw, float *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        float **cpmx1, float **cpmx2, 
                        short **ijp, int icyc, int jcyc )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k;
	char gap[] = "-";
	float wm;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
	if( outgap == 1 )
		;
	else
	{
		wm = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wm )
			{
				wm = lastverticalw[i];
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( lasthorizontalw[j] >= wm )
			{
				wm = lasthorizontalw[j];
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
	}

    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }

	for( i=0; i<icyc; i++ )
	{
		mseq1[i] += lgth1+lgth2;
		*mseq1[i] = 0;
	}
	for( j=0; j<jcyc; j++ )
	{
		mseq2[j] += lgth1+lgth2;
		*mseq2[j] = 0;
	}
	iin = lgth1; jin = lgth2;
	for( k=0; k<=lgth1+lgth2; k++ ) 
	{
		if( ijp[iin][jin] < 0 ) 
		{
			ifi = iin-1; jfi = jin+ijp[iin][jin];
		}
		else if( ijp[iin][jin] > 0 )
		{
			ifi = iin-ijp[iin][jin]; jfi = jin-1;
		}
		else
		{
			ifi = iin-1; jfi = jin-1;
		}
		l = iin - ifi;
		while( --l ) 
		{
			for( i=0; i<icyc; i++ )
				*--mseq1[i] = seq1[i][ifi+l];
			for( j=0; j<jcyc; j++ ) 
				*--mseq2[j] = *gap;
			k++;
		}
		l= jin - jfi;
		while( --l )
		{
			for( i=0; i<icyc; i++ ) 
				*--mseq1[i] = *gap;
			for( j=0; j<jcyc; j++ ) 
				*--mseq2[j] = seq2[j][jfi+l];
			k++;
		}
		if( iin <= 0 || jin <= 0 ) break;
		for( i=0; i<icyc; i++ ) 
			*--mseq1[i] = seq1[i][ifi];
		for( j=0; j<jcyc; j++ ) 
			*--mseq2[j] = seq2[j][jfi];
		k++;
		iin = ifi; jin = jfi;
	}
	return( 0.0 );
}


float A__align( char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	float wm;   /* int ?????? */
	float g;
	float *currentw, *previousw;
#if 1
	float *wtmp;
	short *ijppt;
	float *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static float mi, *m;
	static short **ijp;
	static int mpi, *mp;
	static float *w1, *w2;
	static float *match;
	static float *initverticalw;    /* kufuu sureba iranai */
	static float *lastverticalw;    /* kufuu sureba iranai */
	static char **mseq1;
	static char **mseq2;
	static char **mseq;
	static double *ogcp1;
	static double *ogcp2;
	static double *fgcp1;
	static double *fgcp2;
	static float **cpmx1;
	static float **cpmx2;
	static int **intwork;
	static float **floatwork;
	static int orlgth1 = 0, orlgth2 = 0;

#if 0
	fprintf( stderr, "eff in SA+++align\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
#endif
	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( njob, 1 );
		mseq2 = AllocateCharMtx( njob, 1 );
	}

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeDoubleVec( ogcp1 );
			FreeDoubleVec( ogcp2 );
			FreeDoubleVec( fgcp1 );
			FreeDoubleVec( fgcp2 );

			FreeFloatMtx( cpmx1 );
			FreeFloatMtx( cpmx2 );

			FreeFloatMtx( floatwork );
			FreeIntMtx( intwork );
		}

		ll1 = MAX( (int)(1.1*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.1*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( njob, ll1+ll2 );

		ogcp1 = AllocateDoubleVec( ll1+2 );
		ogcp2 = AllocateDoubleVec( ll2+2 );
		fgcp1 = AllocateDoubleVec( ll1+2 );
		fgcp2 = AllocateDoubleVec( ll2+2 );

		cpmx1 = AllocateFloatMtx( 26, ll1+2 );
		cpmx2 = AllocateFloatMtx( 26, ll2+2 );

		floatwork = AllocateFloatMtx( 26, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( 26, MAX( ll1, ll2 )+2 ); 

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1;
		orlgth2 = ll2;
	}

	for( i=0; i<icyc; i++ ) mseq1[i] = mseq[i];
	for( j=0; j<jcyc; j++ ) mseq2[j] = mseq[icyc+j];


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 )
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 )
		{
			FreeShortMtx( commonIP );
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

#if DEBUG
		fprintf( stderr, "\n\ntrying to allocate %dx%d common matrices ... ", ll1+1, ll2+1 );
#endif

		commonIP = AllocateShortMtx( ll1+10, ll2+10 );

#if DEBUG
		fprintf( stderr, "succeeded\n\n" );
#endif

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;

	cpmx_calc( seq1, cpmx1, eff1, strlen( seq1[0] ), icyc );
	cpmx_calc( seq2, cpmx2, eff2, strlen( seq2[0] ), jcyc );

	OpeningGapCount( ogcp1, icyc, seq1, eff1 );
	OpeningGapCount( ogcp2, jcyc, seq2, eff2 );
	FinalGapCount( fgcp1, icyc, seq1, eff1 );
	FinalGapCount( fgcp2, jcyc, seq2, eff2 );

	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 1.0 - ogcp1[i];
		fgcp1[i] = 1.0 - fgcp1[i];
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 1.0 - ogcp2[i];
		fgcp2[i] = 1.0 - fgcp2[i];
	}
#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;

	match_calc( initverticalw, cpmx2, cpmx1, 0, lgth1, floatwork, intwork, 1 );
	match_calc( currentw, cpmx1, cpmx2, 0, lgth2, floatwork, intwork, 1 );


	if( outgap == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += penalty * 0.5 * ( ogcp1[0] + fgcp1[i-1] ) ;
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += penalty * 0.5 * ( ogcp2[0] + fgcp2[j-1] ) ;
		}
	}
#if OUTGAP0TRY
	else
	{
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * j / 2.0;
		for( i=1; i<lgth1+1; i++ )
			initverticalw[i] -= offset * i / 2.0;
	}
#endif

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + penalty * 0.5 * ogcp1[1]; mp[j] = 0;
	}

	lastverticalw[0] = currentw[lgth2-1];

	if( outgap ) lasti = lgth1+1; else lasti = lgth1;

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "fcgp\n" );
for( i=0; i<lgth1; i++ ) 
	fprintf( stderr, "fgcp1[%d]=%f\n", i, ogcp1[i] );
for( i=0; i<lgth2; i++ ) 
	fprintf( stderr, "fgcp2[%d]=%f\n", i, ogcp2[i] );
#endif

	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc( currentw, cpmx1, cpmx2, i, lgth2, floatwork, intwork, 0 );
		currentw[0] = initverticalw[i];

#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth1+1; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if 1

		mi = previousw[0] + penalty * 0.5 * ogcp2[1]; mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		for( j=1; j<lgth2+1; j++ )
		{
			wm = *prept;
			*ijppt = 0;

#if XXXXXXX
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + (float)penalty * fgcp2[j-1] * 0.5;
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
			}
			g = *prept + penalty * ogcp2[j] * 0.5;
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += penalty_ex;
#endif

			g = *mjpt + (float)penalty * fgcp1[i-1] * 0.5;
#if XXXXXXX 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			g = *prept + penalty * ogcp1[i] * 0.5;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += penalty_ex;
#endif

#if XXXXXXX
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;
		}
		lastverticalw[i] = currentw[lgth2-1];
#else
        floatncpy( previousw, currentw, lgth2+1 );
        previousw[0] = initverticalw[i-1];

        match_calc( currentw, cpmx1, cpmx2, i, lgth2, floatwork, intwork, 0 );
        currentw[0] = initverticalw[i];

        mi = previousw[0] + penalty * 0.5 * ogcp2[1]; mpi = 0;
        for( j=1; j<lgth2+1; j++ )
        {
            wm = previousw[j-1];
            ijp[i][j] = 0;

            g = penalty * fgcp2[j-1] * 0.5;
            if( mi + g > wm )
            {
                wm = mi + g; 
                ijp[i][j] = -( j - mpi );
            }
            g = penalty * ogcp2[j] * 0.5;
            if( mi <= previousw[j-1] + g )
            {
                mi = previousw[j-1] + g;
                mpi = j-1;
            }

            g = penalty * fgcp1[i-1] * 0.5;
            if( m[j] + g > wm )
            {
                wm = m[j] + g;
                ijp[i][j] = +( i - mp[j] );
            }
            g = penalty * ogcp1[i] * 0.5;
            if( m[j] <= previousw[j-1] + g )
            {
                m[j] = previousw[j-1] + g ;
                mp[j] = i-1;
            }
            currentw[j] += wm;
        }
        lastverticalw[i] = currentw[lgth2-1];

#endif
	}

#if OUTGAP0TRY
	if( !outgap )
	{
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * ( lgth2 - j ) / 2.0;
		for( i=1; i<lgth1+1; i++ )
			lastverticalw[i] -= offset * ( lgth1 - i  / 2.0);
	}
#endif
		
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr,"%s\n", seq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr,"%s\n", seq2[j] );
	fprintf( stderr, "====>" );
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], seq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( mseq2[j], seq2[j] );
	*/
	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, cpmx1, cpmx2, ijp, icyc, jcyc );

	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( seq2[j], mseq2[j] );
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
	*/
	return( wm );
}

float translate_and_Calign( char **mseq1, char **mseq2, double *effarr1, double *effarr2, int clus1, int clus2, int alloclen )
{
    int i;
    float wm;
    char **result;
    char *seq, **aseq;
    double *effarr;
    int nseq;
	int resultlen;

    if     ( clus1 == 1 && clus2 != 1 ) 
    {
        seq = mseq1[0]; aseq = mseq2; effarr = effarr2; nseq = clus2+1;
#if 0
		printf( "effarr in transl... = \n" );
		for( i=0; i<clus2; i++ ) printf( "%f ", effarr2[i] );
#endif
    }
    else if( clus1 != 1 && clus2 == 1 ) 
    {
        seq = mseq2[0]; aseq = mseq1; effarr = effarr1; nseq = clus1+1;
    }
    else ErrorExit( "ERROR in translate_and_Calign" );

    result = Calignm1( &wm, aseq, seq, effarr, nseq-2, 0 );

	resultlen = strlen( result[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}
    for( i=0; i<nseq-1; i++ ) strcpy( aseq[i], result[i] );
    strcpy( seq, result[nseq-1] );

    return( 0.0 );
}
