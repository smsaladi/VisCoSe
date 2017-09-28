#include "mltaln.h"
#include "miyata.h"
#include "miyata5.h"
#include "gonnet.h"
#include "DNA.h"

#include "JTT.c"

#define DEBUG 0
#define TEST 0

#define NORMALIZE1 1

int shishagonyuu( double in )
{
	int out;
	if     ( in >  0.0 ) out = ( (int)( in + 0.5 ) );
	else if( in == 0.0 ) out = ( 0 );
	else if( in <  0.0 ) out = ( (int)( in - 0.5 ) );
	else                 out = 0;
	return( out );
}

void constants()
{
	int i, j, x;
	double tmp;


	if( scoremtx == 1 )  /* Gonnet */
	{
		for( i=0; i<26; i++ ) for( j=0; j<26; j++ ) n_dis[i][j] = locn_disd[i][j];
		for( i=0; i<26; i++ ) if( i != 24 ) n_dis[i][24] = n_dis[24][i] = 0.0;
		n_dis[24][24] = 0;
		if( ppenalty == NOTSPECIFIED ) ppenalty = locpenaltyd;
		if( poffset == NOTSPECIFIED ) poffset = locexpenaltyd;
		if( pamN == NOTSPECIFIED ) pamN = 0;
		if( kimuraR == NOTSPECIFIED ) kimuraR = 1;
		penalty = ppenalty;
		offset = poffset;

		for( i=0; i<26; i++ ) amino[i] = locaminod[i];
		for( i=0; i<26; i++ ) amino_grp[amino[i]] = locgrpd[i];
		for( i=0; i<20; i++ ) for( j=0; j<20; j++ )
			n_dis[i][j] -= offset;
#if 1
		fprintf( stderr, "scoreing matrix : \n" );
		for( i=0; i<26; i++ )
		{
			for( j=0; j<26; j++ ) 
			{
				fprintf( stdout, "%#5d", n_dis[i][j] );
			}
			fprintf( stdout, "\n" );
		}
#endif
	}
	else if( scoremtx == 2 ) /* Miyata-Yasunaga */
	{
		for( i=0; i<26; i++ ) for( j=0; j<26; j++ ) n_dis[i][j] = locn_dism[i][j];
		for( i=0; i<26; i++ ) if( i != 24 ) n_dis[i][24] = n_dis[24][i] = exgpm;
		n_dis[24][24] = 0;
		if( ppenalty == NOTSPECIFIED ) ppenalty = locpenaltym;
		if( poffset == NOTSPECIFIED ) poffset = -20;
		if( pamN == NOTSPECIFIED ) pamN = 0;
		if( kimuraR == NOTSPECIFIED ) kimuraR = 1;

		penalty = ppenalty;
		offset = poffset;

		for( i=0; i<26; i++ ) amino[i] = locaminom[i];
		for( i=0; i<26; i++ ) amino_grp[amino[i]] = locgrpm[i];
#if DEBUG
		fprintf( stdout, "scoreing matrix : \n" );
		for( i=0; i<26; i++ )
		{
			for( j=0; j<26; j++ ) 
			{
				fprintf( stdout, "%#5d", n_dis[i][j] );
			}
			fprintf( stdout, "\n" );
		}
#endif
	}
	else if( scoremtx == -1 )  /* DNA */
	{
		double average;
		double **pamx = AllocateDoubleMtx( 11,11 );
		double **pam1 = AllocateDoubleMtx( 4, 4 );

		if( ppenalty == NOTSPECIFIED ) ppenalty = DEFAULTGOP_N;
		if( ppenalty_ex == NOTSPECIFIED ) ppenalty_ex = DEFAULTGEP_N;
		if( poffset == NOTSPECIFIED ) poffset = DEFAULTOFS_N;
		if( pamN == NOTSPECIFIED ) pamN = DEFAULTPAMN;
		if( kimuraR == NOTSPECIFIED ) kimuraR = 2;

		penalty = (int)( 600.0 / 1000.0 * ppenalty + 0.5);
		penalty_ex = (int)( 600.0 / 1000.0 * ppenalty_ex + 0.5);
		offset = (int)( 600.0 / 1000.0 * poffset + 0.5);

		if( kimuraR == 9999 ) 
		{
			for( i=0; i<4; i++ ) for( j=0; j<4; j++ ) 
				pamx[i][j] = (double)locn_disn[i][j];
#if NORMALIZE1
			average = 0.0;
			for( i=0; i<4; i++ ) for( j=0; j<4; j++ ) 
				average += pamx[i][j];
			average /= 16.0;
	
   	     if( disp )
				fprintf( stderr, "average = %f\n", average );
	
			for( i=0; i<4; i++ ) for( j=0; j<4; j++ ) 
				pamx[i][j] -= average;
	
			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
				pamx[i][j] *= 600.0 / average;
			
			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
				pamx[i][j] -= offset; 
#endif
		}
		else
		{
			double f = 0.99;
			double s = (double)kimuraR / ( 2 + kimuraR ) * 0.01;
			double v = (double)1       / ( 2 + kimuraR ) * 0.01;
			pam1[0][0] = f; pam1[0][1] = s; pam1[0][2] = v; pam1[0][3] = v;
			pam1[1][0] = s; pam1[1][1] = f; pam1[1][2] = v; pam1[1][3] = v;
			pam1[2][0] = v; pam1[2][1] = v; pam1[2][2] = f; pam1[2][3] = s;
			pam1[3][0] = v; pam1[3][1] = v; pam1[3][2] = s; pam1[3][3] = f;

			fprintf( stderr, "generating %dPAM scoring matrix for nucleotides ... ", pamN );

	       	if( disp )
   	    	{
   	     		fprintf( stderr, " TPM \n" );
   	        	for( i=0; i<4; i++ )
   		       	{
   	            	for( j=0; j<4; j++ )
   	                	fprintf( stderr, "%+#6.10f", pam1[i][j] );
   	            	fprintf( stderr, "\n" );
   	        	}
   	        	fprintf( stderr, "\n" );
   	     	}


			MtxuntDouble( pamx, 4 );
			for( x=0; x < pamN; x++ ) MtxmltDouble( pamx, pam1, 4 );
			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
				pamx[i][j] /= 1.0 / 4.0;

			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
			{
				if( pamx[i][j] == 0.0 ) 
				{
					fprintf( stderr, "WARNING: pamx[i][j] = 0.0 ?\n" );
					pamx[i][j] = 0.00001; /* by J. Thompson */
				}
				pamx[i][j] = log10( pamx[i][j] ) * 1000.0;
			}

       		if( disp )
       		{
        		fprintf( stderr, " after log\n" );
           		for( i=0; i<4; i++ )
   	       		{
           	    	for( j=0; j<4; j++ )
           	        	fprintf( stderr, "%+#6.10f", pamx[i][j] );
           	    	fprintf( stderr, "\n" );
           		}
           		fprintf( stderr, "\n" );
        	}
			
			average = 0.0;
			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
				average += pamx[i][j] * ( 1.0 / 4.0 ) * ( 1.0 / 4.0 );
			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
				pamx[i][j] -= average;

			average = 0.0;
			for( i=0; i<4; i++ )
				average += pamx[i][i] * 1.0 / 4.0;

			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
				pamx[i][j] *= 600.0 / average;


			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
				pamx[i][j] -= offset;        /* extending gap cost */

			for( i=0; i<4; i++ ) for( j=0; j<4; j++ )
				pamx[i][j] = shishagonyuu( pamx[i][j] );

       		if( disp )
       		{
        		fprintf( stderr, " after shishagonyuu\n" );
           		for( i=0; i<4; i++ )
   	       		{
           	    	for( j=0; j<4; j++ )
           	        	fprintf( stderr, "%+#6.10f", pamx[i][j] );
           	    	fprintf( stderr, "\n" );
           		}
           		fprintf( stderr, "\n" );
        	}
			fprintf( stderr, "done\n" );
		}
	
		for( i=0; i<5; i++ ) 
		{
			pamx[4][i] = pamx[1][i];
			pamx[i][4] = pamx[i][1];
		}	

		for( i=5; i<10; i++ ) for( j=5; j<10; j++ )
			pamx[i][j] = pamx[i-5][j-5];
	
       	if( disp )
       	{
       		fprintf( stderr, " before dis\n" );
          	for( i=0; i<4; i++ )
   	       	{
           	   	for( j=0; j<4; j++ )
           	       	fprintf( stderr, "%+#6.10f", pamx[i][j] );
           	   	fprintf( stderr, "\n" );
           	}
           	fprintf( stderr, "\n" );
        }

       	if( disp )
       	{
        	fprintf( stderr, " score matrix  \n" );
           	for( i=0; i<4; i++ )
   	       	{
               	for( j=0; j<4; j++ )
                   	fprintf( stderr, "%+#6.10f", pamx[i][j] );
               	fprintf( stderr, "\n" );
           	}
           	fprintf( stderr, "\n" );
        }

		for( i=0; i<26; i++ ) amino[i] = locaminon[i];
		for( i=0; i<26; i++ ) amino_grp[amino[i]] = locgrpn[i];
		for( i=0; i<26; i++ ) for( j=0; j<26; j++ ) n_dis[i][j] = 0;
		for( i=0; i<10; i++ ) for( j=0; j<10; j++ ) n_dis[i][j] = shishagonyuu( pamx[i][j] );
        if( disp )
        {
            fprintf( stderr, " score matrix  \n" );
            for( i=0; i<26; i++ )
            {
                for( j=0; j<26; j++ )
                    fprintf( stderr, "%+#6d", n_dis[i][j] );
                fprintf( stderr, "\n" );
            }
            fprintf( stderr, "\n" );
        }
	}
	else         /* JTT */
	{
		double **pam1;
		double **pamx;
		double *freq;
		double average;

		pam1 = AllocateDoubleMtx( 20, 20 );
		pamx = AllocateDoubleMtx( 20, 20 );
		freq = AllocateDoubleVec( 20 );

		if( ppenalty == NOTSPECIFIED ) ppenalty = DEFAULTGOP_J;
		if( ppenalty_ex == NOTSPECIFIED ) ppenalty_ex = DEFAULTGEP_J;
		if( poffset == NOTSPECIFIED ) poffset = DEFAULTOFS_J;
		if( pamN == NOTSPECIFIED )    pamN    = DEFAULTPAMN;
		if( kimuraR == NOTSPECIFIED ) kimuraR = 1;

		penalty = (int)( 600.0 / 1000.0 * ppenalty + 0.5 );
		penalty_ex = (int)( 600.0 / 1000.0 * ppenalty_ex + 0.5 );
		offset = (int)( 600.0 / 1000.0 * poffset + 0.5 );

		JTTmtx( pam1, freq, amino, amino_grp );

		fprintf( stderr, "generating %dPAM scoring matrix for amino acids ... ", pamN );

		MtxuntDouble( pamx, 20 );
		for( x=0; x < pamN; x++ ) MtxmltDouble( pamx, pam1, 20 );

		for( i=0; i<20; i++ ) for( j=0; j<20; j++ ) 
			pamx[i][j] /= freq[i];
        for( i=0; i<20; i++ ) for( j=0; j<20; j++ )
		{
			if( pamx[i][j] == 0.0 ) 
			{
				fprintf( stderr, "WARNING: pamx[%d][%d] = 0.0?\n", i, j );
				pamx[i][j] = 0.00001; /* by J. Thompson */
			}
            pamx[i][j] = log10( pamx[i][j] ) * 1000.0;
		}
 
#if TEST
        average = tmp = 0.0;
        for( i=0; i<20; i++ ) for( j=0; j<20; j++ )
		{
           average += pamx[i][j] * freq[i] * freq[j];
		   tmp += freq[i] * freq[j];
		}
		average /= tmp;
		fprintf( stderr, "Zenbu average = %f, tmp = %f \n", average, tmp );
        average = tmp = 0.0;
        for( i=0; i<20; i++ ) for( j=i; j<20; j++ )
		{
           average += pamx[i][j] * freq[i] * freq[j];
		   tmp += freq[i] * freq[j];
		}
		average /= tmp;
		fprintf( stderr, "Zenbu average2 = %f, tmp = %f \n", average, tmp );
		average = tmp = 0.0;
		for( i=0; i<20; i++ )
		{
			average += pamx[i][i] * freq[i];
			tmp += freq[i];
		}
		average /= tmp;
		fprintf( stderr, "Itch average = %f, tmp = %f \n", average, tmp );
		fprintf( stdout, "raw scoring matrix : \n" );
		for( i=0; i<20; i++ )
		{
			for( j=0; j<=i; j++ ) 
			{
				fprintf( stdout, "%5.0f", pamx[i][j] );
			}
			fprintf( stdout, "\n" );
		}
#endif

#if NORMALIZE1
		average = 0.0;
		for( i=0; i<20; i++ ) for( j=0; j<20; j++ )
			average += pamx[i][j] * freq[i] * freq[j];

		for( i=0; i<20; i++ ) for( j=0; j<20; j++ ) 
			pamx[i][j] -= average;
#if TEST
		fprintf( stdout, "after average substruction : \n" );
		for( i=0; i<20; i++ )
		{
			for( j=0; j<=i; j++ ) 
			{
				fprintf( stdout, "%5.0f", pamx[i][j] );
			}
			fprintf( stdout, "\n" );
		}
#endif
		
		average = 0.0;
		for( i=0; i<20; i++ ) 
			average += pamx[i][i] * freq[i];

		for( i=0; i<20; i++ ) for( j=0; j<20; j++ ) 
			pamx[i][j] *= 600.0 / average;
#if TEST
        fprintf( stdout, "after average division : \n" );
        for( i=0; i<20; i++ )
        {
            for( j=0; j<=i; j++ )
            {
                fprintf( stdout, "%5.0f", pamx[i][j] );
            }
            fprintf( stdout, "\n" );
        }
#endif

		for( i=0; i<20; i++ ) for( j=0; j<20; j++ ) 
			pamx[i][j] -= offset;  
#if TEST
		fprintf( stdout, "after offset substruction (offset = %d): \n", offset );
		for( i=0; i<20; i++ )
		{
			for( j=0; j<=i; j++ ) 
			{
				fprintf( stdout, "%5.0f", pamx[i][j] );
			}
			fprintf( stdout, "\n" );
		}
#endif
#if 0
/* 注意 ！！！！！！！！！！ */
			penalty -= offset;
#endif


		for( i=0; i<20; i++ ) for( j=0; j<20; j++ ) 
			pamx[i][j] = shishagonyuu( pamx[i][j] );

#else

        average = 0.0;
        for( i=0; i<20; i++ ) for( j=0; j<20; j++ )
           average += pamx[i][j];
        average /= 400.0;

        for( i=0; i<20; i++ ) for( j=0; j<20; j++ )
        {
            pamx[i][j] -= average;
            pamx[i][j] = shishagonyuu( pamx[i][j] );
        }
#endif
        if( disp )
        {
            fprintf( stdout, " scoring matrix  \n" );
            for( i=0; i<20; i++ )
            {
				fprintf( stdout, "%c    ", amino[i] );
                for( j=0; j<=20; j++ )
                    fprintf( stdout, "%5.0f", pamx[i][j] );
                fprintf( stdout, "\n" );
            }
			fprintf( stdout, "     " );
            for( i=0; i<20; i++ )
				fprintf( stdout, "    %c", amino[i] );

			average = 0.0;
        	for( i=0; i<20; i++ ) for( j=0; j<20; j++ )
				average += pamx[i][j] * freq[i] * freq[j];
			fprintf( stdout, "average = %f\n", average );

			average = 0.0;
        	for( i=0; i<20; i++ )
				average += pamx[i][i] * freq[i];
			fprintf( stdout, "itch average = %f\n", average );
			fprintf( stderr, "parameters: %d, %d, %d\n", penalty, penalty_ex, offset );

			
			exit( 1 );
        }

		for( i=0; i<26; i++ ) for( j=0; j<26; j++ ) n_dis[i][j] = 0;
		for( i=0; i<20; i++ ) for( j=0; j<20; j++ ) n_dis[i][j] = (int)pamx[i][j];
		for( i=0; i<0x80; i++ ) amino_n[i] = 0;
		for( i=0; i<26; i++ ) amino_n[amino[i]] = i;

		fprintf( stderr, "done.\n" );
	}

#if DEBUG
	fprintf( stderr, "scoremtx = %d\n", scoremtx );
	fprintf( stderr, "amino[] = %s\n", amino );
#endif

	for( i=0; i<0x80; i++ )amino_n[i] = -1;
	for( i=0; i<26; i++)
		amino_n[amino[i]] = i;
    for( i=0; i<0x80; i++ ) for( j=0; j<0x80; j++ ) amino_dis[i][j] = 0;
    for( i=0; i<26; i++) for( j=0; j<26; j++ )
	{
        amino_dis[amino[i]][amino[j]] = n_dis[i][j];
	}


	ppid = 0;


	if( fftThreshold == NOTSPECIFIED )
	{
		fftThreshold = FFT_THRESHOLD;
	}
	if( fftWinSize == NOTSPECIFIED )
	{
		if( scoremtx == -1 ) 
			fftWinSize = FFT_WINSIZE_D;
		else    
			fftWinSize = FFT_WINSIZE_P;
	}


	if( fftscore )
	{
		double av, sd;

		for( i=0; i<20; i++ ) polarity[i] = polarity_[i];
		for( av=0.0, i=0; i<20; i++ ) av += polarity[i];
		av /= 20.0;
		for( sd=0.0, i=0; i<20; i++ ) sd += ( polarity[i]-av ) * ( polarity[i]-av );
		sd /= 20.0; sd = sqrt( sd );
		for( i=0; i<20; i++ ) polarity[i] -= av;
		for( i=0; i<20; i++ ) polarity[i] /= sd;
	
		for( i=0; i<20; i++ ) volume[i] = volume_[i];
		for( av=0.0, i=0; i<20; i++ ) av += volume[i];
		av /= 20.0;
		for( sd=0.0, i=0; i<20; i++ ) sd += ( volume[i]-av ) * ( volume[i]-av );
		sd /= 20.0; sd = sqrt( sd );
		for( i=0; i<20; i++ ) volume[i] -= av;
		for( i=0; i<20; i++ ) volume[i] /= sd;

#if 0
		for( i=0; i<20; i++ ) fprintf( stdout, "amino=%c, pol = %f<-%f, vol = %f<-%f\n", amino[i], polarity[i], polarity_[i], volume[i], volume_[i] );
		for( i=0; i<20; i++ ) fprintf( stdout, "%c %+5.3f %+5.3f\n", amino[i], volume[i], polarity[i] );
#endif
	}
}
