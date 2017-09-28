#include "fft.h"
#include "mtxutl.h"

/***********************************************************
        fft.c -- FFT (��®Fourier�Ѵ�)
***********************************************************/
/*
  �ؿ�{\tt fft()}�β������Ȥ��ƻ��Ѵؿ�ɽ����.
*/
static void make_sintbl(int n, float sintbl[])
{
        int i, n2, n4, n8;
        double c, s, dc, ds, t;

        n2 = n / 2;  n4 = n / 4;  n8 = n / 8;
#if 0
        t = sin(PI / n);
        dc = 2 * t * t;  ds = sqrt(dc * (2 - dc));
        t = 2 * dc;  c = sintbl[n4] = 1;  s = sintbl[0] = 0;
        for (i = 1; i < n8; i++) {
                c -= dc;  dc += t * c;
                s += ds;  ds -= t * s;
                sintbl[i] = s;  sintbl[n4 - i] = c;
        }
        if (n8 != 0) sintbl[n8] = sqrt(0.5);
#else
        sintbl[n4] = 1.0;  
        sintbl[0] = 0.0;
        for (i = 1; i < n4; i++)
            sintbl[i] = sin( (double)i / n * ( 2 * PI ) );
#endif
        for (i = 0; i < n4; i++)
                sintbl[n2 - i] = sintbl[i];
        for (i = 0; i < n2 + n4; i++)
                sintbl[i + n2] = - sintbl[i];
}
/*
  �ؿ�{\tt fft()}�β������Ȥ��ƥӥå�ȿžɽ����.
*/
static void make_bitrev(int n, int bitrev[])
{
        int i, j, k, n2;

        n2 = n / 2;  i = j = 0;
        for ( ; ; ) {
                bitrev[i] = j;
                if (++i >= n) break;
                k = n2;
                while (k <= j) {  j -= k;  k /= 2;  }
                j += k;
        }
}
/*
  ��®Fourier�Ѵ� (Cooley--Tukey�Υ��르�ꥺ��).
  ɸ�����ο� {\tt n} ��2��������˸¤�.
  {\tt x[$k$]} ������, {\tt y[$k$]} ������ ($k = 0$, $1$, $2$,
  \ldots, $|{\tt n}| - 1$).
  ��̤� {\tt x[]}, {\tt y[]} �˾�񤭤����.
  ${\tt n} = 0$ �ʤ�ɽ�Υ�����������.
  ${\tt n} < 0$ �ʤ���Ѵ���Ԥ�.
  ����Ȱۤʤ� $|{\tt n}|$ ���ͤǸƤӽФ���,
  ���Ѵؿ��ȥӥå�ȿž��ɽ���뤿���¿��;ʬ�˻��֤�������.
  ����ɽ�Τ���ε����ΰ�����˼��Ԥ����1���֤� (���ｪλ��
  ������ͤ�0).
  ������ɽ�ε����ΰ���������ˤ� ${\tt n} = 0$ �Ȥ���
  �ƤӽФ� (���ΤȤ��� {\tt x[]}, {\tt y[]} ���ͤ��Ѥ��ʤ�).
*/
int fft_hontai(int n, float x[], float y[])
{
        static int    last_n = 0;    /* ����ƽФ����� {\tt n} */
        static int   *bitrev = NULL; /* �ӥå�ȿžɽ */
        static float *sintbl = NULL; /* ���Ѵؿ�ɽ */
        int i, j, k, ik, h, d, k2, n4, inverse;
        float t, s, c, dx, dy;

        /* ���� */
        if (n < 0) {
                n = -n;  inverse = 1;  /* ���Ѵ� */
        } else inverse = 0;
        n4 = n / 4;
        if (n != last_n || n == 0) {
                last_n = n;
                if (sintbl != NULL) free(sintbl);
                if (bitrev != NULL) free(bitrev);
                if (n == 0) return 0;  /* �����ΰ��������� */
                sintbl = malloc((n + n4) * sizeof(float));
                bitrev = malloc(n * sizeof(int));
                if (sintbl == NULL || bitrev == NULL) {
                        fprintf(stderr, "�����ΰ���­\n");  return 1;
                }
                make_sintbl(n, sintbl);
                make_bitrev(n, bitrev);
        }
        for (i = 0; i < n; i++) {    /* �ӥå�ȿž */
                j = bitrev[i];
                if (i < j) {
                        t = x[i];  x[i] = x[j];  x[j] = t;
                        t = y[i];  y[i] = y[j];  y[j] = t;
                }
        }
        for (k = 1; k < n; k = k2) {    /* �Ѵ� */
                h = 0;  k2 = k + k;  d = n / k2;
                for (j = 0; j < k; j++) {
                        c = sintbl[h + n4];
                        if (inverse) s = - sintbl[h];
                        else         s =   sintbl[h];
                        for (i = j; i < n; i += k2) {
                                ik = i + k;
                                dx = s * y[ik] + c * x[ik];
                                dy = c * y[ik] - s * x[ik];
                                x[ik] = x[i] - dx;  x[i] += dx;
                                y[ik] = y[i] - dy;  y[i] += dy;
                        }
                        h += d;
                }
        }
        if (! inverse)    /* ���Ѵ��Ǥʤ��ʤ�n�ǳ�� */
                for (i = 0; i < n; i++) {  x[i] /= n;  y[i] /= n;  }
        return 0;  /* ���ｪλ */
}

#if 0
int fft( int n, Fukusosuu *in, int disp )
{
	int i;
	static int last_n = 0;
	static float *x, *y;
	if( last_n != abs( n ) )
	{
		
		if( last_n ) 
		{
			free( x );
			free( y );
		}
		last_n = abs( n );
		x = calloc( n, sizeof( float ) );
		y = calloc( n, sizeof( float ) );
	}
	for( i=0; i<n; i++ )
	{
		x[i] = (float)in[i].R;
		y[i] = (float)in[i].I;
	}	
	fft_hontai( n, x, y );
	for( i=0; i<n; i++ )
	{
		(in+i)->R = (double)x[i];
		(in+i)->I = (double)y[i];
	}	
#if 1
if( disp )                                                      
{                                                                    
    FILE *fp;        
    fp = fopen( "outputOfFft", "w" );                                           
    for( i=0; i<n; i++ )                
        fprintf( fp, "%f %f\n", in[i].R, in[i].I );        
    fclose( fp );        
    system( "vi outputOfFft < /dev/tty > /dev/tty " );        
}                
#endif        
}
#else
int fft( int n, Fukusosuu *in, int disp )
{
	int i, m;
	float *x, *y;

	m = abs( n );

	x = calloc( m, sizeof( float ) );
	y = calloc( m, sizeof( float ) );

	for( i=0; i<m; i++ )
	{
		x[i] = (float)in[i].R;
		y[i] = (float)in[i].I;
	}	

#if 0
if( disp )                                                      
{                                                                    
    FILE *fp;        
    fp = fopen( "inputOfFft", "w" );                                           
	fprintf( fp, "m=%d\n", m );
    for( i=0; i<m; i++ )                
        fprintf( fp, "%f %f\n", x[i], y[i] );        
    fclose( fp );        
    system( "vi inputOfFft < /dev/tty > /dev/tty " );        
}                
#endif        
	fft_hontai( n, x, y );
#if 0
if( disp )                                                      
{                                                                    
    FILE *fp;        
    fp = fopen( "outputOfFft", "w" );                                           
	fprintf( fp, "m=%d\n", m );
    for( i=0; i<m; i++ )                
        fprintf( fp, "%f %f\n", x[i], y[i] );        
    fclose( fp );        
    system( "vi outputOfFft < /dev/tty > /dev/tty " );        
}                
#endif        
	for( i=0; i<m; i++ )
	{
		(in+i)->R = (double)x[i];
		(in+i)->I = (double)y[i];
	}	
	free( x );
	free( y );
}
#endif
