#include "mltaln.h"

#define DEBUG 1 

int intcmp( int *str1, int *str2 )
{
	while( *str1 != -1 && *str2 != -1 )
		if( *str1++ != *str2++ ) return( 1 );
	if( *str1 != *str2 ) return( 1 );
	return( 0 );
}

char **arguments( int argc, char *argv[] )
{
    int c;
	
	calledByXced = 0;
	devide = 0;
	fftscore = 1;
	use_fft = 0;
	alg = 'C';
    weight = 0;
    utree = 1;
	tbutree = 0;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'C';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
    scoremtx = NOTSPECIFIED;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;


    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( c = *++argv[0] )
		{
            switch( c )
            {
				case 'Q':
					calledByXced = 1;
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'A':
					alg = 'A';
					break;
				case 'C':
					alg = 'C';
					break;
				case 'S':
					alg = 'S';
					break;
				case 'd':
					disp = 1;
					break;
				case 'o':
					outgap = 0;
					break;
				case 'j':
					tbutree = 1;
					break;
				case 'u':
					tbweight = 0;
					break;
				case 'm':
					tbweight = 3;
					break;
				case 'z':
					fftThreshold = atoi( *++argv );
					--argc;
					goto nextoption;
				case 'w':
					fftWinSize = atoi( *++argv );
					--argc;
					goto nextoption;
				case 'Z':
					checkC = 1;
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc != 2 ) 
    {
        fprintf( stderr, "options: Check source file ! %c ?\n", c );
        exit( 1 );
    }
	fprintf( stderr, "tbitr = %d, tbrweight = %d, tbweight = %d\n", tbitr, tbrweight, tbweight );
//	readOtherOptions( &ppid, &fftThreshold, &fftWinSize );
	return( argv ); 

}

void GroupAlign( int nseq1, int nseq2, char **name, int *nlen, char **seq, char **aseq, char **mseq1, char **mseq2, int ***topol, double **len, double **eff, int alloclen )
{
	int i, j, l;
	int clus1, clus2;
	int s1, s2;
	float pscore;
	FILE *trap;
	static char **name1, **name2;
	double *effarr = NULL;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	static char *indication1, *indication2;
	static int *checkarr; 


	fprintf( stderr, "in GroupAlign fftWinSize   = %d\n", fftWinSize );
	fprintf( stderr, "in GroupAlign fftThreshold = %d\n", fftThreshold );

	if( effarr == NULL ) 
	{
		name1 = AllocateCharMtx( nseq1, B );
		name2 = AllocateCharMtx( nseq2, B );
		checkarr = AllocateIntVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		effarr = AllocateDoubleVec( njob );
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
#if 0
#else
#endif
	}


	if( tbweight == 0 )
	{
		for( i=0; i<njob; i++ ) 
		{
			effarr[i] = 1.0;
		}
	}
	else
	{
		for( i=0; i<njob; i++ ) 
		{
			effarr[i] = eff[i][i];
		}
	}
	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );

	for( i=0; i<nseq1; i++ ) checkarr[i] = i;
	checkarr[nseq1] = -1;
		

/*
	trap = fopen( "pre", "w" );
	if( !trap ) ErrorExit( 1 );
	WriteGapFill( trap, njob, name, nlen, aseq );
	fclose( trap );
*/


	if( tbweight ) 
	{
		for( l=0; l<njob-1; l++ ) 
		{
			for( j=0; j<2; j++ ) 
			{
#if DEBUG
				fprintf( stderr, "l=%d, check -> %d\n", l, intcmp( topol[l][j], checkarr ) );
#endif
				if( !intcmp( topol[l][j], checkarr ) ) goto jump;
			}
		}
		jump:
		if( l == njob-1 ) ErrorExit( "polyphylic" );
	}

	s1 = 0;
	s2 = nseq1;
//	fprintf( stdout, "nseq1 = %d\n", nseq1 );

	clus1 = conjuctionforgaln( 0, nseq1, aseq, mseq1, effarr1, effarr, name, name1, indication1 );
	clus2 = conjuctionforgaln( nseq1, njob, aseq, mseq2, effarr2, effarr, name, name2, indication2 );
	if( checkC )
		for( i=0; i<clus1; i++ ) fprintf( stderr, "eff for mseq1-%d%f\n", i, effarr1[i] );
/*
	fprintf( stderr, "before align all\n" );
	display( aseq, njob );
	fprintf( stderr, "\n" );
	fprintf( stderr, "before align 1 %s \n", indication1 );
	display( mseq1, clus1 );
	fprintf( stderr, "\n" );
	fprintf( stderr, "before align 2 %s \n", indication2 );
	display( mseq2, clus2 );
	fprintf( stderr, "\n" );
*/

	fprintf( trap_g, "group1 = %s\n", indication1 );
	fprintf( trap_g, "group2 = %s\n", indication2 );


	if( use_fft )
	{
		pscore = Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
	}
	else
	{
		if( alg == 'A' )
			pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
		else
		{
   	 	if( outgap && ( clus1 == 1 && clus2 != 1 || clus1 != 1 && clus2 == 1 ) )
				pscore = translate_and_Calign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
			else
				pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
		}
	}
	
/*
	fprintf( stderr, "after align 1 %s \n", indication1 );
	display( mseq1, clus1 );
	fprintf( stderr, "\n" );
	fprintf( stderr, "after align 2 %s \n", indication2 );
	display( mseq2, clus2 );
	fprintf( stderr, "\n" );
*/

	fprintf( stderr, "group-to-group %s /%s     %f\n", indication1, indication2, pscore );
	if( disp ) display( aseq, njob );
	fprintf( stderr, "\n" );

/*
	trap = fopen( "pre", "r+" );
	if( !trap ) ErrorExit( 1 );
	WriteGapFill( trap, njob, name, nlen, aseq );
	fclose( trap );
	fprintf( stdout, "nseq1 = %d\n", nseq1 );
*/
}

static void WriteOptions( FILE *fp )
{
	fprintf( fp, "tree-base method\n" );
	if( tbweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbweight == 3 ) fprintf( fp, "reversely weighted\n" );
	if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
	else if( scoremtx ==  1 ) fprintf( fp, "Dayhoff( machigai ga aru )\n" );
	else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	else if( scoremtx == -1 ) fprintf( fp, "DNA\n" );

	if( scoremtx )
		fprintf( fp, "Gap Penalty = %d, %d\n", penalty, offset );
	else
		fprintf( fp, "Gap Penalty = %d\n", penalty );
}
	 

int main( int argc, char *argv[] )
{
	static int  *nlen;	
	static char **name, **seq;
	static char **seq1, **seq2;
	static char **mseq1, **mseq2;
	static char **aseq;
	static char **bseq;
	static double **pscore;
	static double **node0, **node1;
	int i, j;
	int identity;
	static int ***topol;
	static double **len;
	FILE *prep;
	FILE *gp1, *gp2;
	char c;
	int nlenmax1, nlenmax2, nseq1, nseq2;
	int alloclen;

	argv = arguments( argc, argv );

	fprintf( stderr, "####### in galn\n" );

	initFiles();

	fprintf( stderr, "file1 = %s\n", argv[0] );
	fprintf( stderr, "file2 = %s\n", argv[1] );

	gp1 = fopen( argv[0], "r" ); if( !gp1 ) ErrorExit( "cannot open file1" );
	gp2 = fopen( argv[1], "r" ); if( !gp2 ) ErrorExit( "cannot open file2" );

#if 0
	PreRead( gp1, &nseq1, &nlenmax1 );
	PreRead( gp2, &nseq2, &nlenmax2 );
#else
    getnumlen( gp1 );
	nseq1 = njob; nlenmax1 = nlenmax;
    getnumlen( gp2 );
	nseq2 = njob; nlenmax2 = nlenmax;
#endif

	njob = nseq1 + nseq2;
	nlenmax = MAX( nlenmax1, nlenmax2 );

	rewind( gp1 );
	rewind( gp2 );

    if( nlenmax > N ) 
	{
		fprintf( stderr, "ERROR in main\n" );
	}


	name = AllocateCharMtx( njob, B );
	nlen = AllocateIntVec( njob );
	seq1 = AllocateCharMtx( nseq1, nlenmax*3 );
	seq2 = AllocateCharMtx( nseq2, nlenmax*3 );
	seq  = AllocateCharMtx( njob, 1 );
	aseq = AllocateCharMtx( njob, nlenmax*3 );
	bseq = AllocateCharMtx( njob, nlenmax*3 );
	mseq1 = AllocateCharMtx( njob, 1 );
	mseq2 = AllocateCharMtx( njob, 1 );
	alloclen = nlenmax * 3;

	if( tbweight )
	{
		topol = AllocateIntCub( njob, 2, njob );
		len = AllocateDoubleMtx( njob, 2 );
		pscore = AllocateDoubleMtx( njob, njob );
		node0 = AllocateDoubleMtx( njob, njob );
		node1 = AllocateDoubleMtx( njob, njob );
	}

#if 0
    njob=nseq2; FRead( gp2, name+nseq1, nlen+nseq1, seq2 );
	njob=nseq1; FRead( gp1, name, nlen, seq1 );
#else
    njob=nseq2; readDataforgaln( gp2, name+nseq1, nlen+nseq1, seq2 );
	njob=nseq1; readDataforgaln( gp1, name, nlen, seq1 );
#endif
	njob = nseq1 + nseq2;
	pamN = NOTSPECIFIED;


	commongappick( nseq1, seq1 );
	commongappick( nseq2, seq2 );

	for( i=0; i<nseq1; i++ ) seq[i] = seq1[i];
	for( i=nseq1; i<njob; i++ ) seq[i] = seq2[i-nseq1];
/*
	Write( stdout, njob, name, nlen, seq );
*/

    constants();

    WriteOptions( trap_g );

    c = seqcheck( seq );
    if( c )
    {
        fprintf( stderr, "Illeagal character %c\n", c );
        exit( 1 );
    }
    for( i=1; i<nseq1; i++ ) 
    {
        if( nlen[i] != nlen[0] ) 
            ErrorExit( "group1 is not aligned." );
    }
    for( i=nseq1+1;  i<njob; i++ ) 
    {
        if( nlen[i] != nlen[nseq1] ) 
            ErrorExit( "group2 is not aligned." );
    }
	if( tbweight ) 
	{
	    if( tbutree == 0 )
   	    {
			if( nlen[0] != nlen[nseq1] ) ErrorExit( "Input pre-aligned subgroups." );
   	        for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
   	        pscore[i][j] = (double)substitution_hosei( seq[i], seq[j] );
     	}
   	 	else
		{
#if 0
   	   	     prep = fopen( "hat2", "r" );
   	         if( prep == NULL ) ErrorExit( "Make hat2." );
   	         readhat2( prep, njob, name, pscore );
   	         fclose( prep );
#endif
   	    }
		/*
  	    upg2( njob, pscore, topol, len );
  	    upg2( njob, pscore, topol, len );
		*/
  	    spg( njob, pscore, topol, len );
   	    countnode( njob, topol, node0 );

        weight = tbweight;
        utree = 0; counteff( njob, topol, len, pscore ); utree = 1;
    }
    else
    {
#if 0
        for( i=0; i<njob; i++ ) pscore[i][i] = 1.0;
#endif
    }

	GroupAlign( nseq1, nseq2, name, nlen, seq, aseq, mseq1, mseq2, topol, len, pscore, alloclen );

#if 0
	writePre( njob, name, nlen, aseq, 1 );
#else
	writeDataforgaln( stdout, njob, name, nlen, aseq );
#endif

	if( tbweight )
	{
		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
		{
			pscore[i][j] = (double)substitution_hosei( aseq[i], aseq[j] );
		}
		spg( njob, pscore, topol, len );
		countnode( njob, topol, node1 );
	
		identity = 1;
		for( i=0; i<njob-1 && identity; i++ ) for( j=i+1; j<njob && identity; j++ ) 
			if ( node1[i][j] != node0[i][j] ) identity = 0;

		fprintf( trap_g, "done\n" );

		if( identity )
		{
			fprintf( stderr, "topology not changed\n" );
			exit( 0 );
		}
		else 
		{
			fprintf( stderr, "topology changed\n" );
			exit( 0 );
		}
	}
	else
	{
		fprintf( trap_g, "done\n" );
		return( 0 );
	}
}
