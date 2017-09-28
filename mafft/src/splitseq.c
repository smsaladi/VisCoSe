#include <stdio.h>

main( int ac, char **av )
{
	int i, m, c, b, count, size;
	char name[100];
	FILE *ip, *op;


	if( ac == 3 )
	{
		size = (int)atoi( av[2] );
		if( size == 0 ) size = 5000;
	}
	else if( ac == 2 )
	{
		size = 5000;
	}
	else
	{
		fprintf( stderr, "Usage : %s file size\n", av[0] );
		exit( 1 );
	}

	ip = fopen( av[1], "r" );
	if( !ip )
	{
		fprintf( stderr, "Cannot open %s\n", av[1] );
		exit( 1 );
	}


	b = '\n';
	while( ( c = fgetc( ip ) ) != EOF )
	{
		if( b == '\n' && ( c == '>' || c == '=' ) )
			break;
	}
	ungetc( c, ip );

	m = 0;
	b = '\n';
	count = size;
	op = NULL;
	while( ( c = getc( ip ) ) != EOF )
	{
		if( b == '\n' && ( c == '>' || c == '=' ) )
		{
			if( count == size )
			{
				count = 0;
				if( op ) fclose( op );
				sprintf( name, "sp-%d", ++m );
				op = fopen( name, "w" );
				if( !op )
				{
					fprintf( stderr, "Cannot open %s\n", name );
					exit( 1 );
				}
				fprintf( stderr, "Writing %s\n", name );
			}
			count++;
		}
		fputc( c, op );
		b = c;
	}
	if( op ) fclose( op );
#if 1
	op = fopen( "sp-count", "w" );
	fprintf( op, "%d\n", m );
	fclose( op );
#endif
}
