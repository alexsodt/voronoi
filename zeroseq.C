#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main( int argc, char **argv )
{
	if( argc < 3 )
	{
		printf("0\n");
		return 0;
	}

	int num = atoi(argv[1]);
	int spc = atoi(argv[2]);

	char buffer[256];

	sprintf( buffer, "%d", num );

	if( strlen(buffer) < spc )
	{
		for( int x = strlen(buffer); x < spc; x++ )
			printf("0");
	}
	printf("%s\n", buffer );

	return 0;
}
