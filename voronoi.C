// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "dcd.h"
#include "ctype.h"
#include "alignSet.h"
#include "geometry.h"
#include "binary.h"
#include "voronoi.h"

void indexToColor( int i, double *color )
{
	color[0] = 0;
	color[1] = 0;
	color[2] = 0;

	int base_color_range = 1 + i % 6;

	if( base_color_range & (1 << 2 ) )
		color[0] += 1.0;
	if( base_color_range & (1 << 1 ) )
		color[1] += 1.0;
	if( base_color_range & (1 << 0 ) )
		color[2] += 1.0;

	double mix[3] = { 0.0, 0.0, 0.0 };

	double f_mix = ((i / 6) % 6) / 6.0;

	while( f_mix > 1 ) f_mix -= 1.0; 

	color[0] = color[0] * ( 1 - f_mix) + mix[0] * f_mix;	
	color[1] = color[1] * ( 1 - f_mix) + mix[1] * f_mix;	
	color[2] = color[2] * ( 1 - f_mix) + mix[2] * f_mix;	

/*	
	double mix2[3] = { 0.0, 1.0, 0.0 };
	double f_mix_2 = (i / 36.0);
	
	color[0] = color[0] * ( 1 - f_mix_2) + mix2[0] * f_mix_2;	
	color[1] = color[1] * ( 1 - f_mix_2) + mix2[1] * f_mix_2;	
	color[2] = color[2] * ( 1 - f_mix_2) + mix2[2] * f_mix_2;	
*/
}


void printDigits( int f, int spc, char *buf );
double writeFrame( double *atoms_in, int nat_in, double *colors, int seq, const char *unique, double Lx, double Ly, char *fileName, int *atCode, int doPS, int *borders, int *nborders, double *borderLens, double *borderThetas, int *solvationShell, double *area, struct specialColoring *coloring, int ncolor, double color_cut, int colorByOrderP, double orderPMin, double orderPMax, int chlAtCode, double *orderPOut );

void setColor( double *color, int seqn )
{
	color[0] = 0;
	color[1] = 0;
	color[2] = 0;

	switch( seqn % 3 )
	{
		case 0:
			color[0] = 1; // red
			break;
		case 2:
			color[1] = 1; // green
			break;
		case 1:
			color[2] = 1; // blue
			break;
		case 3:
			color[0] = 1; // yellow
			color[1] = 1; 
			break;
		case 4:
			color[0] = 1; // purple
			color[2] = 1; 
			break;
	}
}

int nChains( char *resname )
{
	if( !strcasecmp( resname, "PSM") ) return 2;
	if( !strcasecmp( resname, "LSM") ) return 2;
	if( !strcasecmp( resname, "NSM") ) return 2;
	if( !strcasecmp( resname, "PLAS") ) return 2;
	if( !strcasecmp( resname, "PLPC") ) return 2;
	if( !strncasecmp( resname, "CHL", 3 ) ) return 1;

	if( strlen( resname ) == 4 )
	{
		// some reasonable standard for it probably being a lipid, override however you like.

		int isPC = !strncasecmp( resname + 2, "PC",2 );
		int isPE = !strncasecmp( resname + 2, "PE",2 );
		int isPS = !strncasecmp( resname + 2, "PS",2 );

		int isDi = resname[0] == 'D';

		int isOleoyl = resname[1] == 'O';

		int isPalm = resname[1] == 'P';
		
		int isStear = resname[0] == 'S';

		if( (isPC || isPE || isPS) && (isDi || isOleoyl || isPalm || isStear ) )
			return 2;
	}
	return 0; // not a lipid.
}

int whichChain (char *resName, char *atomName )
{
	if( atomName[0] != 'C' ) return -1;
		
	int sl = strlen(atomName)-1;

	if( !strcasecmp( resName, "PSM" ) )
	{
		if( atomName[sl] == 'S' ) return 0;	
		if( atomName[sl] == 'F' ) return 1;	
	}
	else if( !strncasecmp( resName, "CHL", 3 ) )
	{
		return 0;
	}
	else 
	{	
		if( sl > 1 )
		{
			if( atomName[1] == '2' ) return 0;
			if( atomName[1] == '3' ) return 1;
		}
	}
	return -1; 
}



int main( int argc, char **argv )
{
	char buffer[4096];

	if( argc < 3)
	{
		printf("Syntax: voronoi psf dcd [options]\n");	
		printf("Options:\n");
		printf("TOP or BOTTOM : perform the analysis on the top or bottom leaflet [default TOP].\n");
		printf("PS=# : generate a postscript file (debug.ps) of the #th frame (starting from frame zero, will still loop through everything).\n"); 
		printf("BIN=fileName : dump a binary of the data to file fileName\n"); 
		printf("LOG=fileName : log the result to fileName [default output.log]\n"); 
		printf("WHOLE : perform the voronoi analysis for whole lipid instead of chain (default: chain)\n"); 
		printf("VERBOSE : verbose log.\n");
		printf("COLOR:RESNAME:RESID:CHAIN:RED:GREEN:BLUE : special coloring for this chain (zero or one) (RED:GREEN:BLUE from zero to one)\n");
		printf("BCOLOR:RESNAME:NBORDERS:RED:GREEN:BLUE : special coloring for this chain (zero or one) (RED:GREEN:BLUE from zero to one)\n");
		printf("CUTOFF=VALUE : value to use as the cutoff for the number of borders\n");
		printf("ORDERP:MIN:MAX : Color this by the order parameter, with white the max and black the min\n");
		printf("CLUSTERFILE=fileName : load clusters from the post-processed file.\n");
		printf("MOVIE : generate a PS file for each frame of the input.\n");
		printf("STEPS=INTEGER : do a movie frame for each (integer) number of dcd frames.\n");
		printf("INTERP=INTEGER : interpolate (integer) number of frames between each dcd frame used.\n");
		printf("--------\n");
		printf("Example: ./voronoi system.psf traj.dcd BOTTOM PS=0 COLOR:PSM:10:0:0.5:0:0.5 BIN=run.bin LOG=run.log\n");

		return 0;
	}

	int doTop = 1;
	char logFile[256];
	sprintf(logFile, "output.log" );
	int doBin = 0;
	char binFile[256];
	sprintf(binFile, "output.bin" );
	int psForFrame = -1;
	int verbose = 0;
	int doWhole = 0;
	int ncolor = 0;
	int ncolorSpace = 10;
	int doMovie = 0;
	int nmovie = 1;
	int ninterp = 5;
	int orderp = 0;
	double orderPMin = 0;
	double orderPMax = 0;
	double color_cut = 0;
	int chlAtCode = 0;
	double *orderPOut = NULL;
	char clusterFile[256];
	int colorByCluster = 0;

	specialColoring *scolors = (specialColoring *)malloc( sizeof(specialColoring) * ncolorSpace );

	for( int c = 3; c < argc; c++ )
	{
		if( !strncasecmp( argv[c], "TOP", 3 ) )
			doTop = 1;
		else if( !strncasecmp( argv[c], "BOTTOM", 6 ) )
			doTop = 0;
		else if( !strncasecmp( argv[c], "VERBOSE", 7 ) )
			verbose = 1;
		else if( !strncasecmp( argv[c], "WHOLE", 5 ) )
			doWhole = 1;
		else if( !strncasecmp( argv[c], "MOVIE", 5 ) )
			doMovie = 1;
		else if( !strncasecmp( argv[c], "INTERP=", 7 ) && strlen(argv[c]) > 7 )
			ninterp = atoi(argv[c]+7);
		else if( !strncasecmp( argv[c], "CUTOFF=", 7 ) && strlen(argv[c]) > 7 )
			color_cut = atoi(argv[c]+7);
		else if( !strncasecmp( argv[c], "STEPS=", 6 ) && strlen(argv[c]) > 6 )
			nmovie = atoi(argv[c]+6);
		else if( !strncasecmp( argv[c], "PS=", 3 ) && strlen(argv[c]) > 3 )
			psForFrame = atoi(argv[c]+3);
		else if( !strncasecmp( argv[c], "BIN=", 4) && strlen(argv[c]) > 4 )
		{
			doBin = 1;
			strcpy( binFile, argv[c]+4 );
		}
		else if( !strncasecmp( argv[c], "LOG=", 4) && strlen(argv[c]) > 4 )
			strcpy( logFile, argv[c]+4 );
		else if( !strncasecmp( argv[c], "CLUSTERFILE=", 12) && strlen(argv[c]) > 12 )
		{
			colorByCluster = 1;
			strcpy( clusterFile, argv[c]+12 );
		}
		else if( !strncasecmp( argv[c], "ORDERP", 6 )  )
		{
			orderp = 1;
			char * t = argv[c]+7;
			if( ! *t )
			{
				printf("Problem with order parameter option '%s'.\n", argv[c] );
				exit(1);
			}

			orderPMin = atof( t);

			while( *t != ':' && *t ) t += 1;

			t += 1;
			
			if( ! *t )
			{
				printf("Problem with order parameter option '%s'.\n", argv[c] );
				exit(1);
			}

			orderPMax = atof( t );

			printf("Coloring by the order parameter only between %lf and %lf.\n", orderPMin, orderPMax );
		}
		else if( !strncasecmp( argv[c], "COLOR", 5 ) || !strncasecmp( argv[c], "BCOLOR",6 ) )
		{
			if( ncolor == ncolorSpace )
			{
				ncolorSpace *= 2;
				scolors = (specialColoring *)realloc( scolors, sizeof(specialColoring) * ncolorSpace );
			}
			char *t = argv[c] + 5;
			int isBColor = 0;
			if( !strncasecmp( argv[c], "BCOLOR", 6) ) isBColor = 1;
			if( isBColor ) t += 1;		

			if( *t != ':' )
			{
				printf("Couldn't make any sense of first color argument '%s'.\n", 				
					argv[c] );
				return 0;
			}
			t+=1;
			int ix = 0;
			// the residue name;
			while( *t != ':' && *t )
			{
				scolors[ncolor].resname[ix] =*t;
				scolors[ncolor].resname[ix+1] = '\0';
				t += 1;
				ix++;
			}
			scolors[ncolor].nborders = 0;
			t += 1;

			if( !isBColor )
			{
				if( *t )
				{
					scolors[ncolor].res = atoi( t ); 
				}	
				else
				{
					printf("Couldn't make any sense of color argument '%s' (residue number).\n", 				
						argv[c] );
					return 0;
				}
				while( *t && *t != ':' ) t+= 1;
				t+=1;
				if( *t )
				{
					scolors[ncolor].chain = atoi( t ); 
				}	
				else
				{
					printf("Couldn't make any sense of color argument '%s' (chain number).\n", 				
						argv[c] );
					return 0;
				}
			}
			else
			{
				scolors[ncolor].nborders = atoi( t );

			}					
			while( *t != ':' && *t ) t += 1;	
			t+=1;

			if( *t )
			{
				scolors[ncolor].color[0] = atof( t ); 
			}	
			else
			{
				printf("Couldn't make any sense of color argument '%s' (red).\n", 				
					argv[c] );
				return 0;
			}
			
			while( *t && *t != ':' ) t+= 1;
			t+=1;

			if( *t )
			{
				scolors[ncolor].color[1] = atof( t ); 
			}	
			else
			{
				printf("Couldn't make any sense of color argument '%s' (green).\n", 				
					argv[c] );
				return 0;
			}
			
			while( *t && *t != ':' ) t+= 1;
			t+=1;

			if( *t )
			{
				scolors[ncolor].color[2] = atof( t ); 
			}	
			else
			{
				printf("Couldn't make any sense of color argument '%s' (green).\n", 				
					argv[c] );
				return 0;
			}

			scolors[ncolor].code = -1;

			ncolor++;
		}
		else
		{
			printf("Couldn't interpret command line option '%s'.\n", argv[c] );
		} 
	}
	
	FILE *psfFile = fopen(argv[1], "r" );
	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[1] );
		return 0;
	}
	if( !strcasecmp( argv[1] + strlen(argv[1])-3, "pdb" ) )
		loadPSFfromPDB( psfFile );	
	else
		loadPSF( psfFile );

	fclose(psfFile);	
	
	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );

	int init_done = 0;

	int nat = curNAtoms();
	
	double *cur_pos = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	memset( cur_pos, 0, sizeof(double) * curNAtoms() * 3 );
	
	FILE *theLog = fopen( logFile, "w");

	if( !theLog )
	{
		printf("Failed to open file '%s' for writing.\n", logFile );
		return 0;
	}

	struct indicator
	{
		char resname[256];
		char segid[256];
		int type;
		int nat;
		int res;
		int chain;
		int natSpace;
		int *atList;
		int isTop;	
	};

	int nind = 0;
	int nindSpace = 100;
	indicator *inds = (indicator *)malloc( sizeof(indicator) * nindSpace );

	double *pos = NULL;
	double *leaflet_pos = NULL;
	int *leaflet_link = NULL;
	int nleaflet = 0;
	int *leaflet_type = NULL;
	int *borders = NULL;
	double *borderLens = NULL;
	double *borderThetas = NULL;
	int *nborders = NULL;
	double *colors = NULL;
	double *lpos1=NULL, *lpos2=NULL, *lpos3=NULL, *lpos4=NULL;
	int *cl1 = NULL, *cl2 = NULL, *cl3 = NULL, *cl4 = NULL;
	double *luw=NULL;
	int is1=0,is2=0,is3=0,is4=0;
	FILE *binary = NULL;

	if( doBin )
	{
		binary = fopen(binFile, "wb");
		if( !binary )
		{
			printf("Failed to open file '%s' for binary writing.\n", binFile );
			return 0;
		}
	}

	int movieFrame = 0;

	// keep track of which code corresponds to cholesterol so that we don't color it with the order parameter
	chlAtCode = -1; 

	int *theCluster = NULL;
	FILE *cFile = NULL;
	if( colorByCluster )
	{
		cFile = fopen(clusterFile,"r");
		if( !cFile )
		{
			printf("Could not open cluster file '%s'.\n", clusterFile );
			return 0;
		}
	}

	for( int c = 2; c <= 2; c++ )
	{
		FILE *dcdFile = fopen(argv[c], "r");

		if( ! dcdFile )
		{
			printf("Couldn't open dcd file '%s'.\n", argv[c] );
			return 0;
		}
			
		readDCDHeader(dcdFile);
		setSymmetric();
		int nframes = curNFrames();

		for( int f = 0; f < nframes; f++ )
		{
			if( verbose ) 
			{
				fprintf(theLog, "Computing frame %d.\n", f );
				fflush(theLog);
			}

			double cur_com[2] = {0,0};
	
			double La, Lb, Lc;
			double alpha,beta,gamma;
			
			loadFrame( dcdFile, at );
			if( !DCDsuccess() )
			{
				nframes = f;
				break;
			}

			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );
	
			if( !init_done )
			{
				for( int a = 0; a < curNAtoms(); a++ )
				{
					int nc = nChains( at[a].resname );

					if( nc == 0 ) continue;

					if( doWhole ) nc = 1;

					int wc = whichChain( at[a].resname, at[a].atname );

					if( wc == -1 ) continue;
					
					if( doWhole ) wc = 0;

					int gotInd = -1;

					for( int i = 0; i < nind; i++ )
					{
						if( at[a].res == inds[i].res )
						{
							if( !strcmp( at[a].resname, inds[i].resname ) &&
							    !strcmp( at[a].segid, inds[i].segid ) &&
								wc == inds[i].chain )
							{
								gotInd = i;
								break;
							}
						} 
					}	

					if( gotInd == -1 )
					{
						if( nind == nindSpace )
						{
							nindSpace *= 2;
							inds = (indicator *)realloc( inds, sizeof(indicator) * nindSpace );
						} 

						gotInd = nind;

						inds[gotInd].natSpace = 10;
						inds[gotInd].atList = (int *)malloc( sizeof(int) * inds[gotInd].natSpace );
						strcpy( inds[gotInd].resname, at[a].resname );
						strcpy( inds[gotInd].segid, at[a].segid );
						inds[gotInd].res = at[a].res;
						inds[gotInd].nat = 0;
						inds[gotInd].chain = wc;	

						int gotType = -1;
						int maxType = -1;
						for( int i = 0; i < nind; i++ )
						{
							if( inds[i].type > maxType ) maxType = inds[i].type;

							if( !strcasecmp( inds[i].resname, inds[gotInd].resname ) )
								gotType = i; 
						}

						if( gotType != -1 )
							inds[gotInd].type = inds[gotType].type;
						else
						{
							if( !strncasecmp( inds[gotInd].resname, "CHL", 3 ) )
								chlAtCode = maxType+1;

							inds[gotInd].type = maxType + 1;
						}
						nind++;
					}	

					if( inds[gotInd].nat == inds[gotInd].natSpace )
					{
						inds[gotInd].natSpace *= 2;
						inds[gotInd].atList = (int *)realloc( inds[gotInd].atList, sizeof(int) * inds[gotInd].natSpace );
					}

					inds[gotInd].atList[inds[gotInd].nat] = a;
					inds[gotInd].nat += 1;
				}

				pos = (double *)malloc( sizeof(double) * nind * 3 );
				leaflet_pos = (double *)malloc( sizeof(double) * nind * 3 );
				leaflet_link = (int *)malloc( sizeof(int) * nind );
				leaflet_type = (int *)malloc( sizeof(int) * nind );
				borders = (int *)malloc( sizeof(int) * nind * nind );
				borderLens = (double *)malloc (sizeof(double) * nind * nind );
				borderThetas = (double *)malloc (sizeof(double) * nind * nind );
				nborders = (int *)malloc( sizeof(int) * nind );
				colors = (double*)malloc( sizeof(double) * 3 * nind );
				orderPOut = (double *)malloc( sizeof(double) * nind );
				theCluster = (int *)malloc( sizeof(int) * nind );

				cl1 = (int *)malloc( sizeof(int) * nind );
				cl2 = (int *)malloc( sizeof(int) * nind );
				cl3 = (int *)malloc( sizeof(int) * nind );
				cl4 = (int *)malloc( sizeof(int) * nind );

				lpos1 = (double *)malloc( sizeof(double) * 3 * nind );
				lpos2 = (double *)malloc( sizeof(double) * 3 * nind );
				lpos3 = (double *)malloc( sizeof(double) * 3 * nind );
				lpos4 = (double *)malloc( sizeof(double) * 3 * nind );
				luw = (double *)malloc( sizeof(double) * 3 * nind );
			}

			double avz = 0;

			for( int i = 0; i < nind; i++ )
			{
				pos[3*i+0] = 0;
				pos[3*i+1] = 0;
				pos[3*i+2] = 0;

				for( int ax = 0; ax < inds[i].nat; ax++ )
				{
					int a = inds[i].atList[ax];

					pos[3*i+0] += at[a].x;
					pos[3*i+1] += at[a].y;
					pos[3*i+2] += at[a].z;
				}

				pos[3*i+0] /= inds[i].nat;
				pos[3*i+1] /= inds[i].nat;
				pos[3*i+2] /= inds[i].nat;

//				avz += pos[3*i+2];
			}			

			if( init_done )
			{	
				// unwrap.
				for( int i = 0; i < nind; i++ )
				{
					while( pos[3*i+0] - luw[3*i+0] < -La/2 ) pos[3*i+0] += La;		
					while( pos[3*i+0] - luw[3*i+0] >  La/2 ) pos[3*i+0] -= La;		
					while( pos[3*i+1] - luw[3*i+1] < -Lb/2 ) pos[3*i+1] += Lb;		
					while( pos[3*i+1] - luw[3*i+1] >  Lb/2 ) pos[3*i+1] -= Lb;		
				}
			}



			double sub_com[3] = {0,0,0};

			for( int i = 0; i < nind; i++ )
			{
				sub_com[0] += pos[3*i+0];
				sub_com[1] += pos[3*i+1];
				sub_com[2] += pos[3*i+2];
			}

			sub_com[0] /= nind;
			sub_com[1] /= nind;
			sub_com[2] /= nind;
			
			for( int i = 0; i < nind; i++ )
			{
				pos[3*i+0]-=sub_com[0];
				pos[3*i+1]-=sub_com[1];
				pos[3*i+2]-=sub_com[2];
			}

			memcpy( luw, pos, sizeof(double) * 3 * nind );
//			avz /= nind;
			
			for( int i = 0; i < nind; i++ )
			{

				if( !init_done ) {
					leaflet_link[i] = -1;
					if (pos[3*i+2] > 0 ) 
					{
						if( doTop ) { leaflet_link[i] = nleaflet; leaflet_type[nleaflet] = inds[i].type; nleaflet++; } 
						inds[i].isTop = 1;
					}
					else
					{
						if( !doTop ) { leaflet_link[i] = nleaflet; leaflet_type[nleaflet] = inds[i].type; nleaflet++; } 
						inds[i].isTop = 0;
					}
				}

				while( pos[3*i+0] > La/2 )
					pos[3*i+0] -= La;
				while( pos[3*i+0] < -La/2 )
					pos[3*i+0] += La;

				while( pos[3*i+1] > Lb/2 )
					pos[3*i+1] -= Lb;
				while( pos[3*i+1] < -Lb/2 )
					pos[3*i+1] += Lb;
			}
	

			if( !init_done )
			{
				for( int i = 0; i < nind; i++ )
				{
					if( leaflet_link[i] == -1 ) continue;
					setColor( colors+leaflet_link[i]*3, inds[i].type ); 

					for( int s = 0; s < ncolor; s++ )
					{
						if( !strcasecmp( scolors[s].resname, inds[i].resname ) )
							scolors[s].code = inds[i].type; 

						if( !strcasecmp( scolors[s].resname, inds[i].resname ) &&
							scolors[s].chain == inds[i].chain &&
						        scolors[s].res == inds[i].res )
						{
							colors[leaflet_link[i]*3+0] = scolors[s].color[0];
							colors[leaflet_link[i]*3+1] = scolors[s].color[1];
							colors[leaflet_link[i]*3+2] = scolors[s].color[2];
						}
					}
				}

				if( binary ) 
				{
					double anticipatedSize = sizeof(binaryHeader) + nind * (double)sizeof(elementDescriptor) + nind * (double)sizeof(binaryElement) * (double)nframes;
					if( anticipatedSize > pow( 2.0, 31.0 ) )
						printf("WARNING: the anticipated binary output file size is greater than 2 GB.\n");
	
					binaryHeader header;

					header.formatCheck = FORMAT_CHECK;
					header.nind = nleaflet;
					header.nframes = nframes;
					
					fwrite( &header, 1, sizeof(binaryHeader), binary );
					fflush(binary);

					for( int i = 0; i < nind; i++ )
					{
						if( leaflet_link[i] == -1 ) continue;

						elementDescriptor elem;
						strcpy(elem.resname, inds[i].resname );
						strcpy(elem.segid, inds[i].segid );
						elem.res = inds[i].res;
						elem.chain = inds[i].chain;

						fwrite( &elem, 1, sizeof(elementDescriptor), binary );
					}
				}
			}
	
			init_done = 1;

			for( int i = 0; i < nind; i++  )
			{
				if( leaflet_link[i] == -1 ) continue;
				
				nborders[leaflet_link[i]] = 0;

				leaflet_pos[leaflet_link[i]*3+0] = pos[i*3+0];
				leaflet_pos[leaflet_link[i]*3+1] = pos[i*3+1];
				leaflet_pos[leaflet_link[i]*3+2] = pos[i*3+2];
			}

	

			char psFile[256];
			sprintf(psFile, "debug.ps");
			

			void putBezierPoint( double *p0, double *p1, double *p2, double *p3, double *out, double t );
	
			if( doMovie && (f % nmovie ) == 0 )
			{
				if( colorByCluster )
				{
					fread( theCluster, sizeof(int), nleaflet, cFile );

					memcpy( cl4, cl3, sizeof(int) * nleaflet );
					memcpy( cl3, cl2, sizeof(int) * nleaflet );
					memcpy( cl2, cl1, sizeof(int) * nleaflet );
					memcpy( cl1, theCluster, sizeof(int) * nleaflet );
				}

				memcpy( lpos4, lpos3,  sizeof(double) * nleaflet*3 );
				is4 = is3;
				memcpy( lpos3, lpos2,  sizeof(double) * nleaflet*3 );
				is3 = is2;
				memcpy( lpos2, lpos1,  sizeof(double) * nleaflet*3 );
				is2 = is1;
				memcpy( lpos1, leaflet_pos,  sizeof(double) * nleaflet *3);
				is1 = 1;

				if( is1 && is2 && is3 && is4 )
				{
					double *movie_write = (double *)malloc( sizeof(double) * 3 * nind );

					for( int f = 0; f < ninterp; f++ )
					{
						double t = ( f / (double)ninterp );

						for( int i = 0; i < nleaflet; i++ )
						{
							double p0[3] = { lpos4[3*i+0], lpos4[3*i+1], lpos4[3*i+2] };
							double p1[3] = { lpos3[3*i+0], lpos3[3*i+1], lpos3[3*i+2] };
							double p2[3] = { lpos2[3*i+0], lpos2[3*i+1], lpos2[3*i+2] };
							double p3[3] = { lpos1[3*i+0], lpos1[3*i+1], lpos1[3*i+2] };

							while( p2[0] - p3[0] < -La/2 ) p2[0] += La;
							while( p2[1] - p3[1] < -Lb/2 ) p2[1] += Lb;
							while( p2[0] - p3[0] >  La/2 ) p2[0] -= La;
							while( p2[1] - p3[1] >  Lb/2 ) p2[1] -= Lb;
							
							while( p1[0] - p2[0] < -La/2 ) p1[0] += La;
							while( p1[1] - p2[1] < -Lb/2 ) p1[1] += Lb;
							while( p1[0] - p2[0] >  La/2 ) p1[0] -= La;
							while( p1[1] - p2[1] >  Lb/2 ) p1[1] -= Lb;
							
							while( p0[0] - p1[0] < -La/2 ) p0[0] += La;
							while( p0[1] - p1[1] < -Lb/2 ) p0[1] += Lb;
							while( p0[0] - p1[0] >  La/2 ) p0[0] -= La;
							while( p0[1] - p1[1] >  Lb/2 ) p0[1] -= Lb;

							putBezierPoint( p0, p1, p2, p3, movie_write+3*i, t ); 
						}
				
						if( colorByCluster ) {
						for( int i = 0; i < nleaflet; i++ )
						{
							double color0[3], color1[3], color2[3], color3[3];

							if( cl4[i] >= 0 )
								indexToColor( cl4[i], color0 );
							else
							{
								if( inds[i].type == chlAtCode )
								{
									color0[0] = 0.0;
									color0[1] = 0.0;
									color0[2] = 0.0;
								} 
								else
								{
									color0[0] = 1.0;
									color0[1] = 1.0;
									color0[2] = 1.0;
								}
							}	
							
							if( cl3[i] >= 0 )
								indexToColor( cl3[i], color1 );
							else
							{
								if( inds[i].type == chlAtCode )
								{
									color1[0] = 0.0;
									color1[1] = 0.0;
									color1[2] = 0.0;
								} 
								else
								{
									color1[0] = 1.0;
									color1[1] = 1.0;
									color1[2] = 1.0;
								}
							}	
							
							if( cl2[i] >= 0 )
								indexToColor( cl2[i], color2 );
							else
							{
								if( inds[i].type == chlAtCode )
								{
									color2[0] = 0.0;
									color2[1] = 0.0;
									color2[2] = 0.0;
								} 
								else
								{
									color2[0] = 1.0;
									color2[1] = 1.0;
									color2[2] = 1.0;
								}
							}	
							
							if( cl1[i] >= 0 )
								indexToColor( cl2[i], color3 );
							else
							{
								if( inds[i].type == chlAtCode )
								{
									color3[0] = 0.0;
									color3[1] = 0.0;
									color3[2] = 0.0;
								} 
								else
								{
									color3[0] = 1.0;
									color3[1] = 1.0;
									color3[2] = 1.0;
								}
							}
							
							putBezierPoint( color0, color1, color2, color3, colors+3*i, t ); 
						}
						}
						char movieFilename[256];
						char numb[256];
						printDigits( movieFrame, 5, numb );
						sprintf( movieFilename, "%s.ps", numb );
						writeFrame( movie_write, nleaflet,  colors, f, "hi", La, Lb, movieFilename, leaflet_type, 1, borders, nborders, borderLens, borderThetas, NULL, NULL, scolors, ncolor, color_cut, orderp, orderPMin, orderPMax, chlAtCode, orderPOut ); 
						movieFrame++;
					}
					free(movie_write);
				}			
			}
						
			double *areas = (double *)malloc( sizeof(double) * nleaflet );
			memset( areas, 0, sizeof(double) * nleaflet );
			writeFrame( leaflet_pos, nleaflet,  colors, f, "hi", La, Lb, psFile, leaflet_type, ( f == psForFrame), borders, nborders, borderLens, borderThetas, NULL, areas, scolors, ncolor, color_cut, orderp, orderPMin, orderPMax, chlAtCode, orderPOut  ); 
	
			if( binary ) 
			{
				int formatCheck = FORMAT_CHECK;
				fwrite( &formatCheck, 1, sizeof(int), binary );
				binaryElement anElement;

				for( int i = 0; i < nind; i++ )
				{
					if( leaflet_link[i] == -1 ) continue;

					binaryElement elem;
					int off = leaflet_link[i];

					elem.type = inds[i].type;
					elem.op = orderPOut[off];
					elem.nborders = nborders[off];
					elem.area = areas[off];
					for( int b = 0; b < nborders[off]; b++ )
					{
						elem.borderLen[b] = borderLens[off*nleaflet+b];
						elem.borderTheta[b] = borderThetas[off*nleaflet+b];
						elem.borderID[b] = (short)borders[off*nleaflet+b];
					}	
					fwrite( &elem, 1, sizeof(binaryElement), binary );
				}
				fflush(binary);
			}
			free(areas);
			
			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
		}
		fclose(dcdFile); 
	}
	
	fclose(theLog);
	
}

void putBezierPoint( double *p0, double *p1, double *p2, double *p3, double *out, double t )
{
	double ctl0[3] = { p0[0], p0[1], p0[2] };	

	double ctl1[3] = { 
					(-5*p0[0] + 18 * p1[0] - 9 * p2[0] + 2 * p3[0])/6, 
					(-5*p0[1] + 18 * p1[1] - 9 * p2[1] + 2 * p3[1])/6, 
					(-5*p0[2] + 18 * p1[2] - 9 * p2[2] + 2 * p3[2])/6 }; 
	double ctl2[3] = { 
					(2*p0[0]  - 9* p1[0] +18 * p2[0] -5 *p3[0])/6, 
					(2*p0[1]  - 9* p1[1] +18 * p2[1] -5 *p3[1])/6, 
					(2*p0[2]  - 9* p1[2] +18 * p2[2] -5 *p3[2])/6 }; 

	double ctl3[3] = { p3[0], p3[1], p3[2] };	

	double use_t = 1.0 / 3.0 + t * (1.0/3.0);

	for( int d = 0; d < 3; d++ )
	{
		double a = ctl3[d] - 3 * ctl2[d] + 3 * ctl1[d] - ctl0[d];
		double b = 3 * ctl2[d] - 6 * ctl1[d] + 3 * ctl0[d];
		double c = 3 * ctl1[d] - 3 * ctl0[d];
		double pd = ctl0[d];

		out[d] =  a * use_t*use_t*use_t + b * use_t*use_t + c * use_t + pd;
	}
}

void printDigits( int f, int spc, char *buf )
{
	char tbuf[256];
	sprintf( tbuf, "%d", f );
	char *prt = buf;

	for( int x = strlen(tbuf); x < spc; x++ )
	{
		*prt = '0';
		prt += 1;
	}
	strcpy( prt, tbuf );
}
