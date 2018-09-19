#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "alignSet.h"
#ifdef MAC
#include "/System/Library/Frameworks/vecLib.framework/Headers/clapack.h"
#else
//#include "/v/apps/x86-64/pgi/linux86-64/6.0/include/acml.h"
//#include "/opt/acml3.6.0/pathscale64/include/acml.h"
extern "C" void dsytrf_( char *, int *, double *, int *, int *,  double *, int *, int *);
extern "C" void dsytri_( char *, int *, double *, int *, int *,  double *, int *);
extern "C" void dsyev_( char *, char *, int *, double *, int *, double *, double *, int *, int *);

#endif

struct DCDheader1
{
	char length;
	char hdr_string[7];
	int nfile;
	int npriv;
	int nsavc;
	int nstep; // icntrl 4
	int unk1[4];
	int fixed; // icntrl 9
	int unk2[1];
	int qcrys;
	int dim4;
	int qcg;
	int unk3[7]; // icntrl12-20
	int itemp;
	int unk4[1];
	int ntitl;
};

struct DCDheader
{
	DCDheader1 header1;
	char *title;
	int natom;	
};


static char **res_name;
static char **at_name;
static char **seg_name;
static int  *at_num;
static int  *res_num;
static int success = 0;

static int nwat_psf = 0;
static int psf_natoms = 0;
static int crd_type = 0;
static DCDheader curHeader;
static double FracINV[9];

double *mass = NULL;
int **blist = NULL;
int *nblist = NULL;

double loadTransform( double FINV[9] )
{
	memcpy( FracINV, FINV, sizeof(double) * 9 );
}
double saveTransform( double FINV[9] )
{
	memcpy( FINV, FracINV, sizeof(double) * 9 );
}

int DCDsuccess(void)
{
	return success;
}

int nhBonds( void )
{
	int nhb = 0;
	int is_h;
	for( int a = 0; a < psf_natoms; a++ )	
	{
		if( fabs(mass[a]) < 1.5 ) is_h = 1;

		for( int b = 0; b < nblist[a]; b++ )
		{
			int at_2 = blist[a][b];

			if( !is_h && fabs(mass[at_2]) < 1.5 )
				nhb += 1;
			if( is_h && fabs(mass[at_2]) < 1.5 && a < at_2 )
				nhb += 1;
		}
	}
	return nhb;
}

void getMasses( double *mass_out )
{
	memcpy( mass_out, mass, sizeof(double) * psf_natoms );
}

void loadPSF( FILE *theFile )
{
	char buffer[1024];

	getLine( theFile, buffer );	
	getLine( theFile, buffer );	
	getLine( theFile, buffer );	
	int nlines;
	sscanf( buffer, " %d ", &nlines );

	for( int x = 0; x < nlines; x++ )
		getLine( theFile, buffer );

	getLine( theFile, buffer );	
	getLine( theFile, buffer );	
	
	int local_nat = 0;
	sscanf( buffer, "%d ", &local_nat );
//	printf("Reading %d atoms from the PSF.\n", local_nat );

	res_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_name = (char **)malloc( sizeof( char *) * local_nat );	
	seg_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_num = (int *)malloc( sizeof(int) * local_nat );
	res_num = (int *)malloc( sizeof(int) * local_nat );

	mass = (double *)malloc( sizeof(double) * local_nat );
	blist = (int **)malloc( sizeof(int*) * local_nat );
	nblist = (int *)malloc( sizeof(int) * local_nat );

	for( int x = 0; x < local_nat; x++ )
	{
		nblist[x] = 0;
		blist[x] = (int *)malloc( sizeof(int) * 8 );
	}

	psf_natoms = local_nat;

	nwat_psf = 0;

	for( int x = 0; x < local_nat; x++ )
	{
		getLine( theFile, buffer );
		
		char l_segname[16];
		char l_resname[16];
		int l_resnum;
		int l_atnum;
		char l_atname[16]; 
		double charge, tmass;
		int at_id_code;
		sscanf( buffer, " %d %s %d %s %s %d %lf %lf ",
			&l_atnum, l_segname, &l_resnum, l_resname, l_atname, &at_id_code, &charge, &tmass );

		res_name[x] = (char *)malloc( sizeof(char) * ( strlen(l_resname) +1 ) );
		at_name[x] = (char *)malloc( sizeof(char) * ( strlen(l_atname) +1 ) );
		seg_name[x] = (char *)malloc( sizeof(char) * ( strlen(l_segname) +1 ) );
		mass[x] = tmass;
		

		strcpy( res_name[x], l_resname );	
		strcpy( at_name[x], l_atname );	
		strcpy( seg_name[x], l_segname );	

		if( !strncasecmp(res_name[x],"TIP",3) && at_name[x][0] == 'O' )
			nwat_psf += 1;

		at_num[x] = l_atnum;
		res_num[x] = l_resnum;
	}	
	
	getLine(theFile, buffer ); // blank line

	while( !feof(theFile) ) 
	{
		getLine(theFile, buffer ); // the read line

		int nunits;
		char unitName[256];
		int nr = sscanf(buffer, " %d !%s\n", &nunits, unitName );

		if( nr != 2 ) break;
		
		if( !strncasecmp( unitName, "NBOND", 5 ) )
		{
			// dihedrals.
			while( !feof(theFile) )
			{
				getLine( theFile, buffer );

				int units[8];
				int nr = sscanf(buffer, " %d %d %d %d %d %d %d %d",
					units+0,
					units+1,
					units+2,
					units+3,
					units+4,
					units+5,
					units+6,
					units+7 );

				for( int p = 0; p < nr; p++ )
					units[p] -= 1;

				for( int p = 0; p < nr/2; p++ )
				{
					int at_1 = units[2*p+0];
					int at_2 = units[2*p+1];

					// only do bond to previous list atom.
					
					blist[at_1][nblist[at_1]] = at_2;
					nblist[at_1] += 1;
					blist[at_2][nblist[at_2]] = at_1;
					nblist[at_2] += 1;
				}

				if( feof(theFile) || strlen(buffer) < 3 ) break;
			}
		}
		else
		{
			while( !feof(theFile) )
			{
				getLine( theFile, buffer );

				if( feof(theFile) || strlen(buffer) < 3 ) break;
			}
		}		
	}
}

static double celv = 0;

double CellVolume( void )
{
	return celv;	
}

double cstats[6];

int nwatPSF( void )
{
	return nwat_psf;
}

double PBCD( double *Lx, double *Ly, double *Lz, double *alpha, double *beta, double *gamma )
{
	*Lx = cstats[0];
	*Ly = cstats[1];
	*Lz = cstats[2];
	*alpha = cstats[3];
	*beta = cstats[4];
	*gamma = cstats[5];
}

void TransformFractional( double *dr )
{
	double tdr[3] = {0,0,0};

	tdr[0] = FracINV[0] * dr[0] + FracINV[1] * dr[1] + FracINV[2] * dr[2];
	tdr[1] = FracINV[3] * dr[0] + FracINV[4] * dr[1] + FracINV[5] * dr[2];
	tdr[2] = FracINV[6] * dr[0] + FracINV[7] * dr[1] + FracINV[8] * dr[2];
	
	dr[0] = tdr[0];
	dr[1] = tdr[1];
	dr[2] = tdr[2];
}

int curNFrames( void )
{
	return curHeader.header1.nstep / curHeader.header1.nsavc;
}

int curNAtoms( void )
{
	return psf_natoms;
}

void setFractional( void )
{
	crd_type = 0;
}

void setSymmetric( void )
{
	crd_type = 1;
}

void setAligned( void )
{
	crd_type = 2;
}


void readDCDHeader( FILE *theFile )
{
	fread( &(curHeader.header1), sizeof(DCDheader1), 1, theFile ); 
//	printf("QCRYS: %d DIM4: %d QCG: %d\n", curHeader.header1.qcrys, curHeader.header1.dim4, curHeader.header1.qcg );	
//	printf("title len: %d\n", curHeader.header1.ntitl );

	curHeader.header1.ntitl *= 80;

	curHeader.title = (char *)malloc( sizeof(char) * (curHeader.header1.ntitl+1) );

	fread( (curHeader.title), sizeof(char), curHeader.header1.ntitl, theFile );
	curHeader.title[curHeader.header1.ntitl] = '\0';
//	printf("title: %s\n", curHeader.title );	
	int junk;
	fread( &junk, sizeof(int), 1, theFile );
	fread( &junk, sizeof(int), 1, theFile );
	fread( &(curHeader.natom), sizeof(int), 1, theFile );
//	printf("natoms: %d\n", curHeader.natom );
	
	if( curHeader.header1.qcrys == 0 )
	{
//		fread( &junk, sizeof(int), 1, theFile ); //and eight more bytes for some reason. ??
	}
}

void loadFrame( FILE *theFile, struct atom_rec *at )
{
	double a[3];
	double b[3];
	double c[3];
	
	int junk;

	if( curHeader.header1.qcrys  )	fread( &junk, sizeof(int), 1, theFile );
//	printf("j1: %d\n", junk );	
	if( curHeader.header1.qcrys  )  fread( &junk, sizeof(int), 1, theFile );
//	printf("j2: %d\n", junk );	
	
	
	if( curHeader.header1.qcrys == 1 )
	{
		fread( &(a[0]), sizeof(double), 1, theFile );
		fread( &(a[1]), sizeof(double), 1, theFile );
		fread( &(b[1]), sizeof(double), 1, theFile );
		fread( &(a[2]), sizeof(double), 1, theFile );
		fread( &(b[2]), sizeof(double), 1, theFile );
		fread( &(c[2]), sizeof(double), 1, theFile );
		
	}
	else
	{
		a[0] = 1;
		a[1] = 0;
		b[1] = 1;
		b[2] = 0;
		c[2] = 1;
	}	
	
	double La,Lb,Lc;
	b[0] = a[1];
	c[0] = a[2];
	c[1] = b[2];
	La = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	Lb = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	Lc = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

	double gamma = (180.0/M_PI)*acos( a[1] );
//	double gamma = (180.0/M_PI)*acos( (a[0]*b[0] + a[1]*b[1]+a[2]*b[2]) / (La*Lb) );
	double alpha = (180.0/M_PI)*acos( (c[0]*b[0] + c[1]*b[1]+c[2]*b[2]) / (Lc*Lb) );
	double beta  = (180.0/M_PI)*acos( (a[0]*c[0] + a[1]*c[1]+a[2]*c[2]) / (La*Lc) );

	cstats[0] = La;
	cstats[1] = Lb;
	cstats[2] = Lc;
	cstats[3] = alpha;
	cstats[4] = beta;
	cstats[5] = gamma;
	
	celv = La*Lb*Lc*sin((M_PI/180.0)*gamma); // currently for hexagonal or orthorhombic.

	double SM[9];
	double SM_aligned[9];

	FracINV[0] = SM[0] = a[0];
	FracINV[1] = SM[1] = a[1];
	FracINV[2] = SM[2] = a[2];
	FracINV[3] = SM[3] = b[0];
	FracINV[4] = SM[4] = b[1];
	FracINV[5] = SM[5] = b[2];
	FracINV[6] = SM[6] = c[0];
	FracINV[7] = SM[7] = c[1];
	FracINV[8] = SM[8] = c[2];

	char uplo = 'U';
	int N = 3;
	int IPIV[N];
	int lwork = N*N;
	double work[lwork];
	int info;
	//dsytrf_( &uplo, &N, SM, &N, IPIV,  &info ); 
	//dsytri_( &uplo, &N, SM, &N, IPIV,  &info );
	dsytrf_( &uplo, &N, SM, &N, IPIV, work, &lwork, &info ); 
	dsytri_( &uplo, &N, SM, &N, IPIV, work, &info );

	SM[1] = SM[3];
	SM[2] = SM[6];
	SM[5] = SM[7];
/*
	for( int x = 0; x < 3; x++ )
	{
		for( int y = 0; y < 3; y++ )
		{
			printf("%lf ", SM[x*3+y] );
		}
		printf("\n");
	}

	printf("%lf %lf %lf :: %lf %lf %lf\n", La, Lb, Lc, alpha, beta, gamma );
*/
	// to align shape matrix, rotate around Z to diagonalize.

	double theta = atan2( -SM[1], SM[0] );

	SM_aligned[0] = cos(theta);
	SM_aligned[1] = sin(theta);
	SM_aligned[2] = 0;
	
	SM_aligned[3] = -sin(theta);
	SM_aligned[4] =  cos(theta);
	SM_aligned[5] = 0;
	
	SM_aligned[6] = 0;
	SM_aligned[7] = 0;
	SM_aligned[8] = 1;

#define TRAJ_TYPE float
	int natom = curHeader.natom;

	TRAJ_TYPE *loc_x = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_y = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_z = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	success = 1;
	fread( &junk, sizeof(int), 1, theFile );
	fread( &junk, sizeof(int), 1, theFile );
	size_t lread = fread( loc_x, sizeof(TRAJ_TYPE), natom, theFile );
        if( lread != natom )
		success = 0;
	fread( &junk, sizeof(int), 1, theFile );
	fread( &junk, sizeof(int), 1, theFile );
	lread = fread( loc_y, sizeof(TRAJ_TYPE), natom, theFile );
        if( lread != natom )
		success = 0;
	fread( &junk, sizeof(int), 1, theFile );
	fread( &junk, sizeof(int), 1, theFile );
	lread = fread( loc_z, sizeof(TRAJ_TYPE), natom, theFile );
        if( lread != natom )
		success = 0;

	for( int x = 0; x < natom; x++ )
	{
		double lx = loc_x[x];
		double ly = loc_y[x];
		double lz = loc_z[x];

		at[x].bead = at_num[x];
		at[x].res = res_num[x];
		at[x].resname = (char *)malloc( sizeof(char) * (strlen(res_name[x])+1) );
		at[x].atname = (char *)malloc( sizeof(char) * (strlen(at_name[x])+1) );
		at[x].segid = (char *)malloc( sizeof(char) * (strlen(seg_name[x])+1) );
		at[x].aux = 0;

		strcpy( at[x].resname, res_name[x] );
		strcpy( at[x].atname, at_name[x] );
		strcpy( at[x].segid, seg_name[x] );
		at[x].altloc = ' ';
		at[x].chain = ' ';
		if( crd_type == 0 )
		{ // fractional
			at[x].x = SM[0] * lx + SM[1] * ly + SM[2] * lz;
			at[x].y = SM[3] * lx + SM[4] * ly + SM[5] * lz;
			at[x].z = SM[6] * lx + SM[7] * ly + SM[8] * lz;
		}
		else if( crd_type == 1 )
		{ // symmetric
			at[x].x = lx;
			at[x].y = ly;
			at[x].z = lz;
		}	
		else if( crd_type == 2 )
		{
			at[x].x = (SM_aligned[0] * lx + SM_aligned[1] * ly + SM_aligned[2] * lz);
			at[x].y = (SM_aligned[3] * lx + SM_aligned[4] * ly + SM_aligned[5] * lz);
			at[x].z = (SM_aligned[6] * lx + SM_aligned[7] * ly + SM_aligned[8] * lz);
		}
	}

	free(loc_x);
	free(loc_y);
	free(loc_z);
	
}

void copyDCDHeader( FILE *theFile1, FILE *theFile2, int activateQCRYS )
{
	int nr = fread( &(curHeader.header1), sizeof(DCDheader1), 1, theFile1 ); 

	if( activateQCRYS) 
		curHeader.header1.qcrys = 1;

	fwrite( &(curHeader.header1), sizeof(DCDheader1), 1, theFile2 ); 

	curHeader.title = (char *)malloc( sizeof(char) * (curHeader.header1.ntitl*80+1) );
	nr = fread( (curHeader.title), sizeof(char), curHeader.header1.ntitl*80, theFile1 );
	fwrite( (curHeader.title), sizeof(char), curHeader.header1.ntitl*80, theFile2 );

	int junk;
	 nr = fread( &junk, sizeof(int), 1, theFile1 );
	printf("header junk1: %d\n", junk );
	fwrite( &junk, sizeof(int), 1, theFile2);
	nr = fread( &junk, sizeof(int), 1, theFile1);
	printf("header junk2: %d\n", junk );
	fwrite( &junk, sizeof(int), 1, theFile2 );

	nr= fread( &(curHeader.natom), sizeof(int), 1, theFile1 );
	fwrite( &(curHeader.natom), sizeof(int), 1, theFile2 );
	
	nr = fread( &junk, sizeof(int), 1, theFile1 ); // natom footer
	printf("header junk3: %d\n", junk );
	fwrite( &junk, sizeof(int), 1, theFile2 );
	
}

void loopATOM( 	int cur_at, TRAJ_TYPE *loc_x, TRAJ_TYPE *loc_y, TRAJ_TYPE *loc_z, int *bdone )
{	
	for( int b = 0; b < nblist[cur_at]; b++ )
	{
		if( bdone[blist[cur_at][b]] == 0 )
		{
			// match target to self.
			
			while( loc_x[blist[cur_at][b]] - loc_x[cur_at] < -0.5 ) loc_x[blist[cur_at][b]] += 1.0;
			while( loc_y[blist[cur_at][b]] - loc_y[cur_at] < -0.5 ) loc_y[blist[cur_at][b]] += 1.0;
			while( loc_z[blist[cur_at][b]] - loc_z[cur_at] < -0.5 ) loc_z[blist[cur_at][b]] += 1.0;
			
			while( loc_x[blist[cur_at][b]] - loc_x[cur_at] > 0.5 ) loc_x[blist[cur_at][b]] -= 1.0;
			while( loc_y[blist[cur_at][b]] - loc_y[cur_at] > 0.5 ) loc_y[blist[cur_at][b]] -= 1.0;
			while( loc_z[blist[cur_at][b]] - loc_z[cur_at] > 0.5 ) loc_z[blist[cur_at][b]] -= 1.0;
		}
	}
		
	bdone[cur_at] = 1; // all my partners are matched.
	
	for( int b = 0; b < nblist[cur_at]; b++ )
	{
		if( bdone[blist[cur_at][b]] == 0 )
			loopATOM( blist[cur_at][b], loc_x, loc_y, loc_z, bdone );
	}		
}

void copyVelocitiesNAMDtoCHARMM( FILE *theFile1, FILE *theFile2 )
{
	static int cntr = 0;

	int natom = curHeader.natom;

	TRAJ_TYPE *loc_x = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_y = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_z = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );


	double src_time_unit = (1e-12); // seconds
	double dst_time_unit = 4.888821E-14; // seconds.

	int junk1 = 48;
//	if( cntr > 0 ) junk1 = 110732;
//	cntr++;
//	fwrite( &junk1, sizeof(int), 1, theFile2 );
		
	fwrite( &junk1, sizeof(int), 1, theFile2 ); // the header
	double xtlabc[6] = {1,0,1,0,0,1};
	fwrite( xtlabc, sizeof(double), 6, theFile2 );
	fwrite( &junk1, sizeof(int), 1, theFile2 ); // the footer

	int junk;

	fread( &junk, sizeof(int), 1, theFile1 ); // the header
	fwrite( &junk, sizeof(int), 1, theFile2 );
	size_t lread = fread( loc_x, sizeof(TRAJ_TYPE), natom, theFile1 );
	fwrite( loc_x, sizeof(TRAJ_TYPE), natom, theFile2 );
	fread( &junk, sizeof(int), 1, theFile1 ); // the footer
	fwrite( &junk, sizeof(int), 1, theFile2 );

	fread( &junk, sizeof(int), 1, theFile1 ); // the header
	fwrite( &junk, sizeof(int), 1, theFile2 );
	lread = fread( loc_y, sizeof(TRAJ_TYPE), natom, theFile1 );
	fwrite( loc_y, sizeof(TRAJ_TYPE), natom, theFile2 );
	fread( &junk, sizeof(int), 1, theFile1 );
	fwrite( &junk, sizeof(int), 1, theFile2 ); // the footer


	fread( &junk, sizeof(int), 1, theFile1 ); // the header
	fwrite( &junk, sizeof(int), 1, theFile2 );
	lread = fread( loc_z, sizeof(TRAJ_TYPE), natom, theFile1 );
	fwrite( loc_z, sizeof(TRAJ_TYPE), natom, theFile2 );
	fread( &junk, sizeof(int), 1, theFile1 );
	fwrite( &junk, sizeof(int), 1, theFile2 ); // the footer

	free(loc_x);
	free(loc_y);
	free(loc_z);
	
//	junk = 0;
//	fwrite( &junk, sizeof(int), 1, theFile2 );
}	

void copyFrameNAMDtoCHARMM( FILE *theFile1, FILE *theFile2, int insertFakeQCRYS, int reimage, int watCen)
{
	double a[3];
	double b[3];
	double c[3];
	
	int junk;

	int nr ;
	for( int x = 0; x < sizeof(int); x++ )
	{
		char cjunk;
		fread( &cjunk, sizeof(char), 1, theFile1 );
//		printf("read value %d.\n", (int)(cjunk) );
		fwrite( &cjunk, sizeof(char), 1, theFile2 );
	}
	double La=1.0;
	double Lb=1.0;	
	double Lc=1.0;

	double gamma = 0;
	double beta  = 0;
	double alpha = 0;

	// read in the crystal dimension in NAMD format.
	if( insertFakeQCRYS ) // then don't read it..
	{
	}
	else
	{
		// the length of vec_A
		fread( &(a[0]), sizeof(double), 1, theFile1 );
		// the angle between vec_A and vec_B
		fread( &(a[1]), sizeof(double), 1, theFile1 );
		// the length of vec_B
		fread( &(b[1]), sizeof(double), 1, theFile1 );
		// the angle between ? .. let's call it vec_A and vec_C.
		fread( &(a[2]), sizeof(double), 1, theFile1 );
		// the angle between ? .. let's call it vec_B and vec_C.
		fread( &(b[2]), sizeof(double), 1, theFile1 );
		// the length of vec_C
		fread( &(c[2]), sizeof(double), 1, theFile1 );
	
		La = a[0];
		Lb = b[1];
		Lc = c[2];
	
		gamma = acos( a[1] );
		beta  = acos( a[2] );
		alpha = acos( b[2] );
	}

	// convert these to the symmetrized matrix.

	double xvec[3] = { La,              0,               0 };
	double yvec[3] = { Lb * cos(gamma), Lb * sin(gamma), 0 }; // rotate this around z into y.
	double zvec[3] = { 0,              0,                Lc };

	double hcur[9] =
	{	xvec[0], xvec[1], xvec[2],
		yvec[0], yvec[1], yvec[2],
		zvec[0], zvec[1], zvec[2] };
	

	double H[9];

	memset( H, 0, sizeof(double) * 9 );

	// h, htranspose.
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			H[i*3+j] = 0;
			for( int k = 0; k < 3; k++ )
				H[i*3+j] += hcur[i*3+k] * hcur[j*3+k];
		}
	}

	char jobz = 'V';
	char uplo = 'U';
	int order =3;
	double ev[3];	
	
	int info = 0;
	double owork = 0;

//	dsyev_(&jobz,&uplo, &order,H,&order,ev,&info);
	{
		int lwork = -1;	
		dsyev_(&jobz,&uplo,&order,H,&order,ev,&owork,&lwork,&info);
		double *work = (double *)malloc( sizeof(int)*(int)lround(owork));
	
		lwork = lround(owork);
		dsyev_(&jobz,&uplo,&order,H,&order,ev,work,&lwork,&info);
	
		free(work);
	}
	double A = sqrt(ev[0]);
	double B = sqrt(ev[1]);
	double C = sqrt(ev[2]);

	double xtlabc[6];

	xtlabc[0] = A * H[0*3+0] * H[0*3+0] + B*H[1*3+0] * H[1*3+0] + C*H[2*3+0]*H[2*3+0];
	xtlabc[2] = A * H[0*3+1] * H[0*3+1] + B*H[1*3+1] * H[1*3+1] + C*H[2*3+1]*H[2*3+1];
	xtlabc[5] = A * H[0*3+2] * H[0*3+2] + B*H[1*3+2] * H[1*3+2] + C*H[2*3+2]*H[2*3+2];
	xtlabc[1] = A * H[0*3+0] * H[0*3+1] + B*H[1*3+0] * H[1*3+1] + C*H[2*3+0]*H[2*3+1];
	xtlabc[3] = A * H[0*3+0] * H[0*3+2] + B*H[1*3+0] * H[1*3+2] + C*H[2*3+0]*H[2*3+2];
	xtlabc[4] = A * H[0*3+1] * H[2*3+2] + B*H[1*3+1] * H[1*3+2] + C*H[2*3+1]*H[2*3+2];

	a[0] = xtlabc[0];
	a[1] = xtlabc[1];
	b[1] = xtlabc[2];
	a[2] = xtlabc[3];
	b[2] = xtlabc[4];
	c[2] = xtlabc[5];	

	b[0] = a[1];
	c[0] = a[2];
	c[1] = b[2];

	double SM[9];
	double SM_aligned[9];

	FracINV[0] = SM[0] = a[0];
	FracINV[1] = SM[1] = a[1];
	FracINV[2] = SM[2] = a[2];
	FracINV[3] = SM[3] = b[0];
	FracINV[4] = SM[4] = b[1];
	FracINV[5] = SM[5] = b[2];
	FracINV[6] = SM[6] = c[0];
	FracINV[7] = SM[7] = c[1];
	FracINV[8] = SM[8] = c[2];

	uplo = 'U';
	int N = 3;
	int IPIV[N];
	int lwork = N*N;
	double work[lwork];
	info=0;
	dsytrf_( &uplo, &N, SM, &N, IPIV, work, &lwork, &info ); 
	dsytri_( &uplo, &N, SM, &N, IPIV, work, &info );
	//dsytrf_( &uplo, &N, SM, &N, IPIV,  &info ); 
	//dsytri_( &uplo, &N, SM, &N, IPIV,  &info );

	SM[1] = SM[3];
	SM[2] = SM[6];
	SM[5] = SM[7];
	
	printf("xtlabc: %lf %lf %lf %lf %lf %lf\n",
		xtlabc[0], xtlabc[1], xtlabc[2], xtlabc[3], xtlabc[4], xtlabc[5] );	
	fwrite( xtlabc, sizeof(double), 6, theFile2 );
	
#define TRAJ_TYPE float

	int natom = curHeader.natom;

	TRAJ_TYPE *loc_x = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_y = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_z = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );


	if( reimage )
	{
		double cen[3] = {0.5,0.5,0.5};
		int junk[3][2];
		
		fread( &(junk[0][0]), sizeof(int), 1, theFile1 );
		fread( &(junk[0][1]), sizeof(int), 1, theFile1 );
		size_t lread = fread( loc_x, sizeof(TRAJ_TYPE), natom, theFile1 );
		fread( &(junk[1][0]), sizeof(int), 1, theFile1 );
		fread( &(junk[1][1]), sizeof(int), 1, theFile1 );
		lread = fread( loc_y, sizeof(TRAJ_TYPE), natom, theFile1 );
		fread( &(junk[2][0]), sizeof(int), 1, theFile1 );
		fread( &(junk[2][1]), sizeof(int), 1, theFile1 );
		lread = fread( loc_z, sizeof(TRAJ_TYPE), natom, theFile1 );

		if( watCen )
		{
			cen[0] = 0;
			cen[1] = 0;
			cen[2] = 0;
			int nc = 0;
			for( int a = 0; a < natom; a++ )
			{
				double lx = loc_x[a];
				double ly = loc_y[a];
				double lz = loc_z[a];
		
				double tx  = 0, ty = 0, tz = 0;
	
				tx = SM[0] * lx + SM[1] * ly + SM[2] * lz;
				ty = SM[3] * lx + SM[4] * ly + SM[5] * lz;
				tz = SM[6] * lx + SM[7] * ly + SM[8] * lz;
			
				while( tx < -0.5) tx += 1.0; 
				while( ty >  0.5) ty -= 1.0; 
				while( tz < -0.5) tz += 1.0; 
				while( tx >  0.5) tx -= 1.0; 
				while( ty < -0.5) ty += 1.0; 
				while( tz >  0.5) tz -= 1.0; 

				if( !strcasecmp( at_name[a], "OH2" ) )	
				{	
					cen[0] += tx;
					cen[1] += ty;
					cen[2] += tz;
					nc += 1;
				}		
			}
			cen[0] /= nc;
			cen[1] /= nc;
			cen[2] /= nc;
		} 

		for( int a = 0; a < natom; a++ )
		{
			double lx = loc_x[a];
			double ly = loc_y[a];
			double lz = loc_z[a];
	
			double tx  = 0, ty = 0, tz = 0;

			tx = SM[0] * lx + SM[1] * ly + SM[2] * lz - cen[0];
			ty = SM[3] * lx + SM[4] * ly + SM[5] * lz - cen[1];
			tz = SM[6] * lx + SM[7] * ly + SM[8] * lz - cen[2];

			while( tx < -0.5) tx += 1.0; 
			while( ty >  0.5) ty -= 1.0; 
			while( tz < -0.5) tz += 1.0; 
			while( tx >  0.5) tx -= 1.0; 
			while( ty < -0.5) ty += 1.0; 
			while( tz >  0.5) tz -= 1.0; 

	
			double nx = 0, ny = 0, nz = 0;

			loc_x[a] = tx + cen[0];			
			loc_y[a] = ty + cen[1];			
			loc_z[a] = tz + cen[2];			
		}
			
		int *bdone = (int *)malloc( sizeof(int) * natom );
		memset( bdone, 0, sizeof(int) * natom );

		for( int a = 0; a < natom; a++ )
		{
			if( bdone[a] ) continue;
			
			loopATOM( a, loc_x, loc_y, loc_z, bdone );
		}

		free(bdone);

		for( int a = 0; a < natom; a++ )
		{
			double tx = loc_x[a];
			double ty = loc_y[a];
			double tz = loc_z[a];

			double nx = FracINV[0] * tx + FracINV[1] * ty + FracINV[2] * tz;
			double ny = FracINV[3] * tx + FracINV[4] * ty + FracINV[5] * tz;
			double nz = FracINV[6] * tx + FracINV[7] * ty + FracINV[8] * tz;			

			loc_x[a] = nx;
			loc_y[a] = ny;
			loc_z[a] = nz;
			
		}

		fwrite( &(junk[0][0]), sizeof(int), 1, theFile2 );
		fwrite( &(junk[0][1]), sizeof(int), 1, theFile2 );
		fwrite( loc_x, sizeof(TRAJ_TYPE), natom, theFile2 );
		fwrite( &(junk[1][0]), sizeof(int), 1, theFile2 );
		fwrite( &(junk[1][1]), sizeof(int), 1, theFile2 );
		fwrite( loc_y, sizeof(TRAJ_TYPE), natom, theFile2 );
		fwrite( &(junk[2][0]), sizeof(int), 1, theFile2 );
		fwrite( &(junk[2][1]), sizeof(int), 1, theFile2 );
		fwrite( loc_z, sizeof(TRAJ_TYPE), natom, theFile2 );
	}
	else
	{
		fread( &junk, sizeof(int), 1, theFile1 );
		fwrite( &junk, sizeof(int), 1, theFile2 );
		fread( &junk, sizeof(int), 1, theFile1 );
		fwrite( &junk, sizeof(int), 1, theFile2 );
		size_t lread = fread( loc_x, sizeof(TRAJ_TYPE), natom, theFile1 );
		fwrite( loc_x, sizeof(TRAJ_TYPE), natom, theFile2 );
		fread( &junk, sizeof(int), 1, theFile1 );
		fwrite( &junk, sizeof(int), 1, theFile2 );
		fread( &junk, sizeof(int), 1, theFile1 );
		fwrite( &junk, sizeof(int), 1, theFile2 );
		lread = fread( loc_y, sizeof(TRAJ_TYPE), natom, theFile1 );
		fwrite( loc_y, sizeof(TRAJ_TYPE), natom, theFile2 );
		fread( &junk, sizeof(int), 1, theFile1 );
		fwrite( &junk, sizeof(int), 1, theFile2 );
		fread( &junk, sizeof(int), 1, theFile1 );
		fwrite( &junk, sizeof(int), 1, theFile2 );
		lread = fread( loc_z, sizeof(TRAJ_TYPE), natom, theFile1 );
		fwrite( loc_z, sizeof(TRAJ_TYPE), natom, theFile2 );
	}

	free(loc_x);
	free(loc_y);
	free(loc_z);
	
	for( int x = 0; x < sizeof(int); x++ )
	{
		char cjunk;
		nr = fread( &cjunk, sizeof(char), 1, theFile1 );
//		printf("read value %d.\n", (int)(cjunk) );
		fwrite( &cjunk, sizeof(char), 1, theFile2 );
	}
}

void copyFrameNewCoords( FILE *theFile1, FILE *theFile2, double *coords )
{
	static double use_abc[6];
	static int abc_init = 0;

	double a[3];
	double b[3];
	double c[3];
	
	int junk;

	int nr ;
	for( int x = 0; x < sizeof(int); x++ )
	{
		char cjunk;
		fread( &cjunk, sizeof(char), 1, theFile1 );
		fwrite( &cjunk, sizeof(char), 1, theFile2 );
	}
	double La=1.0;
	double Lb=1.0;	
	double Lc=1.0;

	double gamma = 0;
	double beta  = 0;
	double alpha = 0;

	// the length of vec_A
	fread( &(a[0]), sizeof(double), 1, theFile1 );
	// the angle between vec_A and vec_B
	fread( &(a[1]), sizeof(double), 1, theFile1 );
	// the length of vec_B
	fread( &(b[1]), sizeof(double), 1, theFile1 );
	// the angle between ? .. let's call it vec_A and vec_C.
	fread( &(a[2]), sizeof(double), 1, theFile1 );
	// the angle between ? .. let's call it vec_B and vec_C.
	fread( &(b[2]), sizeof(double), 1, theFile1 );
	// the length of vec_C
	fread( &(c[2]), sizeof(double), 1, theFile1 );

	double xtlabc[6] = { a[0], a[1], b[1], a[2], b[2], c[2] };
	
	if( abc_init == 0 )
		memcpy( use_abc, xtlabc, sizeof(double) * 6 );
	abc_init = 1;
	memcpy( xtlabc, use_abc, sizeof(double) * 6 );
	fwrite( xtlabc, sizeof(double), 6, theFile2 );
	
#define TRAJ_TYPE float

	int natom = curHeader.natom;

	TRAJ_TYPE *loc_x = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_y = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_z = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	{	
	int junk[3][2];
	
	fread( &(junk[0][0]), sizeof(int), 1, theFile1 );
	fread( &(junk[0][1]), sizeof(int), 1, theFile1 );
	size_t lread = fread( loc_x, sizeof(TRAJ_TYPE), natom, theFile1 );
	fread( &(junk[1][0]), sizeof(int), 1, theFile1 );
	fread( &(junk[1][1]), sizeof(int), 1, theFile1 );
	lread = fread( loc_y, sizeof(TRAJ_TYPE), natom, theFile1 );
	fread( &(junk[2][0]), sizeof(int), 1, theFile1 );
	fread( &(junk[2][1]), sizeof(int), 1, theFile1 );
	lread = fread( loc_z, sizeof(TRAJ_TYPE), natom, theFile1 );

	for( int a = 0; a < natom; a++ )
	{
		loc_x[a] = (float)coords[3*a+0]; 
		loc_y[a] = (float)coords[3*a+1]; 
		loc_z[a] = (float)coords[3*a+2]; 
	}
		
	fwrite( &(junk[0][0]), sizeof(int), 1, theFile2 );
	fwrite( &(junk[0][1]), sizeof(int), 1, theFile2 );
	fwrite( loc_x, sizeof(TRAJ_TYPE), natom, theFile2 );
	fwrite( &(junk[1][0]), sizeof(int), 1, theFile2 );
	fwrite( &(junk[1][1]), sizeof(int), 1, theFile2 );
	fwrite( loc_y, sizeof(TRAJ_TYPE), natom, theFile2 );
	fwrite( &(junk[2][0]), sizeof(int), 1, theFile2 );
	fwrite( &(junk[2][1]), sizeof(int), 1, theFile2 );
	fwrite( loc_z, sizeof(TRAJ_TYPE), natom, theFile2 );
	}
	free(loc_x);
	free(loc_y);
	free(loc_z);
	
	for( int x = 0; x < sizeof(int); x++ )
	{
		char cjunk;
		nr = fread( &cjunk, sizeof(char), 1, theFile1 );
		fwrite( &cjunk, sizeof(char), 1, theFile2 );
	}
}

void loadPSFfromPDB( FILE *theFile )
{
	char buffer[1024];
	memset( buffer, 0, sizeof(char) * 1024 );

	int local_nat = 0;

	rewind(theFile);

	while( !feof(theFile) ) 
	{
		getLine( theFile, buffer );
		if( feof(theFile) ) break;

		if( !strncasecmp( buffer, "ATOM",4) ) 
			local_nat++;
	}

	rewind(theFile);

	res_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_name = (char **)malloc( sizeof( char *) * local_nat );	
	seg_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_num = (int *)malloc( sizeof(int) * local_nat );
	res_num = (int *)malloc( sizeof(int) * local_nat );

	psf_natoms = local_nat;

	int ncarbonSpace = 0;

	int x=0;
	while( !feof(theFile) ) 
	{
		getLine( theFile, buffer );
		
		if( feof(theFile) ) break;	
	
		if( !strncasecmp( buffer, "ATOM",4) )
		{ 
			struct atom_rec lat;

			readATOM( buffer, &lat );

		
	
			res_name[x] = (char *)malloc( sizeof(char) * ( strlen(lat.resname) +1 ) );
			at_name[x] = (char *)malloc( sizeof(char) * ( strlen(lat.atname) +1 ) );
			if( lat.segid )
				seg_name[x] = (char *)malloc( sizeof(char) * ( strlen(lat.segid) +1 ) );
			else
				seg_name[x] = (char *)malloc( sizeof(char) * ( strlen("NONE") + 1 ) );
	
			strcpy( res_name[x], lat.resname );	
			strcpy( at_name[x], lat.atname );	
			if( lat.segid ) 
				strcpy( seg_name[x], lat.segid );	
			else
				sprintf(seg_name[x], "NONE");
			at_num[x] = lat.bead;
			res_num[x] = lat.res;

			lat.zap();
			x++;
		}
	}	
}

