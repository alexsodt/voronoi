#ifndef __dcd_h__
#define __dcd_h__

#include "pdb.h"

void getDoubleBonds( int **dbonds, int *nbonds );
int getNPSFCarbons( void );
void putPSFCarbons( int *carbon_buffer );
void loadPSF( FILE *theFile );
void loadPSFfromPDB( FILE *theFile );
double PBCD( double *Lx, double *Ly, double *Lz, double *alpha, double *beta, double *gamma );
int curNAtoms( void );
int curNFrames( void );
void setFractional( void );
void setSymmetric( void );
void setAligned( void );
void readDCDHeader( FILE *theFile );
void loadFrame( FILE *theFile, struct atom_rec *at );
void TransformFractional( double *dr );
double CellVolume( void );
int DCDsuccess(void);
double loadTransform( double *);
double saveTransform( double *);
int getNPSFDihedrals( void );
void putPSFDihedrals( int *dihe_buffer );
void loadCRD( FILE *theFile, struct atom_rec *at);
void printCRD( FILE *theFile, struct atom_rec *at, int nat);
double getXTLABC( double SMOUT[9] );
void getBonds( int **bonds, int *nbonds );
#endif
