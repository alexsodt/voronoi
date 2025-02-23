/*<html><pre>  -<a                             href="../libqhull/qh-qhull.htm"
  >-------------------------------</a><a name="TOP">-</a>

   qvoronoi.c
     compute Voronoi diagrams and furthest-point Voronoi
     diagrams using qhull

   see unix.c for full interface

   Copyright (c) 1993-2011, The Geometry Center
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "libqhull.h"
#include "mem.h"
#include "qset.h"
typedef struct
{
	int nvertices;
	int site1, site2;
	int *vertices;
} ridge;

ridge *ridgeStorage = NULL;
double *vertStorage = NULL;
int nridges = 0;
int nverts  = 0;

int nridge_space = 0;
int nvert_space = 0;
int write_to_globals = 1;

int *triStorage = NULL;
int ntri = 0;
int ntri_space = 0;

#if __MWERKS__ && __POWERPC__
#include <SIOUX.h>
#include <Files.h>
#include <console.h>
#include <Desk.h>

#elif __cplusplus
extern "C" {
  int isatty(int);
}



#elif _MSC_VER
#include <io.h>
#define isatty _isatty
int _isatty(int);

#else
int isatty(int);  /* returns 1 if stdin is a tty
                   if "Undefined symbol" this can be deleted along with call in main() */
#endif

/*-<a                             href="../libqhull/qh-qhull.htm#TOC"
  >-------------------------------</a><a name="prompt">-</a>

  qh_prompt
    long prompt for qhull

  notes:
    restricted version of libqhull.c

  see:
    concise prompt below
*/

/* duplicated in qvoron_f.htm and qvoronoi.htm
   QJ and Qt are deprecated, but allowed for backwards compatibility
*/
char hidden_options[]=" d n m v H U Qb QB Qc Qf Qg Qi Qm Qr QR Qv Qx TR E V Fa FA FC Fp FS Ft FV Pv Gt Q0 Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 ";

char qh_prompta[]= "\n\
qvoronoi- compute the Voronoi diagram\n\
    http://www.qhull.org  %s\n\
\n\
input (stdin):\n\
    first lines: dimension and number of points (or vice-versa).\n\
    other lines: point coordinates, best if one point per line\n\
    comments:    start with a non-numeric character\n\
\n\
options:\n\
    Qu   - compute furthest-site Voronoi diagram\n\
\n\
Qhull control options:\n\
    Qz   - add point-at-infinity to Voronoi diagram\n\
%s%s%s%s";  /* split up qh_prompt for Visual C++ */
char qh_promptb[]= "\
    Qs   - search all points for the initial simplex\n\
    QGn  - Voronoi vertices if visible from point n, -n if not\n\
    QVn  - Voronoi vertices for input point n, -n if not\n\
\n\
";
char qh_promptc[]= "\
Trace options:\n\
    T4   - trace at level n, 4=all, 5=mem/gauss, -1= events\n\
    Tc   - check frequently during execution\n\
    Ts   - statistics\n\
    Tv   - verify result: structure, convexity, and in-circle test\n\
    Tz   - send all output to stdout\n\
    TFn  - report summary when n or more facets created\n\
    TI file - input data from file, no spaces or single quotes\n\
    TO file - output results to file, may be enclosed in single quotes\n\
    TPn  - turn on tracing when point n added to hull\n\
     TMn - turn on tracing at merge n\n\
     TWn - trace merge facets when width > n\n\
    TVn  - stop qhull after adding point n, -n for before (see TCn)\n\
     TCn - stop qhull after building cone for point n (see TVn)\n\
\n\
Precision options:\n\
    Cn   - radius of centrum (roundoff added).  Merge facets if non-convex\n\
     An  - cosine of maximum angle.  Merge facets if cosine > n or non-convex\n\
           C-0 roundoff, A-0.99/C-0.01 pre-merge, A0.99/C0.01 post-merge\n\
    Rn   - randomly perturb computations by a factor of [1-n,1+n]\n\
    Wn   - min facet width for non-coincident point (before roundoff)\n\
\n\
Output formats (may be combined; if none, produces a summary to stdout):\n\
    s    - summary to stderr\n\
    p    - Voronoi vertices\n\
    o    - OFF format (dim, Voronoi vertices, and Voronoi regions)\n\
    i    - Delaunay regions (use 'Pp' to avoid warning)\n\
    f    - facet dump\n\
\n\
";
char qh_promptd[]= "\
More formats:\n\
    Fc   - count plus coincident points (by Voronoi vertex)\n\
    Fd   - use cdd format for input (homogeneous with offset first)\n\
    FD   - use cdd format for output (offset first)\n\
    FF   - facet dump without ridges\n\
    Fi   - separating hyperplanes for bounded Voronoi regions\n\
    FI   - ID for each Voronoi vertex\n\
    Fm   - merge count for each Voronoi vertex (511 max)\n\
    Fn   - count plus neighboring Voronoi vertices for each Voronoi vertex\n\
    FN   - count and Voronoi vertices for each Voronoi region\n\
    Fo   - separating hyperplanes for unbounded Voronoi regions\n\
    FO   - options and precision constants\n\
    FP   - nearest point and distance for each coincident point\n\
    FQ   - command used for qvoronoi\n\
    Fs   - summary: #int (8), dimension, #points, tot vertices, tot facets,\n\
                    for output: #Voronoi regions, #Voronoi vertices,\n\
                                #coincident points, #non-simplicial regions\n\
                    #real (2), max outer plane and min vertex\n\
    Fv   - Voronoi diagram as Voronoi vertices between adjacent input sites\n\
    Fx   - extreme points of Delaunay triangulation (on convex hull)\n\
\n\
";
char qh_prompte[]= "\
Geomview options (2-d only)\n\
    Ga   - all points as dots\n\
     Gp  -  coplanar points and vertices as radii\n\
     Gv  -  vertices as spheres\n\
    Gi   - inner planes only\n\
     Gn  -  no planes\n\
     Go  -  outer planes only\n\
    Gc   - centrums\n\
    Gh   - hyperplane intersections\n\
    Gr   - ridges\n\
    GDn  - drop dimension n in 3-d and 4-d output\n\
\n\
Print options:\n\
    PAn  - keep n largest Voronoi vertices by 'area'\n\
    Pdk:n - drop facet if normal[k] <= n (default 0.0)\n\
    PDk:n - drop facet if normal[k] >= n\n\
    Pg   - print good Voronoi vertices (needs 'QGn' or 'QVn')\n\
    PFn  - keep Voronoi vertices whose 'area' is at least n\n\
    PG   - print neighbors of good Voronoi vertices\n\
    PMn  - keep n Voronoi vertices with most merges\n\
    Po   - force output.  If error, output neighborhood of facet\n\
    Pp   - do not report precision problems\n\
\n\
    .    - list of all options\n\
    -    - one line descriptions of all options\n\
";
/* for opts, don't assign 'e' or 'E' to a flag (already used for exponent) */

/*-<a                             href="../libqhull/qh-qhull.htm#TOC"
  >-------------------------------</a><a name="prompt2">-</a>

  qh_prompt2
    synopsis for qhull
*/
char qh_prompt2[]= "\n\
qvoronoi- compute the Voronoi diagram.  Qhull %s\n\
    input (stdin): dimension, number of points, point coordinates\n\
    comments start with a non-numeric character\n\
\n\
options (qvoronoi.htm):\n\
    Qu   - compute furthest-site Voronoi diagram\n\
    Tv   - verify result: structure, convexity, and in-circle test\n\
    .    - concise list of all options\n\
    -    - one-line description of all options\n\
\n\
output options (subset):\n\
    s    - summary of results (default)\n\
    p    - Voronoi vertices\n\
    o    - OFF file format (dim, Voronoi vertices, and Voronoi regions)\n\
    FN   - count and Voronoi vertices for each Voronoi region\n\
    Fv   - Voronoi diagram as Voronoi vertices between adjacent input sites\n\
    Fi   - separating hyperplanes for bounded regions, 'Fo' for unbounded\n\
    G    - Geomview output (2-d only)\n\
    QVn  - Voronoi vertices for input point n, -n if not\n\
    TO file- output results to file, may be enclosed in single quotes\n\
\n\
examples:\n\
rbox c P0 D2 | qvoronoi s o         rbox c P0 D2 | qvoronoi Fi\n\
rbox c P0 D2 | qvoronoi Fo          rbox c P0 D2 | qvoronoi Fv\n\
rbox c P0 D2 | qvoronoi s Qu Fv     rbox c P0 D2 | qvoronoi Qu Fo\n\
rbox c G1 d D2 | qvoronoi s p       rbox c P0 D2 | qvoronoi s Fv QV0\n\
\n\
";
/* for opts, don't assign 'e' or 'E' to a flag (already used for exponent) */

/*-<a                             href="../libqhull/qh-qhull.htm#TOC"
  >-------------------------------</a><a name="prompt3">-</a>

  qh_prompt3
    concise prompt for qhull
*/
char qh_prompt3[]= "\n\
Qhull %s.\n\
Except for 'F.' and 'PG', upper-case options take an argument.\n\
\n\
 OFF_format     p_vertices     i_delaunay     summary        facet_dump\n\
\n\
 Fcoincident    Fd_cdd_in      FD_cdd_out     FF-dump-xridge Fi_bounded\n\
 Fxtremes       Fmerges        Fneighbors     FNeigh_region  FOptions\n\
 Fo_unbounded   FPoint_near    FQvoronoi      Fsummary       Fvoronoi\n\
 FIDs\n\
\n\
 Gvertices      Gpoints        Gall_points    Gno_planes     Ginner\n\
 Gcentrums      Ghyperplanes   Gridges        Gouter         GDrop_dim\n\
\n\
 PArea_keep     Pdrop d0:0D0   Pgood          PFacet_area_keep\n\
 PGood_neighbors PMerge_keep   Poutput_forced Pprecision_not\n\
\n\
 QG_vertex_good Qsearch_1st    Qupper_voronoi QV_point_good  Qzinfinite\n\
 T4_trace       Tcheck_often   Tstatistics    Tverify        Tz_stdout\n\
 TFacet_log     TInput_file    TPoint_trace   TMerge_trace   TOutput_file\n\
 TWide_trace    TVertex_stop   TCone_stop\n\
\n\
 Angle_max      Centrum_size   Random_dist    Wide_outside\n\
";

/*-<a                             href="../libqhull/qh-qhull.htm#TOC"
  >-------------------------------</a><a name="main">-</a>

  main( argc, argv )
    processes the command line, calls qhull() to do the work, and exits

  design:
    initializes data structures
    reads points
    finishes initialization
    computes convex hull and other structures
    checks the result
    writes the output
    frees memory
*/
int voronoi( double *pts, int npts, const char *unique ) {
  int curlong, totlong; /* used !qh_NOmem */
  int exitcode, numpoints, dim;
  coordT *points;
  boolT ismalloc;
  char filename[256];

  sprintf(filename, "qvo.%s.run", unique );
  FILE *dummy_out = fopen(filename,"w"); 
  sprintf(filename, "qvo.%s.err", unique );
  FILE *dummy_err = fopen(filename,"w"); 

  if( ridgeStorage != NULL )
  {
       int x;
	for( x = 0; x < nridges; x++ )
		free( ridgeStorage[x].vertices );
	free(ridgeStorage);
  }
  ridgeStorage = NULL;
  nridge_space = 0;
  nridges = 0;

  if( vertStorage != NULL )
	free(vertStorage);
  vertStorage = NULL;
  nvert_space = 100;
  nverts = 1;

  vertStorage = (double *)malloc( sizeof(double) * 3 * nvert_space );

  vertStorage[0] = -10.101;
  vertStorage[1] = -10.101;
  vertStorage[2] = -10.101;

  // I guess for now I want: s o Fv TO qvo2.result 

  int fake_argc = 6;

  char *fake_args[6];
  int x = 0;
  for( x = 0; x < 6; x++ ) fake_args[x] = (char *)malloc( sizeof(char) * 128 );

	sprintf(fake_args[0], "qvoronoi" );
	sprintf(fake_args[1], "s" );
	sprintf(fake_args[2], "o" );
	sprintf(fake_args[3], "Fv" );
	sprintf(fake_args[4], "TO" );
	sprintf(fake_args[5], "qvo2.%s.result", unique );

  qh_init_A(NULL,dummy_out, dummy_err, fake_argc, fake_args);  /* sets qh qhull_command */
  exitcode= setjmp(qh errexit); /* simple statement for CRAY J916 */
  if (!exitcode) {
    qh_option("voronoi  _bbound-last  _coplanar-keep", NULL, NULL);
    qh DELAUNAY= True;     /* 'v'   */
    qh VORONOI= True;
    qh SCALElast= True;    /* 'Qbb' */
    qh_checkflags(qh qhull_command, hidden_options);
    qh_initflags(qh qhull_command);
//    points= qh_readpoints(&numpoints, &dim, &ismalloc);
    dim = 4;
    numpoints = npts;

    points = (coordT *)malloc( sizeof(coordT) * dim * (npts+1) );

    for( x = 0; x < npts; x++ )
	{
		points[4*x+0] = pts[3*x+0];
		points[4*x+1] = pts[3*x+1];
		points[4*x+2] = pts[3*x+2];
		points[4*x+3] = pts[3*x+0]*pts[3*x+0]+pts[3*x+1]*pts[3*x+1]+pts[3*x+2]*pts[3*x+2];
	}

    if (dim >= 5) {
      qh_option("_merge-exact", NULL, NULL);
      qh MERGEexact= True; /* 'Qx' always */
    }
    qh_init_B(points, numpoints, dim, ismalloc);
    qh_qhull();
    qh_check_output();
    qh_produce_output();
    if (qh VERIFYoutput && !qh FORCEoutput && !qh STOPpoint && !qh STOPcone)
      qh_check_points();
    exitcode= qh_ERRnone;
  }
  qh NOerrexit= True;  /* no more setjmp */
#ifdef qh_NOmem
  qh_freeqhull( True);
#else
  qh_freeqhull( False);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    fprintf(stderr, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
       totlong, curlong);
#endif
  fclose(dummy_out);
  fclose(dummy_err);
  return exitcode;
} /* main */

int voronoi2D( double *pts, int npts, const char *unique ) {
  int curlong, totlong; /* used !qh_NOmem */
  int exitcode, numpoints, dim;
  coordT *points;
  boolT ismalloc;
  char filename[256];

  sprintf(filename, "qvo.%s.run", unique );
  FILE *dummy_out = fopen(filename,"w"); 
  sprintf(filename, "qvo.%s.err", unique );
  FILE *dummy_err = fopen(filename,"w"); 

  if( ridgeStorage != NULL )
  {
       int x;
	for( x = 0; x < nridges; x++ )
		free( ridgeStorage[x].vertices );
	free(ridgeStorage);
  }
  ridgeStorage = NULL;
  nridge_space = 0;
  nridges = 0;

  if( vertStorage != NULL )
	free(vertStorage);
  vertStorage = NULL;
  nvert_space = 100;
  nverts = 1;

  vertStorage = (double *)malloc( sizeof(double) * 3 * nvert_space );

  vertStorage[0] = -10.101;
  vertStorage[1] = -10.101;
  vertStorage[2] = -10.101;

  // I guess for now I want: s o Fv TO qvo2.result 

  int fake_argc = 6;

  char *fake_args[6];
  int x = 0;
  for( x = 0; x < 6; x++ ) fake_args[x] = (char *)malloc( sizeof(char) * 128 );

	sprintf(fake_args[0], "qvoronoi" );
	sprintf(fake_args[1], "s" );
	sprintf(fake_args[2], "o" );
	sprintf(fake_args[3], "Fv" );
	sprintf(fake_args[4], "TO" );
	sprintf(fake_args[5], "qvo2.%s.result", unique );

  qh_init_A(NULL,dummy_out, dummy_err, fake_argc, fake_args);  /* sets qh qhull_command */
  exitcode= setjmp(qh errexit); /* simple statement for CRAY J916 */
  if (!exitcode) {
    qh_option("voronoi  _bbound-last  _coplanar-keep", NULL, NULL);
    qh DELAUNAY= True;     /* 'v'   */
    qh VORONOI= True;
    qh SCALElast= True;    /* 'Qbb' */
    qh_checkflags(qh qhull_command, hidden_options);
    qh_initflags(qh qhull_command);
//    points= qh_readpoints(&numpoints, &dim, &ismalloc);
    dim = 3;
    numpoints = npts;

    points = (coordT *)malloc( sizeof(coordT) * dim * (npts+1) );

    for( x = 0; x < npts; x++ )
	{
		points[3*x+0] = pts[2*x+0];
		points[3*x+1] = pts[2*x+1];
		points[3*x+2] = pts[2*x+0]*pts[2*x+0]+pts[2*x+1]*pts[2*x+1];
	}

    if (dim >= 5) {
      qh_option("_merge-exact", NULL, NULL);
      qh MERGEexact= True; /* 'Qx' always */
    }
    qh_init_B(points, numpoints, dim, ismalloc);
    qh_qhull();
    qh_check_output();
    qh_produce_output();
    if (qh VERIFYoutput && !qh FORCEoutput && !qh STOPpoint && !qh STOPcone)
      qh_check_points();
    exitcode= qh_ERRnone;
  }
  qh NOerrexit= True;  /* no more setjmp */
#ifdef qh_NOmem
  qh_freeqhull( True);
#else
  qh_freeqhull( False);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    fprintf(stderr, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
       totlong, curlong);
#endif
  fclose(dummy_out);
  fclose(dummy_err);
  return exitcode;
} /* main */
