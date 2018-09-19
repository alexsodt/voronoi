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
#include "sys/time.h"
#include "voronoi.h"

double op_color_low[3] = { 0, 0, 0.01 }; // black, slightly blue otherwise there are color problems with JPGs
double op_color_high[3] = { 1, 1, 1 }; // white

#define WIDTH 612.0
#define HEIGHT 792.0

const char *ps_header = 
"%%!PS-Adobe-3.0\n"
"%%%%BoundingBox: 0 0 %lf %lf\n"
"%%%%Orientation: Portrait\n"
"%%%%Pages: 1\n"
"%%%%EndComments\n"
"%%%%Page: 1 1\n";

typedef struct
{
	int nvertices;
	int site1, site2;
	int *vertices;
} ridge;
extern ridge *ridgeStorage;
extern double *vertStorage;
extern int nridges;
extern int nverts ;

extern "C" int voronoi2D( double *pts, int npts, const char *unique );
void giftwrap( double *pts_in, double *pts_out, int npts, int *nptsout );
double loopArea( double *pts, int npts );

struct link
{
	int site1;
	int site2;

	double l;
	double edgeL;
};


#define EDGE_FUDGE	20.0
#define PS_WIDTH	WIDTH


double writeFrame( double *atoms_in, int nat_in, double *colors, int seq, const char *unique, double Lx, double Ly, char *fileName, int *atCode, int doPS, int *borders, int *nborders, double *borderLens, double *borderThetas, int *solvationShell, double *area, struct specialColoring *coloring, int nspec_col, double color_cut, int orderP, double orderPMin, double orderPMax, int codeIsChl, double *orderPOut )
{
	for( int x = 0; x < nat_in; x++ )
		orderPOut[x] = -1;


	if( area )
		memset( area, 0, sizeof(double) * nat_in );

	double offset_x = (WIDTH-PS_WIDTH)/2;
	double offset_y = (HEIGHT-PS_WIDTH)/2;
	double scale_factor = (PS_WIDTH)/Lx;

	if( nborders )
	{
		for( int x = 0; x < nat_in; x++ )
			nborders[x] = 0;
	}

	// atoms near the edge need to be periodically replicated.
	int nat = nat_in;

	// PBC wrap.
	for( int x = 0; x < nat_in; x++ )
	{
			while( atoms_in[3*x+0] < 0 ) atoms_in[3*x+0] += Lx; 
		while( atoms_in[3*x+0] > Lx ) atoms_in[3*x+0] -= Lx; 
		while( atoms_in[3*x+1] < 0 ) atoms_in[3*x+1] += Ly; 
		while( atoms_in[3*x+1] > Ly ) atoms_in[3*x+1] -= Ly; 
	}


	for( int x = 0; x < nat_in; x++ )
	{
		int pbc_x = 0;
		int pbc_y = 0;
		if( atoms_in[3*x+0] < EDGE_FUDGE )
			pbc_x = 1;

		if( atoms_in[3*x+0] > Lx- EDGE_FUDGE )
			pbc_x = 1;

		if( atoms_in[3*x+1] < EDGE_FUDGE )
			pbc_y = 1;

		if( atoms_in[3*x+1] > Ly- EDGE_FUDGE )
			pbc_y = 1;

		if( pbc_x ) nat += 1;
		if( pbc_y ) nat += 1;
		if( pbc_y && pbc_x ) nat += 1;
	}
	
	double *atoms = (double *)malloc( sizeof(double) * 2 * nat );
	int *at_link = (int *)malloc( sizeof(int) * nat );

	nat = 0;

	for( int x = 0; x < nat_in; x++ )
	{
		at_link[x] = x;
		atoms[2*nat+0] = atoms_in[3*x+0];
		atoms[2*nat+1] = atoms_in[3*x+1];

		nat++;
	}
	

	for( int x = 0; x < nat_in; x++ )
	{
		double del_x = 0;
		double del_y = 0;
		int pbc_x = 0;
		int pbc_y = 0;

		if( atoms_in[3*x+0] < EDGE_FUDGE )
		{
			del_x = Lx;
			pbc_x = 1;
		}

		if( atoms_in[3*x+0] > Lx- EDGE_FUDGE )
		{
			del_x = -Lx;
			pbc_x = 1;
		}

		if( atoms_in[3*x+1] < EDGE_FUDGE )
		{
			del_y = Ly;
			pbc_y = 1;
		}

		if( atoms_in[3*x+1] > Ly- EDGE_FUDGE )
		{
			del_y = -Ly;
			pbc_y = 1;
		}

		if( pbc_x ) 
		{
			atoms[2*nat+0] = atoms_in[3*x+0] + del_x;
			atoms[2*nat+1] = atoms_in[3*x+1];
			at_link[nat] = x;
			nat += 1;
		}

		if( pbc_y )
		{
			atoms[2*nat+0] = atoms_in[3*x+0];
			atoms[2*nat+1] = atoms_in[3*x+1] + del_y;
			at_link[nat] = x;
			nat += 1;
		}

		if( pbc_y && pbc_x )
		{
			atoms[2*nat+0] = atoms_in[3*x+0] + del_x;
			atoms[2*nat+1] = atoms_in[3*x+1] + del_y;
			at_link[nat] = x;
			nat += 1;
		}
	}

	voronoi2D( atoms, nat, unique );

	// for each atom, do it's volume and shade it.

	// set up the post script file.

	ridge *ridges = ridgeStorage;
	double *allVertices = vertStorage;
			
	FILE *psFile = NULL;

	if( doPS )
	{	
		psFile = fopen(fileName, "w");
	
		fprintf(psFile, ps_header, WIDTH, HEIGHT );
	
		fprintf(psFile, "%lf %lf translate\n", offset_x, offset_y );
	
		fprintf(psFile, "newpath\n");
		fprintf(psFile, "%lf %lf moveto\n", 0., 0. );
		fprintf(psFile, "%lf %lf lineto\n", scale_factor * Lx, 0. );
		fprintf(psFile, "%lf %lf lineto\n", scale_factor * Lx, scale_factor * Ly );
		fprintf(psFile, "%lf %lf lineto\n", 0., scale_factor * Ly );
		fprintf(psFile, "closepath\n");
		fprintf(psFile, "clip\n");
	}

	double *orderpc = (double *)malloc( sizeof(double ) * nat );
	double *orderps = (double *)malloc( sizeof(double ) * nat );
	double *orderpn = (double *)malloc( sizeof(double ) * nat );
	memset( orderpc, 0, sizeof(double) * nat );
	memset( orderps, 0, sizeof(double) * nat );
	memset( orderpn, 0, sizeof(double) * nat );

	for( int a = 0; a < nat; a++ )
	{
		int npts_convex = 0;

		int bad_ad = 0;
				
		if( borders && nborders )
		{
			for( int x = 0; x < nridges; x++ )
			{
				int nvertsl = ridges[x].nvertices;
				int s1 = ridges[x].site1;
				int s2 = ridges[x].site2;
				int *verts = ridges[x].vertices;
			
				if( s1 == a || s2 == a )
				{
					int l = at_link[s1];
					int l2 = at_link[s2];
	
					if( at_link[s1] != s1 && at_link[s2] != s2 ) continue;
	
	                                double *v1 = vertStorage+3*verts[0];    
	                                double *v2 = vertStorage+3*verts[1];    
	                                double dv[2] = { v1[0] - v2[0], v1[1] - v2[1] };
					double dr = sqrt(dv[0]*dv[0]+dv[1]*dv[1]);

					double bv[2] = { atoms[2*s1+0] - atoms[2*s2+0],
							  atoms[2*s1+1] - atoms[2*s2+1] };

					double theta = atan2( bv[1], bv[0] );	

					int gotit = 0;
					for( int p = 0; p < nborders[l]; p++ )
					{
						if( borders[l*nat_in+p] == l2 )
							gotit = 1; 	
					}
					if( !gotit )
					{
						if(borderLens ) borderLens[l*nat_in+nborders[l]] = dr;
						if(borderThetas ) borderThetas[l*nat_in+nborders[l]] = theta;
						borders[l*nat_in+nborders[l]] = l2;
						nborders[l] += 1;
						
						if( dr > color_cut )
						{
							orderpc[l] += cos( 6 * theta );
							orderps[l] += sin( 6 * theta );
							orderpn[l] += 1;
						}
					}
					gotit = 0;
					for( int p = 0; p < nborders[l2]; p++ )
					{
						if( borders[l2*nat_in+p] == l )
							gotit = 1; 	
					}
					if( !gotit )
					{
						if(borderLens ) borderLens[l2*nat_in+nborders[l2]] = dr;
						if(borderThetas ) borderThetas[l2*nat_in+nborders[l2]] = theta;
						borders[l2*nat_in+nborders[l2]] = l;
						nborders[l2] += 1;

						if( dr > color_cut )
						{
							orderpc[l2] += cos( 6 * theta );
							orderps[l2] += sin( 6 * theta );
							orderpn[l2] += 1;
						}
					}
				}
			}
		}
	}

	for( int a = 0; a < nat; a++ )
	{
		int npts_convex = 0;

		int bad_ad = 0;
				
	/*	if( borders && nborders )
		{
			for( int x = 0; x < nridges; x++ )
			{
				int nvertsl = ridges[x].nvertices;
				int s1 = ridges[x].site1;
				int s2 = ridges[x].site2;
				int *verts = ridges[x].vertices;
			
				if( s1 == a || s2 == a )
				{
					int l = at_link[s1];
					int l2 = at_link[s2];
	
					if( at_link[s1] != s1 && at_link[s2] != s2 ) continue;
	
	                                double *v1 = vertStorage+3*verts[0];    
	                                double *v2 = vertStorage+3*verts[1];    
	                                double dv[2] = { v1[0] - v2[0], v1[1] - v2[1] };
					double dr = sqrt(dv[0]*dv[0]+dv[1]*dv[1]);
	
					int gotit = 0;
					for( int p = 0; p < nborders[l]; p++ )
					{
						if( borders[l*nat_in+p] == l2 )
							gotit = 1; 	
					}
					if( !gotit )
					{
						if(borderLens ) borderLens[l*nat_in+nborders[l]] = dr;
						borders[l*nat_in+nborders[l]] = l2;
						nborders[l] += 1;
					}
					gotit = 0;
					for( int p = 0; p < nborders[l2]; p++ )
					{
						if( borders[l2*nat_in+p] == l )
							gotit = 1; 	
					}
					if( !gotit )
					{
						if(borderLens ) borderLens[l2*nat_in+nborders[l2]] = dr;
						borders[l2*nat_in+nborders[l2]] = l;
						nborders[l2] += 1;
					}
				}
			}
		} */
		if( area && at_link[a] == a  )
		{		
			npts_convex = 0;

			int bad_ad = 0;	
			
			for( int x = 0; x < nridges; x++ )
			{
				int nvertsl = ridges[x].nvertices;
				int s1 = ridges[x].site1;
				int s2 = ridges[x].site2;
				int *verts = ridges[x].vertices;
		
				if( s1 == a || s2 == a )
				{
					npts_convex += nvertsl;
					for( int v = 0; v < nvertsl; v++ )
					{
						if( verts[v] == 0 )
							bad_ad = 1;
					}	 
				}
			}

			if( !bad_ad )
			{
				double *pts = (double *)malloc( sizeof(double ) * 3 * npts_convex );
				
				int ctr = 0;
				for( int x = 0; x < nridges; x++ )
				{
					int nvertsl = ridges[x].nvertices;
					int s1 = ridges[x].site1;
					int s2 = ridges[x].site2;
					int *verts = ridges[x].vertices;
			
					if( s1 == a || s2 == a )
					{
						int l = at_link[s1];
						int l2 = at_link[s2];
	
						for( int v = 0; v < nvertsl; v++ )
						{
							pts[3*ctr+0] = allVertices[verts[v]*3+0]; 
							pts[3*ctr+1] = allVertices[verts[v]*3+1];
							pts[3*ctr+2] = 0;
		
							ctr++; 
						}	 
					}
				}
		
				double *pts_out = (double *)malloc( sizeof(double) * 3 * npts_convex );
		
				int ngf = 0;
		
				giftwrap(pts, pts_out, npts_convex, &ngf );
		
				double the_area = loopArea( pts_out, ngf );

				area[at_link[a]] = the_area;

				free(pts_out);
				free(pts);
			}
		}
		
		{
			int oa = at_link[a];
			double val =  orderpc[oa] * orderpc[oa] + 
				      orderps[oa] * orderps[oa];
			val = fabs(val);
			val = sqrt(val);

			val /= orderpn[oa];

			if( ! ( val < 1 ) && !( val >0 ) )
			{
				printf("nan");
			}

			orderPOut[at_link[a]] = val;
		}

		if( doPS )
		{		
			npts_convex = 0;

			int bad_ad = 0;	
			
			for( int x = 0; x < nridges; x++ )
			{
				int nvertsl = ridges[x].nvertices;
				int s1 = ridges[x].site1;
				int s2 = ridges[x].site2;
				int *verts = ridges[x].vertices;
		
				if( s1 == a || s2 == a )
				{
					npts_convex += nvertsl;
					for( int v = 0; v < nvertsl; v++ )
					{
						if( verts[v] == 0 )
							bad_ad = 1;
					}	 
				}
			}

			if( !bad_ad )
			{
				double *pts = (double *)malloc( sizeof(double ) * 3 * npts_convex );
				
				int ctr = 0;
				for( int x = 0; x < nridges; x++ )
				{
					int nvertsl = ridges[x].nvertices;
					int s1 = ridges[x].site1;
					int s2 = ridges[x].site2;
					int *verts = ridges[x].vertices;
			
					if( s1 == a || s2 == a )
					{
						int l = at_link[s1];
						int l2 = at_link[s2];
	
						for( int v = 0; v < nvertsl; v++ )
						{
							pts[3*ctr+0] = allVertices[verts[v]*3+0]; 
							pts[3*ctr+1] = allVertices[verts[v]*3+1];
							pts[3*ctr+2] = 0;
		
							ctr++; 
						}	 
					}
				}
		
				double *pts_out = (double *)malloc( sizeof(double) * 3 * npts_convex );
		
				int ngf = 0;
		
				giftwrap(pts, pts_out, npts_convex, &ngf );
		
				double the_area = loopArea( pts_out, ngf );

				fprintf(psFile, "newpath\n");

				int special_coloring = 0;
				int nbc = 0;
				int oa = at_link[a];

				if( borders && nborders )
				{
					for( int bx = 0; bx < nborders[oa]; bx++ )
					{
						if( borderLens[oa*nat_in+bx] > color_cut )
							nbc++;
					}
				}	

				double use_color[3] = { colors[at_link[a]*3+0], colors[at_link[a]*3+1], colors[at_link[a]*3+2] };
				
				for( int s = 0; s < nspec_col; s++ )
				{
					if( coloring[s].code == atCode[oa] && coloring[s].nborders == nbc )
					{
						use_color[0] = coloring[s].color[0];
						use_color[1] = coloring[s].color[1];
						use_color[2] = coloring[s].color[2];
					}
				}

//				if( atCode[at_link[a]] != codeIsChl )
				{
					int oa = at_link[a];
					double val =  orderpc[oa] * orderpc[oa] + 
						      orderps[oa] * orderps[oa];
					val = fabs(val);
					val = sqrt(val);

					val /= orderpn[oa];

					if( ! ( val < 1 ) && !( val >0 ) )
					{
						printf("nan");
					}

					orderPOut[at_link[a]] = val;

					double t = (val - orderPMin ) / (orderPMax - orderPMin);

					if( t > 0.9 ) t = 1.0;
					else t = 0.0;

					if( t < 0 ) t = 0;
					if( t > 1 ) t = 1;

					if( orderP && atCode[at_link[a]] != codeIsChl )
					{
						use_color[0] = (1-t) * op_color_low[0] + t * op_color_high[0]; 
						use_color[1] = (1-t) * op_color_low[1] + t * op_color_high[1]; 
						use_color[2] = (1-t) * op_color_low[2] + t * op_color_high[2]; 
					}
				}

				if( solvationShell )
				{
					if( solvationShell[at_link[a]] % 2 == 0 )
						fprintf(psFile, "1. 0. 0. setrgbcolor\n");
					else
						fprintf(psFile, "0. 1. 0. setrgbcolor\n" );
					
				}
				else
				//	fprintf(psFile, "%lf %lf %lf setrgbcolor\n", colors[at_link[a]*3+0], colors[at_link[a]*3+1], colors[at_link[a]*3+2] );
					fprintf(psFile, "%lf %lf %lf setrgbcolor\n", use_color[0], use_color[1], use_color[2] );
		
				for( int p = 0; p < ngf; p++ )
				{
					if( p == 0 )
						fprintf(psFile, " %lf %lf moveto\n", pts_out[3*p+0]*scale_factor, pts_out[3*p+1] * scale_factor );
					else
						fprintf(psFile, " %lf %lf lineto\n", pts_out[3*p+0]*scale_factor, pts_out[3*p+1] * scale_factor );				
				}
		
				fprintf(psFile, "closepath\n");
				fprintf(psFile, "fill\n");
				fprintf(psFile, "newpath\n");
				//fprintf(psFile, "%lf %lf %lf setrgbcolor\n", colors[at_link[a]*3+0], colors[at_link[a]*3+1], colors[at_link[a]*3+2] );
				fprintf(psFile, "%lf %lf %lf setrgbcolor\n", use_color[0], use_color[1], use_color[2] );
		
				for( int p = 0; p < ngf; p++ )
				{
					if( p == 0 )
						fprintf(psFile, " %lf %lf moveto\n", pts_out[3*p+0]*scale_factor, pts_out[3*p+1] * scale_factor );
					else
						fprintf(psFile, " %lf %lf lineto\n", pts_out[3*p+0]*scale_factor, pts_out[3*p+1] * scale_factor );				
				}
		
				fprintf(psFile, "closepath\n");
				fprintf(psFile, "0. 0. 0. setrgbcolor\n");
				fprintf(psFile, "1 setlinewidth\n");
				fprintf(psFile, "stroke\n");
		
				free(pts_out);
				free(pts);
			}
		}
	}	
	
	if( doPS )
	{
		for( int p = 0; p < nat; p++ )
		{
			fprintf(psFile, "newpath\n");
			fprintf(psFile, "0. 0. 0. setrgbcolor\n");
			fprintf(psFile, "%lf %lf 3.0 0 360 arc\n", atoms[2*p+0]*scale_factor, atoms[2*p+1] * scale_factor );
			fprintf(psFile, "closepath\nfill\n");
		
			if( solvationShell )
			{
				fprintf(psFile, "%lf %lf moveto\n", atoms[2*p+0]*scale_factor-3, atoms[2*p+1] * scale_factor-15 );
				fprintf(psFile, "/Helvetica findfont\n"  );
				fprintf(psFile, "12.0 scalefont\n"  );
				fprintf(psFile, "setfont\n"  );
				fprintf(psFile, "(%d) show\n", solvationShell[at_link[p]] );
				
				fprintf(psFile, "%lf %lf moveto\n", atoms[2*p+0]*scale_factor-3, atoms[2*p+1] * scale_factor+3 );
				fprintf(psFile, "/Helvetica findfont\n"  );
				fprintf(psFile, "12.0 scalefont\n"  );
				fprintf(psFile, "setfont\n"  );
				fprintf(psFile, "(%d) show\n", at_link[p] );
			}
	
		}

	
		fprintf(psFile, "showpage\n");
		fclose(psFile);
	}

	free(atoms);
	free(at_link);
}

/*
        giftwrap finds the convex hull of points in XY.  I got the algorithm from wikipedia.
*/

void giftwrap( double *pts_in, double *pts_out, int npts, int *nptsout )
{
	double xmin = 1e10;

	int start_pt = -1;

	for( int p = 0; p < npts; p++ )
	{
		if( pts_in[3*p+0] < xmin )
		{
			xmin = pts_in[3*p+0];
			start_pt = p;
		}
	}

	int newpts = 1;

	pts_out[0] = pts_in[3*start_pt+0];
	pts_out[1] = pts_in[3*start_pt+1];
	pts_out[2] = pts_in[3*start_pt+2];

	int use_pt = start_pt;
	double cur_pt[3] = { pts_out[0], pts_out[1], pts_out[2] };
	int next_pt = -1;


	double prev_v[3] = {0,1};

	int pthit[npts];

	memset( pthit, 0, sizeof(int) * npts );

	pthit[start_pt] = 1;

	int is_start_pt = 0;

	do {
		double min_dp= -1e10;

		for( int p = 0; p < npts; p++ )
		{
			double dr[3] = { pts_in[3*p+0] - cur_pt[0],
					 pts_in[3*p+1] - cur_pt[1],
					 pts_in[3*p+2] - cur_pt[2] };

			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2] );
			dr[0] /= r;
			dr[1] /= r;
			dr[2] /= r;

			if( r < 1e-6 )
				continue;

			double dp = prev_v[0] * dr[0] + prev_v[1] * dr[1]; 

			if( dp > min_dp )
			{
				min_dp = dp;
				next_pt = p;
			}
		}

		if( newpts >= npts )
		{
			exit(1);
		}
	
//		fprintf(orig_stdout, "from %lf %lf %lf to %lf %lf %lf from %d to %d.\n",	
//				cur_pt[0], cur_pt[1], cur_pt[2], 
//				pts_in[3*next_pt+0],
//				pts_in[3*next_pt+1],pts_in[3*next_pt+2],	
//				use_pt, next_pt );

		prev_v[0] = pts_in[3*next_pt+0] - cur_pt[0];
		prev_v[1] = pts_in[3*next_pt+1] - cur_pt[1];
		prev_v[2] = pts_in[3*next_pt+2] - cur_pt[2];

		double lpv = sqrt(prev_v[0]*prev_v[0]+prev_v[1]*prev_v[1] +prev_v[2]*prev_v[2]);
		prev_v[0] /= lpv;
		prev_v[1] /= lpv;
		prev_v[2] /= lpv;


		is_start_pt = 0;

		double dr[2] = { pts_in[3*next_pt+0] - pts_in[start_pt*3+0],
				 pts_in[3*next_pt+1] - pts_in[start_pt*3+1] };
			double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
		if( l < 1e-6 )
			is_start_pt = 1;

		if( next_pt != start_pt )
		{
			pts_out[newpts*3+0] = pts_in[3*next_pt+0];
			pts_out[newpts*3+1] = pts_in[3*next_pt+1];
			pts_out[newpts*3+2] = pts_in[3*next_pt+2];
			newpts++;
		}

		use_pt = next_pt;

		cur_pt[0] = pts_in[3*next_pt+0];
		cur_pt[1] = pts_in[3*next_pt+1];
		cur_pt[2] = pts_in[3*next_pt+2];

		if( pthit[next_pt] && (start_pt!=next_pt) && !is_start_pt )
		{
			exit(1);
		}
		pthit[start_pt] = 1;

	} while (next_pt != start_pt && !is_start_pt);


	*nptsout = newpts;
}


/*
	triangle_area gives the area of a triangle in three dimensions.
*/

double triangle_area( double *pt1, double *pt2, double *pt3 )
{
	double dr1[3] = { pt1[0] - pt2[0], pt1[1] - pt2[1], pt1[2] - pt2[2] };
	double l1 = sqrt(dr1[0]*dr1[0]+dr1[1]*dr1[1]+dr1[2]*dr1[2]);

	if( fabs(l1) < 1e-10 ) return 0;

	double dr2[3] = { pt3[0] - pt2[0], pt3[1] - pt2[1], pt3[2] - pt2[2] };

	double nv1[3] = { dr1[0] / l1, dr1[1] / l1, dr1[2] / l1 };

	double dp = dr2[0] * nv1[0] + dr2[1] * nv1[1] + dr2[2] * nv1[2];

	dr2[0] -= dp * nv1[0];
	dr2[1] -= dp * nv1[1];
	dr2[2] -= dp * nv1[2];

	double l2 = sqrt(dr2[0]*dr2[0]+dr2[1]*dr2[1]+dr2[2]*dr2[2]);

	if( !(l1*l2 < 0) && !(l1*l2 > -1) )
	{
		printf("nan\n");
		exit(1);
	}

	return 0.5 * l1 * l2;
}

/*
	loopArea gives the area of a (co-planar) set of points describing a convex polygon, sequentially connected.
*/

double loopArea( double *pts, int npts )
{
	double area = 0;

	for( int x = npts-1; x >= 2; x-- )
		area += triangle_area( pts+0, pts+3*x, pts + 3*(x-1) );
	return area;
}
