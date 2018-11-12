#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "binary.h"


int main ( int argc, char **argv )
{
	if( argc < 2 )
	{
		printf("Syntax: postprocess binary.output [OPTIONS]\n");
		printf("---------------------------------\n");
		printf("options: cutoff=[FLOAT] default 0.0\n");
		printf("options: op_cut=[FLOAT] default 0.7\n");
		printf("options: av_window=[INT] default 1, but 10 seems good\n");
		printf("options: thresh_hit [default off]\n");
		printf("\t\tthresh_hit doesn't count a border until it hits\n");
		printf("\t\tthe cutoff, then it is counted both\n");
		printf("\t\tin the future and the past.\n");
		printf("\t\tThis option requires the entire binary file to be read into memory\n");
		printf("\t\te.g., it may require circa one gigabyte of memory.\n");
		printf("options: cluster_save=[FILENAME]\n");
		printf("\t\tsave cluster assignments to a binary file.\n");
		printf("options: min_cluster_size=[INT] default: 3\n");
		return 0;
	}
	
	double cutoff = 0.0;
	double op_cut = 0.7;	
	int av_window = 1;
	int do_forward_backward = 0;
	double fb_cut = 0;
	char cluster_file[256];
	int save_cluster_file = 0;
	int min_cluster_size = 3;

	for( int c = 2; c < argc; c++ )	
	{
		if( !strncasecmp( argv[c], "CUTOFF=", 7 ) && strlen(argv[c]) > 7 )
			cutoff = atof(argv[c]+7);
		else if( !strncasecmp( argv[c], "op_cut=", 7 ) && strlen(argv[c]) > 7 )
			op_cut = atof(argv[c]+7);
		else if( !strncasecmp( argv[c], "av_window=", 10 ) && strlen(argv[c]) > 10 )
			av_window = atoi(argv[c]+10);
		else if( !strncasecmp( argv[c], "min_cluster_size=", 17 ) && strlen(argv[c]) > 17 )
			min_cluster_size = atoi(argv[c]+17);
		else if( !strncasecmp( argv[c], "cluster_save=", 13) && strlen(argv[c]) > 13 )
		{
			strcpy( cluster_file, argv[c]+13 );
			save_cluster_file = 1;
		}
		else if( !strcasecmp( argv[c], "thresh_hit" ) )
			do_forward_backward = 1;
		else
		{
			printf("Couldn't interpret option '%s'.\n", argv[c] );
			exit(1);
		}
	}
//	printf("cutoff: %lf\n", cutoff );
	FILE *theFile = fopen( argv[1], "rb" );

	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", argv[1] );
		return 0;
	}

	binaryHeader header;

	fread( &header, 1, sizeof(binaryHeader), theFile );	

	if( header.formatCheck != FORMAT_CHECK )
	{
		printf("The file header did not pass a format check.\n");
		return 0;
	}

	int nind = header.nind;
	int nframes = header.nframes;

	elementDescriptor desc[nind];

	fread( desc, nind, sizeof(elementDescriptor), theFile );

	binaryElement readFrame[nind];

	double av_contacts[nind];
	memset( av_contacts, 0, sizeof(double) * nind );
	double av_len[nind];
	memset( av_len, 0, sizeof(double) * nind );

#define N_OP_BINS 100
#define CL_MAX 200
	
	double *sd = (double *)malloc( sizeof(double) * N_OP_BINS );
	memset( sd, 0, sizeof(double) * N_OP_BINS );

	double *op = (double *)malloc( sizeof(double) * N_OP_BINS );
	memset( op, 0, sizeof(double) * N_OP_BINS );
	
	double *op2 = (double *)malloc( sizeof(double) * N_OP_BINS );
	memset( op2, 0, sizeof(double) * N_OP_BINS );

	double *cl_hist = (double *)malloc( sizeof(double) * CL_MAX );
	memset( cl_hist, 0, sizeof(double) * CL_MAX );

	double wsize = sizeof(float) * nind * av_window / (1024.0)/(1024.0);

	if( wsize > 1000 )
	{
		printf("windowing array greater than one gigabyte.  Remove this code if you'd still like to run with this window.\n");
		exit(1);
	}

	double *op_window = (double *)malloc( sizeof(double) * nind * av_window );
	double *cur_op = (double *)malloc( sizeof(double) * nind );
	memset( cur_op, 0, sizeof(double) * nind );	

	binaryElement *wholeTraj = NULL;
	
	int pos = ftell(theFile);

	double *areas = (double *)malloc( sizeof(double) * 10000 );
	double *areas2 = (double *)malloc( sizeof(double) * 10000 );
	double *narea = (double *)malloc( sizeof(double) * 10000 );
	char *names[100];
	for( int t = 0; t < 100; t++ )
		names[t] = NULL;
	memset( areas, 0, sizeof(double) * 10000 );
	memset( areas2, 0, sizeof(double) * 10000 );
	memset( narea, 0, sizeof(double) * 10000 );

	
	
	for( int f = 0; f < nframes; f++)
	{
		int formatCheck = 0;
		int nr = fread( &formatCheck, sizeof(int), 1, theFile );
		wholeTraj = (binaryElement *)malloc( sizeof(binaryElement) * (nframes * nind) );
		nr = fread( wholeTraj+nind*f, sizeof(binaryElement), nind, theFile );
		
		for( int x = 0; x < nind; x++ )
		{
			int type = wholeTraj[nind*f+x].type;
	
			if( !names[type] ) 
			{
				names[type] = (char *)malloc( sizeof(char) * ( 1 +strlen(desc[x].resname) ) );
				strcpy( names[type], desc[x].resname );
			}
	
			double val = wholeTraj[nind*f+x].area;


			areas[type*2+desc[x].chain] += val;
			areas2[type*2+desc[x].chain] += val*val;
			narea[type*2+desc[x].chain] += 1;

//			areas[wholeTraj[nind*f+x].type] += 
//			printf("%s %lf %d\n", desc[x].resname, 
//					    wholeTraj[nind*f+x].area, desc[x], wholeTraj[nind*f+x].type  ); 
		}
		

		if( nr != nind )
		{
			printf("Failed to read frame '%d/%d' (counting from 1).\n", f+1, nframes );
			return 0;
		}
		free(wholeTraj);
	}

	for( int type = 0; type < 100; type++ )
	{
		if( names[type] && narea[2*type] > 0 )
		{
			double av = areas[2*type]/narea[2*type];
			double av2 = areas2[2*type]/narea[2*type];
			double std = sqrt(av2-av*av);
		
			if( narea[2*type+1] > 0 )
			{
				double av_2 = areas[2*type+1]/narea[2*type+1];
				double av2_2 = areas2[2*type+1]/narea[2*type+1];
				double std_2 = sqrt(av2_2-av_2*av_2);
				printf("%s area1 %lf stdev1 %lf area2 %lf stdev2 %lf\n", names[type], av, std, av_2, std_2 );
			}
			else
				printf("%s area %lf stdev %lf area2 N/A stdev2 N/A\n", names[type], av, std );

		}
			
	}

	fseek( theFile, pos, SEEK_SET );

	if( do_forward_backward )
	{
		int formatCheck = 0;
		wholeTraj = (binaryElement *)malloc( sizeof(binaryElement) * (nframes * nind) );

		for( int f = 0; f < nframes; f++)
		{
			int nr = fread( &formatCheck, sizeof(int), 1, theFile );
			nr = fread( wholeTraj+nind*f, sizeof(binaryElement), nind, theFile );
			if( nr != nind )
			{
				printf("Failed to read frame '%d/%d' (counting from 1).\n", f+1, nframes );
				return 0;
			}
		}

		double thresh = cutoff;
		int skipBacktrack = 0;
		int skipForetrack = 0;
		for( int t_master = 0; t_master < nframes; t_master++ )
		{
			for( int i = 0; i < nind; i++ )
			{
				for( int ix = 0; ix < wholeTraj[t_master*nind+i].nborders; ix++ )
				{
					int j = wholeTraj[t_master*nind+i].borderID[ix];
	
					double val = wholeTraj[t_master*nind+i].borderLen[ix];
	
					if( val > thresh )
					{
						wholeTraj[t_master*nind+i].borderLen[ix] = -1;
	
						if( skipBacktrack == 0 )
						{
							int done = 0;
							for( int t_backtrack = t_master-1; !done && t_backtrack >= 0; t_backtrack-- )
							{
								int got_j = 0;
								for( int ix2 = 0; !done && ix2 < wholeTraj[t_backtrack*nind+i].nborders; ix2++ )
								{
									if( wholeTraj[t_backtrack*nind+i].borderID[ix2] == j )
									{
										got_j = 1;
										if( wholeTraj[t_backtrack*nind+i].borderLen[ix2] >= 0 )
											wholeTraj[t_backtrack*nind+i].borderLen[ix2] = -1;
										else 
											done = 1;
									}	
								}
								if( !got_j ) done = 1;
							}
						}
					}
					else if( t_master > 0 && skipForetrack == 0 )
					{
						int t_backtrack = t_master-1;
						for( int ix2 = 0; ix2 < wholeTraj[t_backtrack*nind+i].nborders; ix2++ )
						{
							if( wholeTraj[t_backtrack*nind+i].borderID[ix2] == j )
							{
								if( wholeTraj[t_backtrack*nind+i].borderLen[ix2] < 0  )
									wholeTraj[t_master*nind+i].borderLen[ix] = -1;
							}	
						}
					} 
				}
			}
		}

		for( int f = 0; f < nframes; f++ )
		{
			for( int x = 0; x < nind; x++ )
			{
				for( int bx = 0; bx < wholeTraj[f*nind+x].nborders; bx++ )
				{
					if( wholeTraj[f*nind+x].borderLen[bx] < 0 )
						wholeTraj[f*nind+x].borderLen[bx] = cutoff + 1;
					else
						wholeTraj[f*nind+x].borderLen[bx] = 0;
				}
			}
		}	
	}	

	/* we would like the cluster assignments to make sense -- i.e., your previous cluster index
	   and the new cluster index should be the same if the cluster hasn't changed much.
	   clusters have a majority previous cluster, which they claim.
	   a cluster which has the most of a previous cluster has the best claim.
	   if you do not have the best claim to a previous cluster, you are assigned
	   a completely new index.  */

	int *prev_cluster = (int *)malloc( sizeof(int) * nind );
	int *cur_cluster = (int *)malloc( sizeof(int) * nind );
	for( int x = 0; x < nind; x++ )
		prev_cluster[x] = -1;

	FILE *clusterFile = NULL;

	if( save_cluster_file )
		clusterFile = fopen( cluster_file, "wb" );

	int cur_claim = 0;

	double *curSD = (double *)malloc( sizeof(double) * nind );
		int *clusters = (int *)malloc( sizeof(int) * nind * nind );
		int *clusterSize = (int *)malloc( sizeof(int) * nind );

	for( int f = 0; f < nframes; f++)
	{
		int formatCheck = 0;
		
		if( do_forward_backward )
		{
			memcpy( readFrame, wholeTraj + f * nind , sizeof(binaryElement) * nind);	
		}
		else
		{
			int nr = fread( &formatCheck, sizeof(int), 1, theFile );

			if( nr != 1 || formatCheck != FORMAT_CHECK )
			{
				printf("The file header did not pass a format check at frame %d/%d (counting from 1).\n", f+1, nframes);
				return 0;
			}

			nr = fread( readFrame, sizeof(binaryElement), nind, theFile );

			if( nr != nind )
			{
				printf("Failed to read frame '%d/%d' (counting from 1).\n", f+1, nframes );
				return 0;
			}
		}


		/* calculate the order parameter */

		for( int x = 0; x < nind; x++ )
		{
			double orderpc = 0;
			double orderps = 0;
			double orderpn = 0;

			for( int b = 0; b < readFrame[x].nborders; b++ )
			{
				if( readFrame[x].borderLen[b] < cutoff ) continue; 

				double theta = readFrame[x].borderTheta[b];

				orderpc += cos( 6 * theta ); 
				orderps += sin( 6 * theta );
				orderpn += 1; 
			}

			double val = orderpc * orderpc + orderps * orderps;

			val = fabs(val);
			val = sqrt(val);

//			printf("%d %lf %lf\n", x, readFrame[x].op, val / orderpn );

//			printf("orderpn: %lf val: %lf\n", orderpn, val );

			if( orderpn > 0 )
			{
				val /= orderpn;

				readFrame[x].op = val;
			}
			else
				readFrame[x].op = -1;


			op_window[x*av_window+(f%av_window)] = readFrame[x].op;
		}
		int lim = av_window;

		if( f < av_window ) lim = f+1;
		
		for( int x = 0; x < nind; x++ )
		{
			double av = 0; 
			double av2 = 0;
			double nav = 0;
			for( int bx = 0; bx < readFrame[x].nborders; bx++ )
			{
				av += readFrame[x].borderLen[bx];
				av2 += readFrame[x].borderLen[bx] * readFrame[x].borderLen[bx];
				nav+=1;
			}

			av2 /= nav;
			av /= nav;

			curSD[x] = sqrt(av2-av*av);

			int sdb = N_OP_BINS * curSD[x] /3.0;

			if( sdb < N_OP_BINS )
				sd[sdb] += 1;

//			printf("sd %lf\n", curSD[x] ); 
		}


		for( int x = 0; x < nind; x++ )
		{
			cur_op[x] = 0;
			for( int i = 0; i < lim; i++ )
				cur_op[x] += op_window[x*av_window+i]; 
			cur_op[x] /= lim;
				

//			cur_op[x] *= exp(-curSD[x]);	
//			if( readFrame[x].type )
//				cur_op[x] = 0;
		}


		memset( clusterSize, 0, sizeof(int) * nind );

		for( int x = 0; x < nind; x++ )
		{
			readFrame[x].cluster = -1;

			if( cur_op[x] >= op_cut )
//			if( readFrame[x].op >= op_cut )
			{
//				printf("assigning %d!\n", x );
				readFrame[x].cluster = x;
				clusters[x*nind+0] = x;
				clusterSize[x] = 1;
			}			
		}

		for( int x = 0; x < nind; x++ )
		{
//			printf("x: %d cluster 438 size: %d\n", x, clusterSize[438] );
			if( cur_op[x] < op_cut ) continue;
//			if( readFrame[x].op < op_cut ) continue;

			for( int b = 0; b < readFrame[x].nborders; b++ )
			{
				if( readFrame[x].borderLen[b] > cutoff )
				{
					int y = readFrame[x].borderID[b]; 
			
					if( cur_op[y] < op_cut ) continue;
				//	if( readFrame[y].op < op_cut ) continue;

					int add_cluster = readFrame[y].cluster;
					int my_cluster = readFrame[x].cluster;
	
					if( my_cluster != add_cluster )
					{
						for( int p = 0; p < clusterSize[my_cluster]; p++ )
						{
//							if( clusters[my_cluster*nind+p] == 54 ) printf("Moving 54 to %d!\n", add_cluster );
							clusters[add_cluster*nind+clusterSize[add_cluster]] = clusters[my_cluster*nind+p];
							readFrame[clusters[my_cluster*nind+p]].cluster = add_cluster;
							clusterSize[add_cluster]++;
						}
	
						clusterSize[my_cluster] = 0;
					}
				}
			}
		}


		int nclusters = 0;

		for( int x = 0; x < nind; x++ )
		{
			cur_cluster[x] = -1;
			if( clusterSize[x] >= min_cluster_size) 
				nclusters++;
		}

		if( nclusters > 0 )
		{
			int cluster_id_claim[nclusters];
			int cluster_id_claim_second[nclusters];
			int cluster_n_claim[nclusters];
			int cluster_n_claim_second[nclusters];
			int cluster_map[nclusters];
			int curc = 0;

			int maxv = 0;

			for( int x = 0; x < nind; x++ )
			{
				if( prev_cluster[x] >= maxv ) maxv = prev_cluster[x]+1;
			}

			for( int x = 0; x < nind; x++ )
			{
				if( clusterSize[x] >= min_cluster_size)
				{
					if( clusterSize[x] < CL_MAX )
						cl_hist[clusterSize[x]] += 1;
					else
						clusterSize[CL_MAX-1] += 1;

					int id_claim[maxv];
					memset( id_claim, 0, sizeof(int) * maxv );

					for( int p = 0; p < clusterSize[x]; p++ )
					{
						if( prev_cluster[clusters[x*nind+p]] >= 0 )
							id_claim[prev_cluster[clusters[x*nind+p]]] += 1;
					}

					int best_n = 0;
					int cur_best = -1;
					int best_n2 = 0;
					int cur_best2 = -1;

					for( int p = 0; p < maxv; p++ )
					{
						if( id_claim[p] > best_n  )
						{
							cur_best2 = cur_best;
							best_n2 = best_n;

							cur_best = p;
							best_n = id_claim[p];
						}
					}

					cluster_id_claim[curc] = cur_best;	
					cluster_n_claim[curc] = best_n;
					
					cluster_id_claim_second[curc] = cur_best2;	
					cluster_n_claim_second[curc] = best_n2;

					cluster_map[curc] = x;

					curc++;
				}	
			}

			int max_size = 0;

			for( int m = 0; m < nclusters; m++ )
			{
				if( clusterSize[cluster_map[m]] > max_size ) max_size = clusterSize[cluster_map[m]];

				int claim = cluster_id_claim[m];		
				int nclaim = cluster_n_claim[m];

				int yes = 0;

				if( claim >= 0 )
				{
					yes = 1;
					for( int m2 = 0; m2 < nclusters; m2++ )
					{
						if( m2 == m ) continue;
	
						if( cluster_id_claim[m2] == claim && 
						    cluster_n_claim[m2] > nclaim )
							yes = 0;
						else if( cluster_id_claim[m2] == claim && 
						    cluster_n_claim[m2] == nclaim  && m2 < m )
							yes = 0;
					}	
				}

				if( yes )
				{
					int cl = cluster_map[m];
					for( int p = 0; p < clusterSize[cl]; p++ )
						cur_cluster[clusters[cl*nind+p]] = claim;
				}
				else
				{
					int cl = cluster_map[m];
					for( int p = 0; p < clusterSize[cl]; p++ )
						cur_cluster[clusters[cl*nind+p]] = cur_claim;
					cur_claim++;
//					printf("cur_claim: %d\n", cur_claim );
				}
			}		
//			printf("nclusters: %d max: %d\n", nclusters, max_size );	
		}
		for( int x = 0; x < nind; x++ )
		{
/*
			if( readFrame[x].op > 0.9 )
				cur_cluster[x] = 1;
			else
				cur_cluster[x] = -1;*/

			if( cur_cluster[x] >= 0 )
				prev_cluster[x] = cur_cluster[x];
		}


		if( clusterFile )
		{
			for( int x = 0; x < nind; x++ )
			{
//				if( cur_cluster[x] >= 0 ) printf("res %d cluster %d.\n", x, cur_cluster[x]);
			}
//			for( int x = 0; x < nind; x++ )
//				printf(" %d ", cur_cluster[x] );
//			printf("\n");
			fwrite( cur_cluster, sizeof(int), nind, clusterFile );
		}

		for( int x = 0; x < nind; x++ )		
		{
			if( readFrame[x].op >= 0 )
			{
				int opb = cur_op[x] * N_OP_BINS;

				if( opb >= 0 && opb < N_OP_BINS )
				{
					if( readFrame[x].type )
						op[opb] += 1;
					else
						op2[opb] += 1;
				}
			}
	
			for( int b = 0; b < readFrame[x].nborders; b++ )
			{
				if( readFrame[x].borderLen[b] > cutoff )
				{
					av_contacts[x] += 1;
					av_len[x] += readFrame[x].borderLen[b];
				}
			}
		}
	}
/*
	printf("OP HISTOGRAM\n");
	for( int b = 0; b < N_OP_BINS; b++ )
	{
		printf("%lf %lf %lf\n", (b+0.5) / (double)N_OP_BINS, op[b], op2[b] );
	}
	printf("CLUSTER HISTOGRAM\n");
	for( int b = 0; b < CL_MAX; b++ )
	{
		printf("%d %lf\n", b, cl_hist[b] );
	}
*/	

#if 0
	for( int x = 0; x < nind; x++ )
	{
		printf("Residue %s:%s residue %d chain %d average number of contacts: %lf average length: %lf\n",  
			desc[x].resname, desc[x].segid, desc[x].res, desc[x].chain,
			av_contacts[x] / nframes, av_len[x] / (av_contacts[x] ) );
	}
#endif
}









