#ifndef NETSTATS_H
#define NETSTATS_H

#include "edgetree.h"
#include "model.h"
#include "MHproposal.h"

/* *** don't forget tail -> head, so these functions accept tails first, not heads */

void network_stats_wrapper(int *tails, int *heads, int *timings, int *time, int *lasttoggle, int *dnedges, 
			   int *dn,
				int *y0tails, int *y0heads,int *y0dn,int *y0dnedges,double *y0nodalstatus,

			   int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs, double *nodalstatus,   double *stats);
void SummStats(Edge n_edges, Vertex n_nodes, Vertex *tails, Vertex *heads,
	       Network *nwp, Model *m, double *stats, Network *y0);
#endif
