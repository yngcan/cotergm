#include "netstats.h"
#include "edgetree.h"

/*****************
 void network_stats_wrapper

 Wrapper for a call from R.  Return the change in the statistics when
 we go from an empty graph to the observed graph.  If the empty graph
 has true global values equal to zero for all statistics, then this
 change gives the true global values for the observed graph.
*****************/

/* *** don't forget tail-> head, so this fucntion now accepts tails before heads */

void network_stats_wrapper(int *tails, int *heads, int *timings, int *time, int *lasttoggle, int *dnedges, 
			   int *dn,
				int *y0tails, int *y0heads,int *y0dn,int *y0dnedges,double *y0nodalstatus,

				int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs,   double *nodalstatus, double *stats)
{
  int directed_flag;
  Vertex n_nodes;
  Edge n_edges, y0n_edges;
  Network nw[2];
  Network y0[1];
  Model *m;
  Vertex bip;

  double * negstatus;


/*	     Rprintf("prestart with setup\n"); */
  n_nodes = (Vertex)*dn; 
  Vertex y0n_nodes =  (Vertex)*y0dn;
  n_edges = (Edge)*dnedges;
  y0n_edges = (Edge)*y0dnedges;
  directed_flag = *dflag;
  bip = (Vertex)*bipartite;
  
  negstatus = (double *)malloc( n_nodes * sizeof(double));

  /*set to negative */
  for (int i=0; i<n_nodes; i++)
	  negstatus[i] = 1.0 ;

  if(*lasttoggle == 0) lasttoggle = NULL;

  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms, n_nodes, nodalstatus);

  nw[0]=NetworkInitialize(NULL, NULL, 0,
                          n_nodes, directed_flag, bip, *timings?1:0, *timings?*time:0, *timings?lasttoggle:NULL, negstatus);

  y0[0]=NetworkInitialize(y0tails, y0heads, y0n_edges,
				y0n_nodes, directed_flag, bip, 0, 0, NULL, y0nodalstatus);
  /* Compute the change statistics and copy them to stats for return
     to R.  Note that stats already has the statistics of an empty
     network, so d_??? statistics will add on to them, while s_???
     statistics will simply overwrite them. */
  SummStats(n_edges, n_nodes, tails, heads, nw, m, stats, y0);
  
  ModelDestroy(m);
  NetworkDestroy(nw);
}


/****************
 void SummStats Computes summary statistics for a network. Must be
 passed an empty network (and a possible discordance network) and 
 passed an empty network
*****************/

/* *** don't forget tail-> head, so this fucntion now accepts tails before heads */

void SummStats(Edge n_edges, Vertex n_nodes, Vertex *tails, Vertex *heads,
Network *nwp, Model *m, double *stats, Network *y0){

	Vertex *nodes; /*may need to change*/

  nodes = (Vertex *)malloc( 1 * sizeof(Vertex));

  GetRNGstate();  /* R function enabling uniform RNG */
  
  ShuffleEdges(tails,heads,n_edges); /* Shuffle edgelist. */
  
  for (unsigned int termi=0; termi < m->n_terms; termi++)
    m->termarray[termi].dstats = m->workspace;
  
  /* Doing this one toggle at a time saves a lot of toggles... */
  for(Edge e=0; e<n_edges; e++){
    ModelTerm *mtp = m->termarray;
    double *statspos=stats;
    
    nodes[0] = 0;

    for (unsigned int termi=0; termi < m->n_terms; termi++, mtp++){
      if(!mtp->s_func){
        (*(mtp->d_func))(1, tails+e, heads+e, nodes,
        mtp, nwp, y0);  /* Call d_??? function */
        for (unsigned int i=0; i < mtp->nstats; i++,statspos++)
          *statspos += mtp->dstats[i];
      }else statspos += mtp->nstats;
    }
    ToggleEdge(tails[e],heads[e],nwp);
  }
  

  /* Doing this one toggle at a time saves a lot of toggles... */
  for(int v=0; v<n_nodes; v++){
	  if(m->nodalstatus[v]== 2.0){

		  nodes[0] =  (int) (v+1);

		  ModelTerm *mtp = m->termarray;

		double *statspos=stats;

    for (unsigned int termi=0; termi < m->n_terms; termi++, mtp++){
      if(!mtp->s_func){
        (*(mtp->d_func))(1, tails, heads, nodes,
        mtp, nwp, y0);  /* Call d_??? function */
        for (unsigned int i=0; i < mtp->nstats; i++,statspos++)
          *statspos += mtp->dstats[i];
      }else statspos += mtp->nstats;
    }
    ToggleNode(v+1,nwp);
  }
  }



  ModelTerm *mtp = m->termarray;
  double *dstats = m->workspace;
  double *statspos=stats;
  for (unsigned int termi=0; termi < m->n_terms; termi++, dstats+=mtp->nstats, mtp++ ){
    if(mtp->s_func){
      (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
      for (unsigned int i=0; i < mtp->nstats; i++,statspos++)
        *statspos = mtp->dstats[i];
    }else statspos += mtp->nstats;
  }
  
  PutRNGstate();
}

