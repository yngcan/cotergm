#include "CD.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
 *****************/

/*****************
 void MCMC_wrapper

 Wrapper for a call from R.

 and don't forget that tail -> head
 *****************/
void CD_wrapper(int *dnumnets, int *nedges,
		int *tails, int *heads,
		int *dn,
		int *y0tails, int *y0heads,int *y0dn,int *y0nedges,
		double *y0nodalstatus,
		int *dflag, int *bipartite,
		int *nterms, char **funnames,
		char **sonames,
		char **MHproposaltype, char **MHproposalpackage,
		double *inputs, double *nodalstatus, double *theta0, int *samplesize, int *nsteps,
		double *sample,
		int *fVerbose,
		int *attribs, int *maxout, int *maxin, int *minout,
		int *minin, int *condAllDegExact, int *attriblength,
		int *status){
	int directed_flag;
	Vertex n_nodes, bip, *undotail, *undohead;
	Edge n_networks;
	Network nw[1],y0[1];
	Model *m;
	MHproposal MH;

	n_nodes = (Vertex)*dn;
	Vertex y0n_nodes =  (Vertex)*y0dn;

	n_networks = (Edge)*dnumnets;
	bip = (Vertex)*bipartite;

	GetRNGstate();  /* R function enabling uniform RNG */

	directed_flag = *dflag;

	m=ModelInitialize(*funnames, *sonames, &inputs, *nterms, n_nodes, nodalstatus);

	/* Form the network */
	nw[0]=NetworkInitialize(tails, heads, nedges[0],
			n_nodes, directed_flag, bip, 0, 0, NULL, nodalstatus);

	y0[0]=NetworkInitialize(y0tails, y0heads, y0nedges[0],
			y0n_nodes, directed_flag, bip, 0, 0, NULL, nodalstatus);

	MH_init(&MH,
			*MHproposaltype, *MHproposalpackage,
			inputs,
			*fVerbose,
			nw, attribs, maxout, maxin, minout, minin,
			*condAllDegExact, *attriblength);

	undotail = calloc(MH.ntoggles * *nsteps, sizeof(Vertex));
	undohead = calloc(MH.ntoggles * *nsteps, sizeof(Vertex));

	*status = CDSample(&MH,
			theta0, sample, *samplesize, *nsteps, undotail, undohead,
			*fVerbose, nw, m,y0);

	free(undotail);
	free(undohead);
	MH_free(&MH);

	ModelDestroy(m);
	NetworkDestroy(nw);
	NetworkDestroy(y0);
	PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 MCMCStatus MCMCSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
 *********************/
MCMCStatus CDSample(MHproposal *MHp,
		double *theta, double *networkstatistics, 
		int samplesize, int nsteps, Vertex *undotail, Vertex *undohead, int fVerbose,
		Network *nwp, Model *m, Network *y0){
	int i, j;

	/*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
	 *********************/
	/*for (j=0; j < m->n_stats; j++) */
	/*  networkstatistics[j] = 0.0; */
	/* Rprintf("\n"); */
	/* for (j=0; j < m->n_stats; j++){ */
	/*   Rprintf("j %d %f\n",j,networkstatistics[j]); */
	/* } */
	/* Rprintf("\n"); */

	/* Now sample networks */
	for (i=0; i < samplesize; i++){

		if(CDStep(MHp, theta, networkstatistics, nsteps, undotail, undohead,
				fVerbose, nwp, m, y0)!=MCMC_OK)
			return MCMC_MH_FAILED;

#ifdef Win32
		if( ((100*i) % samplesize)==0 && samplesize > 500){
			R_FlushConsole();
			R_ProcessEvents();
		}
#endif
		networkstatistics += m->n_stats;
	}

	return MCMC_OK;
}

/*********************
 void MetropolisHastings

 In this function, theta is a m->n_stats-vector just as in MCMCSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
 *********************/
MCMCStatus CDStep(MHproposal *MHp,
		double *theta, double *networkstatistics,
		int nsteps, Vertex *undotail, Vertex *undohead,
		int fVerbose,
		Network *nwp,
		Model *m, Network *y0) {

	unsigned int unsuccessful=0, ntoggled=0;

	for(unsigned int step=0; step<nsteps; step++){
		MHp->logratio = 0;
		(*(MHp->func))(MHp, nwp); /* Call MH function to propose toggles */

		if(MHp->toggletail[0]==MH_FAILED){
			if(MHp->togglehead[0]==MH_UNRECOVERABLE)
				error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
			if(MHp->togglehead[0]==MH_IMPOSSIBLE){
				Rprintf("MH Proposal function encountered a configuration from which no toggle(s) can be proposed.\n");
				return MCMC_MH_FAILED;
			}
			if(MHp->togglehead[0]==MH_UNSUCCESSFUL){
				warning("MH Proposal function failed to find a valid proposal.");
				unsuccessful++;
				if(unsuccessful>MH_QUIT_UNSUCCESSFUL){
					Rprintf("Too many MH Proposal function failures.\n");
					return MCMC_MH_FAILED;
				}

				return MH_UNSUCCESSFUL;
			}
		}

		if(fVerbose>=5){
			Rprintf("Proposal: ");
			for(unsigned int i=0; i<MHp->ntoggles; i++)
				Rprintf(" (%d, %d)", MHp->toggletail[i], MHp->togglehead[i]);
			Rprintf("\n");
		}

		/* Calculate change statistics,
       remembering that tail -> head */
		ChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->togglenode,nwp, m, y0);

		if(fVerbose>=5){
			Rprintf("Changes: (");
			for(unsigned int i=0; i<m->n_stats; i++)
				Rprintf(" %f ", m->workspace[i]);
			Rprintf(")\n");
		}

		/* Calculate inner product */
		double ip=0;
		for (unsigned int i=0; i<m->n_stats; i++){
			ip += theta[i] * m->workspace[i];
		}
		/* The logic is to set cutoff = ip+logratio ,
       then let the MH probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
		double cutoff = ip + MHp->logratio;

		if(fVerbose>=5){
			Rprintf("log acceptance probability: %f + %f = %f\n", ip, MHp->logratio, cutoff);
		}

		/* if we accept the proposed network */
		if (cutoff >= 0.0 || log(unif_rand()) < cutoff) {
			if(fVerbose>=5){
				Rprintf("Accepted.\n");
			}

			if(step<nsteps-1){
				/* Make proposed toggles (updating timestamps--i.e., for real this time) */
				for(unsigned int i=0; i < MHp->ntoggles; i++){
					undotail[ntoggled]=MHp->toggletail[i];
					undohead[ntoggled]=MHp->togglehead[i];
					ntoggled++;
					ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);

					if(MHp->discord)
						for(Network **nwd=MHp->discord; *nwd!=NULL; nwd++){
							ToggleEdge(MHp->toggletail[i],  MHp->togglehead[i], *nwd);
						}
				}
			}

			/* record network statistics for posterity */
			for (unsigned int i = 0; i < m->n_stats; i++){
				networkstatistics[i] += m->workspace[i];
			}

		}else{
			if(fVerbose>=5){
				Rprintf("Rejected.\n");
			}
		}
	}

	/* Undo toggles. */
	for(unsigned int i=0; i < ntoggled; i++){
		Vertex t = undotail[i], h = undohead[i];
		ToggleEdge(t, h, nwp);

		if(MHp->discord)
			for(Network **nwd=MHp->discord; *nwd!=NULL; nwd++){
				ToggleEdge(t, h, *nwd);
			}
	}

	return MCMC_OK;
}

