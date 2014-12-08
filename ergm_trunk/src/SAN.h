#ifndef SAN_H
#define SAN_H

#include "edgetree.h"
#include "changestat.h"
#include "MHproposal.h"
#include "model.h"
#include "MCMC.h"

/* *** don't forget tail -> head, so this function accepts tails first, not heads  */


void SAN_wrapper (int *dnumnets, int *nedges,
		  int *tails, int *heads,
		  int *dn,
     		int *y0tails, int *y0heads,int *y0dn,int *y0nedges,double *y0nodalstatus,

		  int *dflag, int *bipartite,
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHproposaltype, char **MHproposalpackage,
		  double *inputs, double *nodalstatus, double *theta0, double *tau, int *samplesize,
		  double *sample, int *burnin, int *interval,  
		  int *newnetworktails, 
		  int *newnetworkheads, 
		  int *newnodalstatus,
		  double *invcov,
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		  int *maxedges,
		  int *status);

MCMCStatus SANSample (MHproposal *MHp,
		double *theta, double *invcov, double *tau, double *networkstatistics, 
		int samplesize, int burnin, 
		int interval, int fVerbose, int nmax,
		Network *nwp, Model *m, Network *y0);
MCMCStatus SANMetropolisHastings (MHproposal *MHp,
			 double *theta, double *invcov, double *tau, double *statistics, 
			 int nsteps, int *staken,
			 int fVerbose,
			 Network *nwp, Model *m, Network *y0);
#endif
