#include "wtMHproposals.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_CompleteOrdering

 Default MH algorithm for ERGM over complete orderings
*********************/
void MH_CompleteOrdering(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex tail, head1, head2;
  
  if(MHp->ntoggles == 0) { // Initialize CompleteOrdering 
    MHp->ntoggles=2;
    return;
  }

  unsigned int trytoggle;

  for(trytoggle=0; trytoggle<MAX_TRIES; trytoggle++){
  
    GetRandDyad(&tail, &head1, nwp);
    
    if(nwp->bipartite){
      head2 = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite - 1);
      if(head2 >= head1) head2++;
    }else{
      head2 = 1 + unif_rand() * (nwp->nnodes - 2);
      if(head2 >= tail) head2++;
      if(head2 >= head1){
	head2++;
	if(head2 == tail) head2++; // Increment if landed on top of tail. Note that head2 must have started out < tail, so there is no possibility of head2 being incremented 3 times.
      }
    }
    
    Mtail[0] = Mtail[1] = tail;
    Mhead[0] = head1;
    Mhead[1] = head2;
    
    Mweight[1] = WtGetEdge(Mtail[0],Mhead[0],nwp);
    Mweight[0] = WtGetEdge(Mtail[1],Mhead[1],nwp);

    if(Mweight[0]!=Mweight[1]) break;
  }
  if(trytoggle>=MAX_TRIES){
    MHp->toggletail[0]=MH_FAILED;
    MHp->togglehead[0]=MH_UNSUCCESSFUL;	
  }
}

/*********************
 void MH_CompleteOrdering

 Default MH algorithm for ERGM over complete orderings
*********************/
void MH_CompleteOrderingEquivalent(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex tail, head1, head2;
  static Vertex negos;
  
  // Initialize
  if(MHp->ntoggles == 0) { 
    negos = nwp->bipartite ? nwp->bipartite : nwp->nnodes;
    MHp->ntoggles=2;

    // This really needs to be moved into R code somewhere.
    // The proposal requires that there be no ties, but at the same
    // time, the dyad values must not violate constraints implied by
    // the data.
    // Therefore, these values must somehow be made unique; the order
    // within the equivalence class doesn't matter, but values in
    // different classes must not cross.
    // The way I am doing this is very simple but very dangerous,
    // since it assumes that we will never ever encounter two values
    // that are legitimately different but sufficiently close to each other.

    const double MAX_JITTER = 0.00001;

    for(tail=1;tail<=negos;tail++){
      for(head1=nwp->bipartite? negos+1 : 1; head1<=nwp->nnodes; head1++){
	if(tail==head1) continue;
	WtSetEdge(tail, head1, WtGetEdge(tail, head1, nwp) + unif_rand()*MAX_JITTER, nwp);
      }
    }
    return;
  }
  
  double *clist;

  do{
     tail = 1 + unif_rand() * negos;
     clist = MHp->inputs + (unsigned int) MHp->inputs[tail-1]; // List of equivalence classes for that i, preceded by their number.
  }while(*clist==0);

  unsigned int c = 1 + (unsigned int) (unif_rand() * (*clist)); // Select an equivalence class. TODO: Oversample bigger classes.
  unsigned int jstart = clist[c], jend = clist[c+1]; // Locate the start and the end of the list of j's within this class.
  head1 = clist[(int)(jstart + unif_rand()*(jend-jstart))]; // Select a j1 at random from this class.
  do{
    head2 = clist[(int)(jstart + unif_rand()*(jend-jstart))]; // Keep trying until you get a new j2.
  }while(head2 == head1);
  
  Mtail[0] = Mtail[1] = tail;
  Mhead[0] = head1;
  Mhead[1] = head2;
  
  Mweight[1] = WtGetEdge(Mtail[0],Mhead[0],nwp);
  Mweight[0] = WtGetEdge(Mtail[1],Mhead[1],nwp);
}


/*********************
 void MH_StdNormal

 Default MH algorithm for a standard-normal-reference ERGM
*********************/
void MH_StdNormal(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  
  if(MHp->ntoggles == 0) { // Initialize StdNormal 
    MHp->ntoggles=1;
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double propsd = 0.2; // This ought to be tunable.

  Mweight[0] = rnorm(oldwt, propsd);    
  
  // Symmetric proposal, but depends on the reference measure
  MHp->logratio += -(Mweight[0]*Mweight[0]-oldwt*oldwt)/2;
}

/*********************
 void MH_StdNormalRank

 Default MH algorithm for a standard-normal-reference ERGM
 subject to the constraint that the ranking of alters for
 each ego is preserved.
*********************/
void MH_StdNormalRank(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  
  if(MHp->ntoggles == 0) { // Initialize StdNormalRank
    MHp->ntoggles=1;
    return;
  }
  
  Mtail[0] = 1 + unif_rand() * nwp->nnodes;
  while ((Mhead[0] = 1 + unif_rand() * nwp->nnodes) == Mtail[0]);
 
  // Note that undirected networks do not make sense for this proposal.

  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  // Evaluate rank constraints for this dyad.
  double ub = +HUGE_VAL, lb = -HUGE_VAL;
  int head_class = MHp->inputs[(Mtail[0]-1)*nwp->nnodes + (Mhead[0]-1)];

  for(unsigned int h=1;h<=nwp->nnodes; h++){
    if(h==Mtail[0] || h==Mhead[0]) continue;
    unsigned int h_class=MHp->inputs[(Mtail[0]-1)*nwp->nnodes + (h-1)];

    if(h_class){
      if(h_class==head_class+1){
	// If this alter bounds the current alter from above...
	double b = WtGetEdge(Mtail[0],h,nwp);
	if(b<ub) ub=b;
      }else if(h_class==head_class-1){
	// If this alter bounds the current alter from below...
	double b = WtGetEdge(Mtail[0],h,nwp);
	if(b>lb) lb=b;
      }
    }
  }

  const double propsd = 1; // This ought to be tunable.
  // Note that for a given dyad, which branch gets taken here is always the same.
  if(ub==+HUGE_VAL && lb==-HUGE_VAL){
    // No constraint -> a normal jump
    Mweight[0] = rnorm(oldwt, propsd);
  }else if(ub==+HUGE_VAL){
    // Only lower bound -> a constrained normal jump
    do{ Mweight[0] = rnorm(oldwt, propsd); }while(Mweight[0]<lb);
      
    MHp->logratio += dnorm(oldwt, Mweight[0], propsd, 1) - pnorm(lb, Mweight[0], propsd, 0, 1);
    MHp->logratio -= dnorm(Mweight[0], oldwt, propsd, 1) - pnorm(lb, oldwt, propsd, 0, 1);
  }else if(lb==-HUGE_VAL){
    // Only upper bound -> a constrained normal jump
    do{ Mweight[0] = rnorm(oldwt, propsd); }while(Mweight[0]>ub);
      
    MHp->logratio += dnorm(oldwt, Mweight[0], propsd, 1) - pnorm(ub, Mweight[0], propsd, 1, 1);
    MHp->logratio -= dnorm(Mweight[0], oldwt, propsd, 1) - pnorm(ub, oldwt, propsd, 1, 1);
  }else{
    // Bounded from both sides -> uniform
    Mweight[0] = runif(lb,ub);
  }

  // Depends on the reference measure
  MHp->logratio += -(Mweight[0]*Mweight[0]-oldwt*oldwt)/2;
}

/*********************
 void MH_Unif

 Default MH algorithm for continuous-uniform-reference ERGM
*********************/
void MH_Unif(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize Unif 
    MHp->ntoggles=1;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    Mweight[0] = runif(a,b);
  }while(Mweight[0]==oldwt);

  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}


/*********************
 void MH_UnifNonObserved

 Missing data MH algorithm for continuous-uniform-reference ERGM.
*********************/
void MH_UnifNonObserved(WtMHproposal *MHp, WtNetwork *nwp)  {  
  static Edge nmissing;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize Unif 
    MHp->ntoggles=1;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    nmissing = MHp->inputs[2];
    return;
  }

  if(nmissing==0){
    *Mtail = MH_FAILED;
    *Mhead = MH_IMPOSSIBLE;
    return;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  Edge rane = 1 + unif_rand() * nmissing;

  Mtail[0]=MHp->inputs[rane+2];
  Mhead[1]=MHp->inputs[nmissing+rane+2];

  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    Mweight[0] = runif(a,b);
  }while(Mweight[0]==oldwt);

  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}


/*********************
 void MH_DiscUnif

 Default MH algorithm for discrete-uniform-reference ERGM
*********************/
void MH_DiscUnif(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize DiscUnif 
    MHp->ntoggles=1;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    Mweight[0] = floor(runif(a,b+1));
  }while(Mweight[0]==oldwt);

  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}


/*********************
 void MH_DiscUnifNonObserved

 Missing data MH algorithm for discrete-uniform-reference ERGM.
*********************/
void MH_DiscUnifNonObserved(WtMHproposal *MHp, WtNetwork *nwp)  {  
  static Edge nmissing;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize DiscUnif 
    MHp->ntoggles=1;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    nmissing = MHp->inputs[2];
    return;
  }

  if(nmissing==0){
    *Mtail = MH_FAILED;
    *Mhead = MH_IMPOSSIBLE;
    return;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  Edge rane = 1 + unif_rand() * nmissing;

  Mtail[0]=MHp->inputs[rane+2];
  Mhead[1]=MHp->inputs[nmissing+rane+2];

  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    Mweight[0] = floor(runif(a,b+1));
  }while(Mweight[0]==oldwt);

  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}
