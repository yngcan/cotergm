#include "MHproposals.h"
#include "edgelist.h"
#include "edgetree.h"
/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mnode (MHp->togglenode)


void MH_randomtoggleList_nodal (MHproposal *MHp, Network *nwp)
{
	static Edge nedges0; static Vertex nnodes0;
	static Vertex *togglenodes;
	static double pdyad = 0.50;
	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=1;
		nedges0 = MHp->inputs[0];
		nnodes0 = MHp->inputs[nedges0*2+1];
		togglenodes = (Vertex *)malloc( nnodes0 * sizeof(Vertex));

		/*Specify which nodes are allowed to toggle */
		for (int i=0;i<nnodes0;i++)
			togglenodes[i] = (int) MHp->inputs[nedges0*2+2+i];
		return;
	}

	if(nedges0==0){ /* Attempting dissolution on a complete graph. */
		Mtail[0]=MH_FAILED;
		Mhead[0]=MH_IMPOSSIBLE;
		return;
	}

	BD_LOOP({
		/* Select a dyad at random that is in the reference graph. (We
	 have a convenient sampling frame.) */
		/* Generate. */
		double  u = unif_rand();
		if(u < pdyad){ /* toggle nodal value */

		for (int j=0;j<100;j++){
	    Vertex index_node = (int) (unif_rand() * nnodes0);
		Mnode[0] = togglenodes[index_node];
		int nplus = 0;
		for(int i=0;i<nnodes0;i++){
			if(i != index_node)
			nplus += (nwp->nodalstatus[togglenodes[i]-1] == 2)? 1:0;
		}
		if(!((nwp->nodalstatus[togglenodes[index_node]-1] ==1.0 & nplus == (nnodes0-1)) |(nwp->nodalstatus[togglenodes[index_node]-1] ==2.0 & nplus == 1))) break;
		}

		Mtail[0]=MH_NODE;
		Mhead[0]=MH_NODE;
		} else{
		Edge rane = 1 + unif_rand() * nedges0;
		Mnode[0]=0;
		Mtail[0]=MHp->inputs[rane];
		Mhead[0]=MHp->inputs[nedges0+rane];}
	});

}




/*********************
 void MH_randomtoggle

 Default MH algorithm
 *********************/
void MH_randomtoggle (MHproposal *MHp, Network *nwp)  {  

	/* *** don't forget tail-> head now */

	if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
		MHp->ntoggles=1;
		return;
	}

	BD_LOOP({
		GetRandDyad(Mtail, Mhead, nwp);
	});
}

/********************
   void MH_TNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
 ***********************/
void MH_TNT (MHproposal *MHp, Network *nwp) 
{
	/* *** don't forget tail-> head now */

	Edge nedges=nwp->nedges;
	static double comp=0.5;
	static double odds;
	static Dyad ndyads;

	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=1;
		odds = comp/(1.0-comp);
		ndyads = DYADCOUNT(nwp->nnodes, nwp->bipartite, nwp->directed_flag);
		return;
	}

	double logratio=0;
	BD_LOOP({
		if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
			GetRandEdge(Mtail, Mhead, nwp);
			/* Thanks to Robert Goudie for pointing out an error in the previous
	   version of this sampler when proposing to go from nedges==0 to nedges==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
			logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
					nedges / (odds*ndyads + nedges)));
		}else{ /* Select a dyad at random */
			GetRandDyad(Mtail, Mhead, nwp);

			if(EdgetreeSearch(Mtail[0],Mhead[0],nwp->outedges)!=0){
				logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
						nedges / (odds*ndyads + nedges)));
			}else{
				logratio = log((nedges==0 ? comp*ndyads + (1.0-comp) :
						1.0 + (odds*ndyads)/(nedges + 1)));
			}
		}
	});
	MHp->logratio += logratio;
}

/********************
   void MH_TriNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
 ***********************/
void MH_TriNT (MHproposal *MHp, Network *nwp) 
{
	/* *** don't forget tail-> head now */

	Edge nedges=nwp->nedges;
	const double comp=0.5, pthird=0.5;
	static double odds;
	static Dyad ndyads;

	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=3;
		odds = comp/(1.0-comp);
		ndyads = DYADCOUNT(nwp->nnodes, nwp->bipartite, nwp->directed_flag);
		return;
	}

	double logratio=0;
	BD_LOOP({
		if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
			GetRandEdge(Mtail, Mhead, nwp);
			/* Thanks to Robert Goudie for pointing out an error in the previous
	   version of this sampler when proposing to go from nedges==0 to nedges==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
			logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
					nedges / (odds*ndyads + nedges)));
		}else{ /* Select a dyad at random */
			GetRandDyad(Mtail, Mhead, nwp);

			if(EdgetreeSearch(Mtail[0],Mhead[0],nwp->outedges)!=0){
				logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
						nedges / (odds*ndyads + nedges)));
			}else{
				logratio = log((nedges==0 ? comp*ndyads + (1.0-comp) :
						1.0 + (odds*ndyads)/(nedges + 1)));
			}
		}

		Vertex third;
		do{
			third = 1 + unif_rand()*nwp->nnodes;
		}while(third==*Mtail || third==*Mhead);

		MHp->ntoggles=1;
		if(unif_rand()<pthird){
			*(Mtail+MHp->ntoggles) = nwp->directed_flag ? *Mtail : MIN(*Mtail, third);
			*(Mhead+MHp->ntoggles) = nwp->directed_flag ? third : MAX(*Mtail, third);
			MHp->ntoggles++;
		}
		if(unif_rand()<pthird){
			*(Mtail+MHp->ntoggles) = nwp->directed_flag ? third : MIN(*Mhead, third);
			*(Mhead+MHp->ntoggles) = nwp->directed_flag ? *Mhead : MAX(*Mhead, third);
			MHp->ntoggles++;
		}
	});

	MHp->logratio += logratio;
}


/********************
   void MH_TNT10
   Attempts to do 10 TNT steps at once, but this seems flawed currently
   because it does not correctly update network quantities like nedges
   after each of the 10 proposed toggles.
 ***********************/
void MH_TNT10 (MHproposal *MHp, Network *nwp) 
{
	/* *** don't forget tail-> head now */

	Edge nedges=nwp->nedges;
	static double comp=0.5;
	static double odds;
	static Dyad ndyads;

	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=10;
		odds = comp/(1.0-comp);
		ndyads = DYADCOUNT(nwp->nnodes, nwp->bipartite, nwp->directed_flag);
		return;
	}

	double logratio;
	BD_LOOP({
		logratio = 0;
		for(unsigned int n = 0; n < 10; n++){
			if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
				GetRandEdge(Mtail, Mhead, nwp);
				logratio += log(nedges  / (odds*ndyads + nedges));
			}else{ /* Select a dyad at random */
				GetRandDyad(Mtail+n, Mhead+n, nwp);
				if(EdgetreeSearch(Mtail[n],Mhead[n],nwp->outedges)!=0){
					logratio += log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
							nedges / (odds*ndyads + nedges)));
				}else{
					logratio += log((nedges==0 ? comp*ndyads + (1.0-comp) :
							1.0 + (odds*ndyads)/(nedges + 1)));
				}
			}
		}
	});
	MHp->logratio += logratio;
}

/*********************
 void MH_constantedges
 propose pairs of toggles that keep number of edges
 the same.  This is done by (a) choosing an existing edge
 at random; (b) repeatedly choosing dyads at random until 
 one is found that does not have an edge; and (c) proposing
 toggling both these dyads.  Note that step (b) will be very 
 inefficient if the network is nearly complete, so this proposal is
 NOT recommended for such networks.  However, most network
 datasets are sparse, so this is not likely to be an issue.
 *********************/
void MH_ConstantEdges (MHproposal *MHp, Network *nwp)  {  
	/* *** don't forget tail-> head now */

	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=2;
		return;
	}
	/* Note:  This proposal cannot be used for full or empty observed graphs.
     If desired, we could check for this at initialization phase. 
     (For now, however, no way to easily return an error message and stop.)*/
	BD_LOOP({
		/* First, select edge at random */
		GetRandEdge(Mtail, Mhead, nwp);
		/* Second, select dyad at random until it has no edge */
		do{
			GetRandDyad(Mtail+1, Mhead+1, nwp);
		}while(EdgetreeSearch(Mtail[1], Mhead[1], nwp->outedges) != 0);
	});
}

/*********************
 void MH_CondDegreeDist
 It used to be called  MH_CondDegDistSwapToggles
 *********************/
void MH_CondDegreeDist (MHproposal *MHp, Network *nwp) {  
	int noutedge=0, ninedge=0, k, fvalid;
	int k0, j0, j1, k1;
	int j0h, j1h;
	int trynode;
	Vertex e, alter, tail=0, head, head1;

	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=2;
		return;
	}

	fvalid = 0;
	trynode = 0;
	while(fvalid==0 && trynode < 500){

		trynode++;
		/* select a node at random */
		while(noutedge+ninedge==0){
			/* select a node at random */
			tail = 1 + unif_rand() * nwp->nnodes;
			ninedge  = nwp->indegree[tail];
			noutedge = nwp->outdegree[tail];
		}

		/* choose a edge of the node at random */
		/* *** don't forget tail-> head now */

		k0 = (int)(unif_rand() * (noutedge+ninedge));
		if (k0 < noutedge){
			k=0;
			for(e = EdgetreeMinimum(nwp->outedges, tail);
					((head = nwp->outedges[e].value) != 0 && k<k0);
					e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
		}else{
			k=0;
			for(e = EdgetreeMinimum(nwp->inedges, tail);
					((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
					e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
		}

		if ( (!nwp->directed_flag && tail > head) ||
				(nwp->directed_flag && k0 >= noutedge) ) {
			Mtail[0] = head;
			Mhead[0] = tail;
		}else{
			Mtail[0] = tail;
			Mhead[0] = head;
		}

		k1=0;
		fvalid=0;
		while(fvalid==0 && k1 < 100){
			while((alter = 1 + unif_rand() * nwp->nnodes) == tail);
			fvalid=1;
			if(alter == head){fvalid=0;}
			if (k0 < noutedge || !nwp->directed_flag){
				for(e = EdgetreeMinimum(nwp->outedges, tail);
						(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
						e = EdgetreeSuccessor(nwp->outedges, e)){
					if(alter==head1){fvalid=0;}}
			}
			if (k0 >= noutedge || !nwp->directed_flag){
				for(e = EdgetreeMinimum(nwp->inedges, tail);
						(fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
						e = EdgetreeSuccessor(nwp->inedges, e)){
					if(alter==head1){fvalid=0;}}
			}
			k1++;
		}

		if (k1 == 100){
			fvalid=0;
			continue;
		}

		if ( (!nwp->directed_flag && alter > tail) ||
				(nwp->directed_flag && k0 < noutedge) )
		{
			Mtail[1] = tail;
			Mhead[1] = alter;
		}else{
			Mtail[1] = alter;
			Mhead[1] = tail;
		}

		if (!nwp->directed_flag){
			/* Check undirected degrees */
			k0 =nwp->outdegree[tail]  + nwp->indegree[tail];
			j0h=nwp->outdegree[head]  + nwp->indegree[head];
			j1h=nwp->outdegree[alter] + nwp->indegree[alter];

			j0=j0h-1;
			j1=j1h+1;

			if( ( (j0==j1h) && (j1==j0h) ) ){
				fvalid = 1;
			}else{
				fvalid = 0;
			}
		}else{
			/* Check directed degrees */
			if(k0 < noutedge){
				/* Check indegrees */
				j0h=nwp->indegree[head];
				j1h=nwp->indegree[alter];
			}else{
				/* Check outdegrees */
				j0h=nwp->outdegree[head];
				j1h=nwp->outdegree[alter];
			}
			j0=j0h-1;
			j1=j1h+1;

			if( ( (j0==j1h) && (j1==j0h) ) ){
				fvalid = 1;
			}else{
				fvalid = 0;
			}
		}

	}

	if (trynode==500){
		Mtail[1] = Mtail[0];
		Mhead[1] = Mhead[0];
	}
}

/*********************
 void MH_CondOutDegreeDist
 *********************/
void MH_CondOutDegreeDist (MHproposal *MHp, Network *nwp) {  
	int noutedge=0, k, fvalid=0;
	int k0, k1;
	int trynode;
	Vertex e, alter, tail=0, head, head1;

	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=2;
		return;
	}

	fvalid = 0;
	trynode = 0;
	while(fvalid==0 && trynode < 1500){

		trynode++;

		while(noutedge==0){
			/* select a node at random */
			tail = 1 + unif_rand() * nwp->nnodes;
			noutedge = nwp->outdegree[tail];
		}

		k0 = (int)(unif_rand() * noutedge);
		k=0;
		for(e = EdgetreeMinimum(nwp->outedges, tail);
				((head = nwp->outedges[e].value) != 0 && k<k0);
				e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
		Mtail[0] = tail;
		Mhead[0] = head;

		k1=0;
		fvalid=0;
		while(fvalid==0 && k1 < 100){
			while((alter = 1 + unif_rand() * nwp->nnodes) == tail);
			fvalid=1;
			if(alter == head){fvalid=0;}
			for(e = EdgetreeMinimum(nwp->outedges, tail);
					(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
					e = EdgetreeSuccessor(nwp->outedges, e)){
				if(alter==head1){fvalid=0;}}
			k1++;
		}
		if (k1 == 100){
			fvalid=0;
			continue;
		}

		Mtail[1] = tail;
		Mhead[1] = alter;
	}

	if(trynode==1500 || !CheckTogglesValid(MHp, nwp)){
		Mtail[0] = 1;
		Mhead[0] = 2;
		Mtail[1] = 1;
		Mhead[1] = 2;
	}


}

/*********************
 void MH_CondInDegreeDist
 *********************/
void MH_CondInDegreeDist (MHproposal *MHp, Network *nwp) {  
	int ninedge=0, k, fvalid=0;
	int k0, k1;
	int trynode;
	Vertex e, alter, tail=0, head, head1;

	/* *** don't forget tail-> head now */


	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=2;
		return;
	}

	fvalid = 0;
	trynode = 0;
	while(fvalid==0 && trynode < 1500){

		trynode++;

		while(ninedge==0){
			/* select a node at random */
			tail = 1 + unif_rand() * nwp->nnodes;
			ninedge = nwp->indegree[tail];
		}

		k0 = (int)(unif_rand() * ninedge);
		k=0;
		for(e = EdgetreeMinimum(nwp->inedges, tail);
				((head = nwp->inedges[e].value) != 0 && k<k0);
				e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
		Mtail[0] = head;
		Mhead[0] = tail;

		k1=0;
		fvalid=0;
		while(fvalid==0 && k1 < 100){
			while((alter = 1 + unif_rand() * nwp->nnodes) == tail);
			fvalid=1;
			if(alter == head){fvalid=0;}
			for(e = EdgetreeMinimum(nwp->inedges, tail);
					(fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
					e = EdgetreeSuccessor(nwp->inedges, e)){
				if(alter==head1){fvalid=0;}}
			k1++;
		}
		if (k1 == 100){
			fvalid=0;
			continue;
		}

		Mtail[1] = alter;
		Mhead[1] = tail;

	}

	if(trynode==1500){
		Mtail[0] = 1;
		Mhead[0] = 2;
		Mtail[1] = 1;
		Mhead[1] = 2;
	}
}

/*********************
 void MH_TwoRandomToggles
 *********************/
void MH_TwoRandomToggles (MHproposal *MHp, Network *nwp) {  
	Vertex tail, head;
	int i;

	/* *** don't forget tail-> head now */

	if(MHp->ntoggles == 0) { /* Initialize OneRandomToggle */
		MHp->ntoggles=2;
		return;
	}

	for (i = 0; i < 2; i++){
		tail = 1 + unif_rand() * nwp->nnodes;
		while ((head = 1 + unif_rand() * nwp->nnodes) == tail);
		if (!nwp->directed_flag && tail > head) {
			Mtail[i] = head;
			Mhead[i] = tail;
		}else{
			Mtail[i] = tail;
			Mhead[i] = head;
		}
	}
}

/*********************
 void MH_RandomNode
 *********************/
void MH_randomnode (MHproposal *MHp, Network *nwp) {

	Vertex root, alter;
	int j;

	if(MHp->ntoggles == 0) { /* Initialize OneRandomToggle */
		MHp->ntoggles= nwp->nnodes - 1;
		return;
	}

	root = 1 + unif_rand() * nwp->nnodes;

	j = 0;
	for (alter = 1; alter <= nwp->nnodes; alter++)
	{
		/* there is never an edge (root, root) */
		if (alter != root) {
			if (!nwp->directed_flag && root > alter) {
				Mtail[j] = alter;
				Mhead[j] = root;
			}else{
				Mtail[j] = root;
				Mhead[j] = alter;
			}
			j++;
		}
	}
}

/********************
   void MH_randomtoggleList
   Propose ONLY edges on a static list
 ***********************/
void MH_randomtoggleList (MHproposal *MHp, Network *nwp) 
{  
	static Edge nedges0;

	if(MHp->ntoggles == 0) { /* Initialize */
		MHp->ntoggles=1;
		nedges0 = MHp->inputs[0];
		return;
	}

	if(nedges0==0){ /* Attempting dissolution on a complete graph. */
		Mtail[0]=MH_FAILED;
		Mhead[0]=MH_IMPOSSIBLE;
		return;
	}

	BD_LOOP({
		/* Select a dyad at random that is in the reference graph. (We
	 have a convenient sampling frame.) */
		/* Generate. */
		Edge rane = 1 + unif_rand() * nedges0;
		Mtail[0]=MHp->inputs[rane];
		Mhead[0]=MHp->inputs[nedges0+rane];
	});
}



/* The ones below have not been tested */

/*********************
 void MH_ConstrainedCondOutDegDist
 *********************/
void MH_ConstrainedCondOutDegDist (MHproposal *MHp, Network *nwp){  
	int noutedge=0, k, fvalid=0;
	int k0, k1;
	Vertex e, alter, tail, head, head1;

	/* *** don't forget tail-> head now */

	while(noutedge==0){
		/* select a node at random */
		tail = 1 + unif_rand() * nwp->nnodes;
		noutedge = nwp->outdegree[tail];
	}

	k0 = (int)(unif_rand() * noutedge);
	k=0;
	for(e = EdgetreeMinimum(nwp->outedges, tail);
			((head = nwp->outedges[e].value) != 0 && k<k0);
			e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
	Mtail[0] = tail;
	Mhead[0] = head;

	k1=0;
	fvalid=0;
	while(fvalid==0 && k1 < 100){
		while((alter = 1 + unif_rand() * nwp->nnodes) == tail);
		fvalid=1;
		if(alter == head){fvalid=0;}
		for(e = EdgetreeMinimum(nwp->outedges, tail);
				(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
				e = EdgetreeSuccessor(nwp->outedges, e)){
			if(alter==head1){fvalid=0;}}
		k1++;
	}
	if (k1 == 100){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}

	Mtail[1] = tail;
	Mhead[1] = alter;

	if (!fvalid){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}

	for(k=0; k < 2; k++){
		if (dEdgeListSearch(Mtail[k], Mhead[k], MHp->inputs)==0){
			Mtail[0] = Mhead[0] = 0;
			Mtail[1] = Mhead[1] = 0;
		}
	}
}


void MH_NodePairedTiesToggles (MHproposal *MHp, Network *nwp) {  
	/* chooses a node and toggles all ties and
	 and toggles an equal number of matching nonties
	 for that node */
	int nedge=0,j,k;
	int fvalid = 1;
	Vertex e, tail, prop;

	/* *** don't forget tail-> head now */

	/* double to integer coercion */
	tail = 1 + unif_rand() * nwp->nnodes;

	for(e = EdgetreeMinimum(nwp->outedges, tail);
			(prop = nwp->outedges[e].value) != 0; /* loop if */
			e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
	{
		Mtail[nedge] = tail;
		Mhead[nedge] = prop;
		++nedge;
	}
	for(e = EdgetreeMinimum(nwp->inedges, tail);
			(prop = nwp->inedges[e].value) != 0; /* loop if */
			e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
	{
		Mhead[nedge] = tail;
		Mtail[nedge] = prop;
		++nedge;
	}

	if(nedge > nwp->nnodes-nedge){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}
	j = 0;
	while (j <=nedge)
	{
		prop = 1 + unif_rand() * nwp->nnodes;
		k=0;
		fvalid=1;
		while(fvalid==1 && k<nedge+j){
			if(EdgetreeSearch( MIN(prop,Mtail[k]),
					MAX(prop,Mtail[k]), nwp->outedges) +
					EdgetreeSearch( MIN(prop,Mhead[k]),
							MAX(prop,Mhead[k]), nwp->outedges)==0
			){++k;
			}else{
				fvalid=0;
			}
		}
		if(prop>tail){
			Mtail[j+nedge] = tail;
			Mhead[j+nedge] = prop;
		}else{
			Mtail[j+nedge] = prop;
			Mhead[j+nedge] = tail;
		}
		++j;
	}

	j = 2*nedge;
	if (!CheckTogglesValid(MHp, nwp))
	{
		*Mtail = *Mhead = 0;
	}
}

/*********************
 void MH_OneRandomTnTNode
 *********************/
void MH_OneRandomTnTNode (MHproposal *MHp, Network *nwp) {  
	Vertex tail=0, head, e, head1;
	int noutedge=0, ninedge=0, k0=0, ndyad, fvalid=0, k;

	/* *** don't forget tail-> head now */

	if ( nwp->directed_flag )
	{
		ndyad = (nwp->nnodes - 1) * nwp->nnodes;
	}else{
		ndyad = (nwp->nnodes - 1) * nwp->nnodes / 2;
	}

	double logratio=0;
	fvalid=0;
	while(fvalid==0){

		if ( unif_rand() < 0.5 && nwp->nedges > 0)
		{

			/* select a tie */
			ninedge=0;
			noutedge=0;
			while(noutedge+ninedge==0){
				/* select a node at random */
				tail = 1 + unif_rand() * nwp->nnodes;
				ninedge = nwp->indegree[tail];
				noutedge = nwp->outdegree[tail];
			}

			k0 = (int)(unif_rand() * (noutedge+ninedge));
			if (k0 < noutedge){
				k=0;
				for(e = EdgetreeMinimum(nwp->outedges, tail);
						((head = nwp->outedges[e].value) != 0 && k<k0);
						e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
			}else{
				k=0;
				for(e = EdgetreeMinimum(nwp->inedges, tail);
						((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
						e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
			}
			if ( (!nwp->directed_flag && tail > head) ||
					(nwp->directed_flag && k0 >= noutedge) )
			{
				Mtail[0] = head;
				Mhead[0] = tail;
			}else{
				Mtail[0] = tail;
				Mhead[0] = head;
			}

			logratio = log(((noutedge+ninedge)*1.0)/(nwp->nnodes-1-noutedge-ninedge-1));
			fvalid =1;
		}else{
			/* Choose random non-tie */

			/* select a node at random */
			ninedge=nwp->nnodes-1;
			noutedge=0;
			while(noutedge+ninedge>=(nwp->nnodes-1)){
				ninedge=0;
				/* select a node at random */
				tail = 1 + unif_rand() * nwp->nnodes;
				ninedge = nwp->indegree[tail];
				noutedge = nwp->outdegree[tail];
			}

			fvalid=0;
			while(fvalid==0){
				while ((head = 1 + unif_rand() * nwp->nnodes) == tail);
				fvalid=1;
				for(e = EdgetreeMinimum(nwp->outedges, tail);
						(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
						e = EdgetreeSuccessor(nwp->outedges, e)){
					if(head==head1){fvalid=0;}}
				if (!(nwp->directed_flag)){
					for(e = EdgetreeMinimum(nwp->inedges, tail);
							(fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
							e = EdgetreeSuccessor(nwp->inedges, e)){
						if(head==head1){fvalid=0;}}
				}
			}

			if ( (!nwp->directed_flag && tail > head) ||
					(nwp->directed_flag && k0 >= noutedge) )
			{
				Mtail[0] = head;
				Mhead[0] = tail;
			}else{
				Mtail[0] = tail;
				Mhead[0] = head;
			}

			if ( nwp->directed_flag )
			{
				logratio = log((nwp->nnodes-1-noutedge-ninedge)/(noutedge+ninedge+1.0));
			}else{
				logratio = log((nwp->nnodes-1-noutedge-ninedge)/(noutedge+ninedge+1.0));
			}
		}
	}
	MHp->logratio += logratio;
}

/*********************
 void MH_ReallocateWithReplacement
 *********************/
void MH_ReallocateWithReplacement (MHproposal *MHp, Network *nwp) {  
	int i;
	Vertex root;
	Vertex* edges;
	int edgecount = 0;

	/* select a node at random */
	root = 1 + unif_rand() * nwp->nnodes;

	edges = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
	for (i = 0; i <= nwp->nnodes; i++)
		edges[i] = NO_EDGE;

	/* count current edges and mark them in an array */
	for (i = 1; i <= nwp->nnodes; i++)
	{
		if (root == i) continue;
		if (EdgetreeSearch(root, i, nwp->outedges) > 0)
		{
			edges[i] = OLD_EDGE;
			edgecount++;
		}
		if (!nwp->directed_flag && (root > i) &&
				(EdgetreeSearch(i, root, nwp->outedges) > 0))
		{
			edges[i] = OLD_EDGE;
			edgecount++;
		}
	}

	/* select edgecount edges to create */
	for (i = 0; i < edgecount; i++)
	{
		Vertex newhead;
		/* get a new edge, neither the root nor something already chosen */
		while ((newhead = 1 + unif_rand() * nwp->nnodes) == root ||
				(edges[newhead] & NEW_EDGE))
			;

		/* if this edge already exists - (OLD_EDGE | NEW_EDGE) == CAN_IGNORE */
		edges[newhead] = edges[newhead] | NEW_EDGE;
	}

	/* index into Mtail/Mhead is  */
	edgecount = 0;

	/* add to toggle list:  anything that is non zero in edges array
     should be toggled, whether on or off. */
	for (i = 0; i <= nwp->nnodes; i++)
	{
		if (edges[i] == NO_EDGE || edges[i] == CAN_IGNORE) continue;

		/* double to integer coercion */
		Mtail[edgecount] = root;
		Mhead[edgecount] = i;

		if (!nwp->directed_flag && (Mtail[edgecount] > Mhead[edgecount]))
		{
			Vertex temp;
			temp = Mtail[edgecount];
			Mtail[edgecount] = Mhead[edgecount];
			Mhead[edgecount] = temp;
		}
		edgecount++;
	}
	free(edges);
}

/*********************
 void MH_AllTogglesForOneNode
 *********************/
void MH_AllTogglesForOneNode (MHproposal *MHp, Network *nwp) {

	int i;
	int j;
	int root;

	root = 1 + unif_rand() * nwp->nnodes;

	j = 0;
	for (i = 1; i <= nwp->nnodes; i++)
	{
		/* probability here only do this with .8? */

		/* there is never an edge (root, root) */
		if (i == root)
			continue;

		/* double to integer coercion */
		Mtail[j] = root;
		Mhead[j] = i;

		if (!nwp->directed_flag && (Mtail[j] > Mhead[j]))
		{
			Vertex temp;
			temp = Mtail[j];
			Mtail[j] = Mhead[j];
			Mhead[j] = temp;
		}
		j++;
	}
}


/*********************
 void MH_SwitchLabelTwoNodesToggles
 *********************/
void MH_SwitchLabelTwoNodesToggles (MHproposal *MHp, Network *nwp) {  
	int nedge1=0, nedge2=0, k, ntoggles;
	Vertex *edges1, *edges2;
	Vertex e, tail2, head2, tail1, head1;

	/* *** don't forget tail-> head now */

	/* select a node at random */
	edges1 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
	edges2 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));

	while(nedge1==0){
		tail1 = 1 + unif_rand() * nwp->nnodes;

		for(e = EdgetreeMinimum(nwp->outedges, tail1);
				(head1 = nwp->outedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
		{
			edges1[nedge1] = head1;
			++nedge1;
		}
		for(e = EdgetreeMinimum(nwp->inedges, tail1);
				(head1 = nwp->inedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
		{
			edges1[nedge1] = head1;
			++nedge1;
		}
	}

	while((tail2 = 1 + unif_rand() * nwp->nnodes) == tail1);

	for(e = EdgetreeMinimum(nwp->outedges, tail2);
			(head2 = nwp->outedges[e].value) != 0; /* loop if */
			e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
	{
		edges2[nedge2] = head2;
		++nedge2;
	}
	for(e = EdgetreeMinimum(nwp->inedges, tail2);
			(head2 = nwp->inedges[e].value) != 0; /* loop if */
			e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
	{
		edges2[nedge2] = head2;
		++nedge2;
	}

	ntoggles = 0;
	for(k=0; k < nedge1; k++){
		if (tail1 > edges1[k])
		{
			Mtail[ntoggles] = edges1[k];
			Mhead[ntoggles] = tail1;
		}
		if (tail1 < edges1[k]){
			Mtail[ntoggles] = tail1;
			Mhead[ntoggles] = edges1[k];
		}
		if(tail1 != edges1[k]) ntoggles++;
	}

	for(k=0; k < nedge2; k++){
		if (tail1 > edges2[k])
		{
			Mtail[ntoggles] = edges2[k];
			Mhead[ntoggles] = tail1;
		}
		if (tail1 < edges2[k]){
			Mtail[ntoggles] = tail1;
			Mhead[ntoggles] = edges2[k];
		}
		if(tail1 != edges2[k]) ntoggles++;
	}

	for(k=0; k < nedge2; k++){
		if (tail2 > edges2[k])
		{
			Mtail[ntoggles] = edges2[k];
			Mhead[ntoggles] = tail2;
		}
		if (tail2 < edges2[k]){
			Mtail[ntoggles] = tail2;
			Mhead[ntoggles] = edges2[k];
		}
		if(tail2 != edges2[k]) ntoggles++;
	}

	for(k=0; k < nedge1; k++){
		if (tail2 > edges1[k])
		{
			Mtail[ntoggles] = edges1[k];
			Mhead[ntoggles] = tail2;
		}
		if (tail2 < edges1[k]){
			Mtail[ntoggles] = tail2;
			Mhead[ntoggles] = edges1[k];
		}
		if(tail2 != edges1[k]) ntoggles++;
	}
	free(edges1);
	free(edges2);
}


/*********************
 void MH_ConstrainedCondDegDist
 *********************/
void MH_ConstrainedCondDegDist (MHproposal *MHp, Network *nwp)  {  
	int noutedge=0, ninedge=0, k, fvalid=0;
	int k0, j0, j1, k1;
	int j0h, j1h;
	Vertex *outedges, *inedges;
	Vertex e, alter, tail=0, head;

	/* *** don't forget tail-> head now */

	/* select a node at random */
	outedges = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
	inedges = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));

	while(noutedge==0 && ninedge==0){
		tail = 1 + unif_rand() * nwp->nnodes;

		for(e = EdgetreeMinimum(nwp->outedges, tail);
				(head = nwp->outedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
		{
			outedges[noutedge] = head;
			++noutedge;
		}
		for(e = EdgetreeMinimum(nwp->inedges, tail);
				(head = nwp->inedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
		{
			inedges[ninedge] = head;
			++ninedge;
		}
	}

	k0 = (int)(unif_rand() * (noutedge+ninedge));
	if (k0 < noutedge){
		head = outedges[k0];
	}else{
		head = inedges[k0-noutedge];
	}
	if ( (!nwp->directed_flag && tail > head) ||
			(  nwp->directed_flag  && k0 >= noutedge) )
	{
		Mtail[0] = head;
		Mhead[0] = tail;
	}else{
		Mtail[0] = tail;
		Mhead[0] = head;
	}

	if (dEdgeListSearch(Mtail[0], Mhead[0], MHp->inputs)==0){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}

	fvalid=0;
	k1=0;
	while(fvalid==0 && k1 < 100){
		while((alter = 1 + unif_rand() * nwp->nnodes) == tail);
		if(alter != head){fvalid=1;}
		fvalid=1;
		if (k0 < noutedge || !(nwp->directed_flag)){
			k=0;
			while(fvalid==1 && noutedge > 0 && k <= noutedge-1){
				if(alter == outedges[k]){fvalid=0;}else{++k;}
			}
		}
		if (k0 >= noutedge || !(nwp->directed_flag)){
			k=0;
			while(fvalid==1 && ninedge > 0 && k <= ninedge-1){
				if(alter == inedges[k]){fvalid=0;}else{++k;}
			}
		}
		k1++;
	}

	if (k1 == 100){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}

	if ( (!nwp->directed_flag && alter > tail) ||
			(nwp->directed_flag && k0 < noutedge) )
	{
		Mtail[1] = tail;
		Mhead[1] = alter;
	}else{
		Mtail[1] = alter;
		Mhead[1] = tail;
	}

	if (dEdgeListSearch(Mtail[1], Mhead[1], MHp->inputs)==0){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}

	free(outedges);
	free(inedges);

	/* Check undirected degrees */

	/* *** don't forget tail-> head now */

	if (!nwp->directed_flag){
		k0=nwp->outdegree[tail]+ nwp->indegree[tail];
		j0h=nwp->outdegree[head]+ nwp->indegree[head];
		j1h=nwp->outdegree[alter]+ nwp->indegree[alter];

		j0=j0h-1;
		j1=j1h+1;

		if( ( (j0==j1h) && (j1==j0h) ) ){
			fvalid = 1;
		}else{
			fvalid = 0;
		}
	}else{
		if(k0 < noutedge){
			/* Check indegrees */
			j0h=nwp->indegree[head];
			j1h=nwp->indegree[alter];
		}else{
			/* Check outdegrees */
			j0h=nwp->outdegree[head];
			j1h=nwp->outdegree[alter];
		}
		j0=j0h-1;
		j1=j1h+1;

		if( ( (j0==j1h) && (j1==j0h) ) ){
			fvalid = 1;
		}else{
			fvalid = 0;
		}
	}

	if (!fvalid){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}
}

void MH_ConstrainedNodePairedTiesToggles (MHproposal *MHp,
		Network *nwp) {
	/* chooses a node and toggles all ties and
     and toggles an equal number of matching nonties
     for that node */
	int nedge=0,j,k;
	int fvalid = 1;
	Vertex e, tail, prop;

	/* *** don't forget tail-> head now */

	/* double to integer coercion */
	tail = 1 + unif_rand() * nwp->nnodes;

	for(e = EdgetreeMinimum(nwp->outedges, tail);
			(prop = nwp->outedges[e].value) != 0; /* loop if */
			e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
	{
		Mtail[nedge] = tail;
		Mhead[nedge] = prop;
		++nedge;
	}
	for(e = EdgetreeMinimum(nwp->inedges, tail);
			(prop = nwp->inedges[e].value) != 0; /* loop if */
			e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
	{
		Mhead[nedge] = tail;
		Mtail[nedge] = prop;
		++nedge;
	}

	if(nedge > nwp->nnodes-nedge){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}
	j = 0;
	while (j <=nedge)
	{
		prop = 1 + unif_rand() * nwp->nnodes;
		k=0;
		fvalid=1;
		while(fvalid==1 && k<nedge+j){
			if(EdgetreeSearch( MIN(prop,Mtail[k]),
					MAX(prop,Mtail[k]), nwp->outedges) +
					EdgetreeSearch( MIN(prop,Mhead[k]),
							MAX(prop,Mhead[k]), nwp->outedges)==0
			){++k;
			}else{
				fvalid=0;}
		}
		if(prop>tail){
			Mtail[j+nedge] = tail;
			Mhead[j+nedge] = prop;
		}else{
			Mtail[j+nedge] = prop;
			Mhead[j+nedge] = tail;
		}
		++j;
	}

	j = 2*nedge;
	if (!CheckConstrainedTogglesValid(MHp, nwp))
	{
		*Mtail = *Mhead = 0;
	}
}

/*********************
 void MH_ConstrainedReallocateWithReplacement
 *********************/
void MH_ConstrainedReallocateWithReplacement (MHproposal *MHp,
		Network *nwp) {
	int i;
	Vertex root;
	Vertex* edges;
	int edgecount = 0;

	/* select a node at random */
	root = 1 + unif_rand() * nwp->nnodes;

	edges = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
	for (i = 0; i <= nwp->nnodes; i++)
		edges[i] = NO_EDGE;

	/* count current edges and mark them in an array */
	for (i = 1; i <= nwp->nnodes; i++)
	{
		if (root == i) continue;
		if (EdgetreeSearch(root, i, nwp->outedges) > 0)
		{
			edges[i] = OLD_EDGE;
			edgecount++;
		}
		if (!nwp->directed_flag && (root > i) &&
				(EdgetreeSearch(i, root, nwp->outedges) > 0))
		{
			edges[i] = OLD_EDGE;
			edgecount++;
		}
	}

	/* select edgecount edges to create */
	for (i = 0; i < edgecount; i++)
	{
		Vertex newhead;

		/* get a new edge, neither the root nor something already chosen */
		while ((newhead = 1 + unif_rand() * nwp->nnodes) == root ||
				(edges[newhead] & NEW_EDGE))
			;

		/* if this edge already exists - (OLD_EDGE | NEW_EDGE) == CAN_IGNORE */
		edges[newhead] = edges[newhead] | NEW_EDGE;
	}

	/* index into Mtail/Mhead is  */
	edgecount = 0;

	/* add to toggle list:  anything that is non zero in edges array
     should be toggled, whether on or off. */
	for (i = 0; i <= nwp->nnodes; i++)
	{
		if (edges[i] == NO_EDGE || edges[i] == CAN_IGNORE) continue;

		/* double to integer coercion */
		Mtail[edgecount] = root;
		Mhead[edgecount] = i;

		if (!nwp->directed_flag && (Mtail[edgecount] > Mhead[edgecount]))
		{
			Vertex temp;
			temp = Mtail[edgecount];
			Mtail[edgecount] = Mhead[edgecount];
			Mhead[edgecount] = temp;
		}
		edgecount++;
	}
	free(edges);
}

/*********************
 void MH_ConstrainedAllTogglesForOneNode
 *********************/
void MH_ConstrainedAllTogglesForOneNode (MHproposal *MHp,
		Network *nwp) {
	int i;
	int j;
	int root;

	root = 1 + unif_rand() * nwp->nnodes;

	j = 0;
	for (i = 1; i <= nwp->nnodes; i++)
	{
		/* probability here only do this with .8? */

		/* there is never an edge (root, root) */
		if (i == root)
			continue;

		/* double to integer coercion */
		Mtail[j] = root;
		Mhead[j] = i;

		if (!nwp->directed_flag && (Mtail[j] > Mhead[j]))
		{
			Vertex temp;
			temp = Mtail[j];
			Mtail[j] = Mhead[j];
			Mhead[j] = temp;
		}
		j++;
	}
}

/*********************
 void MH_ConstrainedTwoRandomToggles
 *********************/
void MH_ConstrainedTwoRandomToggles (MHproposal *MHp,
		Network *nwp) {
	int i;

	for (i = 0; i < 2; i++)
	{
		/* double to integer coercion */
		Mtail[i] = 1 + unif_rand() * nwp->nnodes;
		while ((Mhead[i] = 1 + unif_rand() * nwp->nnodes) == Mtail[i]);

		while(dEdgeListSearch(Mtail[i], Mhead[i], MHp->inputs)==0){
			Mtail[i] = 1 + unif_rand() * nwp->nnodes;
			while ((Mhead[i] = 1 + unif_rand() * nwp->nnodes) == Mtail[i]);
		}
		if (!nwp->directed_flag && Mtail[i] > Mhead[i])
		{
			Vertex temp;
			temp = Mtail[i];
			Mtail[i] = Mhead[i];
			Mhead[i] = temp;
		}
	}

	if (!CheckConstrainedTogglesValid(MHp, nwp))
	{
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}
}

/*********************
 void MH_ConstrainedCondDeg
 *********************/
void MH_ConstrainedCondDeg (MHproposal *MHp,
		Network *nwp) {
	/* WARNING: THIS NEEDS TO BE FIXED */
	int nedge1=0, nedge2=0, k, toomany, fvalid=0;
	Vertex *edges1, *edges2;
	Vertex e, tail2=0, head2, tail1, head1;

	/* select a node at random */
	edges1 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
	edges2 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));

	while(nedge1==0){
		tail1 = 1 + unif_rand() * nwp->nnodes;

		for(e = EdgetreeMinimum(nwp->outedges, tail1);
				(head1 = nwp->outedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
		{
			edges1[nedge1] = head1;
			++nedge1;
		}
		for(e = EdgetreeMinimum(nwp->inedges, tail1);
				(head1 = nwp->inedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
		{
			edges1[nedge1] = head1;
			++nedge1;
		}
	}

	head1 = edges1[(int)(unif_rand() * nedge1)];
	if (tail1 > head1)
	{
		Mtail[0] = head1;
		Mhead[0] = tail1;
	}else{
		Mtail[0] = tail1;
		Mhead[0] = head1;
	}

	toomany = 0;
	while(nedge2==0 && toomany < 100){
		fvalid=0;
		while(fvalid==0){
			while((tail2 = 1 + unif_rand() * nwp->nnodes) == tail1);
			k=0;
			fvalid=1;
			while(fvalid==1 && k < nedge1){
				if(tail2 == edges1[k]){fvalid=0;}else{++k;}
			}
		}

		for(e = EdgetreeMinimum(nwp->outedges, tail2);
				(head2 = nwp->outedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
		{
			edges2[nedge2] = head2;
			++nedge2;
		}
		for(e = EdgetreeMinimum(nwp->inedges, tail2);
				(head2 = nwp->inedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
		{
			edges2[nedge2] = head2;
			++nedge2;
		}
		++toomany;
	}
	if (toomany==100){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}
	toomany=0;
	fvalid=0;
	while(fvalid==0 && toomany < 10){
		while((head2 = edges2[(int)(unif_rand() * nedge2)]) == tail1);
		k=0;
		fvalid=1;
		while(fvalid==1 && k < nedge1){
			if(head2 == edges1[k]){fvalid=0;}else{++k;}
		}
		++toomany;
	}
	if (!fvalid || toomany==10){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
		free(edges1);
		free(edges2);
	}
	if (tail2 > head2)
	{
		Mtail[1] = head2;
		Mhead[1] = tail2;
	}else{
		Mtail[1] = tail2;
		Mhead[1] = head2;
	}
	free(edges1);
	free(edges2);
}

/*********************
 void MH_ConstrainedSwitchLabelTwoNodesToggles
 *********************/
void MH_ConstrainedSwitchLabelTwoNodesToggles (MHproposal *MHp,
		Network *nwp)  {
	int nedge1=0, nedge2=0, k, ntoggles;
	Vertex *edges1, *edges2;
	Vertex e, tail2, head2, tail1, head1;

	/* *** don't forget tail-> head now */

	/* select a node at random */

	edges1 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
	edges2 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));

	while(nedge1==0){
		tail1 = 1 + unif_rand() * nwp->nnodes;

		for(e = EdgetreeMinimum(nwp->outedges, tail1);
				(head1 = nwp->outedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
		{
			edges1[nedge1] = head1;
			++nedge1;
		}
		for(e = EdgetreeMinimum(nwp->inedges, tail1);
				(head1 = nwp->inedges[e].value) != 0; /* loop if */
				e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
		{
			edges1[nedge1] = head1;
			++nedge1;
		}
	}

	while((tail2 = 1 + unif_rand() * nwp->nnodes) == tail1);

	for(e = EdgetreeMinimum(nwp->outedges, tail2);
			(head2 = nwp->outedges[e].value) != 0; /* loop if */
			e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
	{
		edges2[nedge2] = head2;
		++nedge2;
	}
	for(e = EdgetreeMinimum(nwp->inedges, tail2);
			(head2 = nwp->inedges[e].value) != 0; /* loop if */
			e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
	{
		edges2[nedge2] = head2;
		++nedge2;
	}

	ntoggles = 0;
	for(k=0; k < nedge1; k++){
		if (tail1 > edges1[k])
		{
			Mtail[ntoggles] = edges1[k];
			Mhead[ntoggles] = tail1;
		}
		if (tail1 < edges1[k]){
			Mtail[ntoggles] = tail1;
			Mhead[ntoggles] = edges1[k];
		}
		if(tail1 != edges1[k]) ntoggles++;
	}

	for(k=0; k < nedge2; k++){
		if (tail1 > edges2[k])
		{
			Mtail[ntoggles] = edges2[k];
			Mhead[ntoggles] = tail1;
		}
		if (tail1 < edges2[k]){
			Mtail[ntoggles] = tail1;
			Mhead[ntoggles] = edges2[k];
		}
		if(tail1 != edges2[k]) ntoggles++;
	}

	for(k=0; k < nedge2; k++){
		if (tail2 > edges2[k])
		{
			Mtail[ntoggles] = edges2[k];
			Mhead[ntoggles] = tail2;
		}
		if (tail2 < edges2[k]){
			Mtail[ntoggles] = tail2;
			Mhead[ntoggles] = edges2[k];
		}
		if(tail2 != edges2[k]) ntoggles++;
	}

	for(k=0; k < nedge1; k++){
		if (tail2 > edges1[k])
		{
			Mtail[ntoggles] = edges1[k];
			Mhead[ntoggles] = tail2;
		}
		if (tail2 < edges1[k]){
			Mtail[ntoggles] = tail2;
			Mhead[ntoggles] = edges1[k];
		}
		if(tail2 != edges1[k]) ntoggles++;
	}
	free(edges1);
	free(edges2);
}

/*********************
 void MH_ConstantEdgesToggles
 *********************/
void MH_ConstantEdgesToggles (MHproposal *MHp, Network *nwp)  {  
	int noutedge=0, ninedge=0, k, fvalid=0;
	int k0, k1;
	Vertex e, alter, tail, head, head1;

	/* *** don't forget tail-> head now */

	while(noutedge+ninedge==0){
		/* select a node at random */
		tail = 1 + unif_rand() * nwp->nnodes;
		ninedge  = nwp->indegree[tail];
		noutedge = nwp->outdegree[tail];
	}

	k0 = (int)(unif_rand() * (noutedge+ninedge));
	if (k0 < noutedge){
		k=0;
		for(e = EdgetreeMinimum(nwp->outedges, tail);
				((head = nwp->outedges[e].value) != 0 && k<k0);
				e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
	}else{
		k=0;
		for(e = EdgetreeMinimum(nwp->inedges, tail);
				((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
				e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
	}

	if ( (!nwp->directed_flag && tail > head) ||
			(nwp->directed_flag && k0 >= noutedge) )
	{
		Mtail[0] = head;
		Mhead[0] = tail;
	}else{
		Mtail[0] = tail;
		Mhead[0] = head;
	}

	k1=0;
	fvalid=0;
	while(fvalid==0 && k1 < 100){
		while((alter = 1 + unif_rand() * nwp->nnodes) == tail);
		fvalid=1;
		if(alter == head){fvalid=0;}
		if (k0 < noutedge || !(nwp->directed_flag)){
			for(e = EdgetreeMinimum(nwp->outedges, tail);
					(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
					e = EdgetreeSuccessor(nwp->outedges, e)){
				if(alter==head1){fvalid=0;}}
		}
		if (k0 >= noutedge || !(nwp->directed_flag)){
			for(e = EdgetreeMinimum(nwp->inedges, tail);
					(fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
					e = EdgetreeSuccessor(nwp->inedges, e)){
				if(alter==head1){fvalid=0;}}
		}
		k1++;
	}
	if (k1 == 100){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}

	if ( (!nwp->directed_flag && alter > tail) ||
			(nwp->directed_flag && k0 < noutedge) )
	{
		Mtail[1] = tail;
		Mhead[1] = alter;
	}else{
		Mtail[1] = alter;
		Mhead[1] = tail;
	}

	if (!fvalid){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}else{
	}
}

/*********************
 void MH_CondDegSwitchToggles
 *********************/
void MH_CondDegSwitchToggles (MHproposal *MHp, Network *nwp)  {  
	int noutedge, ninedge, i;
	int k, k0, toomany;
	Vertex e, tail, head;

	/* *** don't forget tail-> head now */

	/* select a node at random */
	for (i = 0; i < 2; i++){
		toomany=0;
		noutedge=0;
		ninedge=0;
		while(noutedge==0 && ninedge==0 && toomany < 100){
			tail = 1 + unif_rand() * nwp->nnodes;
			ninedge=0;
			noutedge=0;
			while(noutedge+ninedge==0){
				/* select a node at random */
				tail = 1 + unif_rand() * nwp->nnodes;
				ninedge = nwp->indegree[tail];
				noutedge = nwp->outdegree[tail];
			}
			++toomany;
		}

		if (toomany == 100){
			Mtail[0] = Mhead[0] = 0;
			Mtail[1] = Mhead[1] = 0;
		}

		k0 = (int)(unif_rand() * (noutedge+ninedge));
		if (k0 < noutedge){
			k=0;
			for(e = EdgetreeMinimum(nwp->outedges, tail);
					((head = nwp->outedges[e].value) != 0 && k<k0);
					e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
		}else{
			k=0;
			for(e = EdgetreeMinimum(nwp->inedges, tail);
					((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
					e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
		}
		if ( (!nwp->directed_flag && tail > head) ||
				(nwp->directed_flag && k0 >= noutedge) )
		{
			Mtail[i] = head;
			Mhead[i] = tail;
		}else{
			Mtail[i] = tail;
			Mhead[i] = head;
		}
	}

	if (EdgetreeSearch( Mtail[0],Mhead[1], nwp->outedges) ||
			EdgetreeSearch( Mtail[1],Mhead[0], nwp->outedges) ){
		Mtail[0] = Mhead[0] = 0;
		Mtail[1] = Mhead[1] = 0;
	}

	if ( (!nwp->directed_flag && Mtail[0] > Mhead[1]) )
	{
		Mtail[2] = Mhead[1];
		Mhead[2] = Mtail[0];
	}else{
		Mtail[2] = Mtail[0];
		Mhead[2] = Mhead[1];
	}

	if ( (!nwp->directed_flag && Mtail[1] > Mhead[0]) )
	{
		Mtail[3] = Mhead[0];
		Mhead[3] = Mtail[1];
	}else{
		Mtail[3] = Mtail[1];
		Mhead[3] = Mhead[0];
	}
}


