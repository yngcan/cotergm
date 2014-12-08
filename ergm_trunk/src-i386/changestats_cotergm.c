#include "changestats_cotergm.h"
#include "edgetree.h"

/*****************
 changestat: d_plus
 *****************/
D_CHANGESTAT_FN(d_plus){

	int i;
	Vertex node, e, v;

	ZERO_ALL_CHANGESTATS(i);
	FOR_EACH_TOGGLE(i) {

		node=NODE(i);

		if(node!=0){
			CHANGE_STAT[0] += NW_STATUS[node-1] == 2? -1.0 : 1.0;
		}
		TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
	}
	UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}
/*****************
 changestat: d_nodematch_cotergm
 *****************/
D_CHANGESTAT_FN(d_nodematch_cotergm) {
	double matchval, statusi, matchnomatch = 0.0,
			ismatch =0.0, isnomatch=0.0 , ismatchy0edge=0.0, isnomatchy0edge=0.0, ismatchy0noedge=0.0, isnomatchy0noedge=0.0,
			ismatchy0nodal=0.0, isnomatchy0nodal=0.0, ismatchy0nonodal=0.0, isnomatchy0nonodal=0.0, ismatchy0edgey0nodal=0.0,
			isnomatchy0edgey0nodal=0.0, ismatchy0edgey0nonodal=0.0, isnomatchy0edgey0nonodal=0.0,ismatchy0noedgey0nodal=0.0,
			isnomatchy0noedgey0nodal=0.0, ismatchy0noedgey0nonodal=0.0, isnomatchy0noedgey0nonodal=0.0;
	Vertex tail, head, node, ninputs, e, v;
	int i, j, edgeflag, edgeflag_y0,nodalflag_y0;

	ninputs = N_INPUT_PARAMS - N_NODES - 2;

	/* *** don't forget tail -> head */
	ZERO_ALL_CHANGESTATS(i);
	FOR_EACH_TOGGLE(i) {

		tail=TAIL(i);
		head=HEAD(i);
		node=NODE(i);

		if(node==0){ /*changestat for dyad change*/
			matchval = NW_STATUS[tail-1];
			if (matchval == NW_STATUS[head-1]) { /* We have a match! */
				edgeflag = IS_OUTEDGE(tail, head);

				if (ninputs==0) {/* diff=F in network statistic specification */
					CHANGE_STAT[0] += edgeflag ? -1.0 : 1.0;
				} else{
					edgeflag_y0 = IS_OUTEDGE_Y0(tail,head);
//					edgeflag_y0 = 1;
//					Rprintf("edgeflag_y0=%d \n",edgeflag_y0);
					nodalflag_y0 = Y0_STATUS[tail-1]==Y0_STATUS[head-1];
					for (j=0; j< ninputs; j++) {  //assume edgestats and nodal
						if (matchval == INPUT_PARAM[j]){
							if(INPUT_PARAM[ninputs]==0 & INPUT_PARAM[ninputs+1]==0)//
								CHANGE_STAT[j] += (edgeflag)  ? -1.0 : 1.0;

							if(INPUT_PARAM[ninputs]==0 & INPUT_PARAM[ninputs+1]==1 & nodalflag_y0)//y0nodal==1
								CHANGE_STAT[j] += (edgeflag)  ? -1.0 : 1.0;

							if(INPUT_PARAM[ninputs]==0 & INPUT_PARAM[ninputs+1]==2 & !nodalflag_y0) //y0nodal==2
								CHANGE_STAT[j] += (edgeflag)  ? -1.0 : 1.0;

							if(INPUT_PARAM[ninputs]==1 & INPUT_PARAM[ninputs+1]==0 & edgeflag_y0) //y0edge==1
								CHANGE_STAT[j] += (edgeflag)  ? -1.0 : 1.0;
					\
							if(INPUT_PARAM[ninputs]==1 & INPUT_PARAM[ninputs+1]==1 & edgeflag_y0 & nodalflag_y0) //y0edge==1,y0nodal==1
								CHANGE_STAT[j] += (edgeflag)  ? -1.0 : 1.0;

							if(INPUT_PARAM[ninputs]==1 & INPUT_PARAM[ninputs+1]==2 & edgeflag_y0 & !nodalflag_y0) //y0edge==1,y0nodal==2
								CHANGE_STAT[j] += (edgeflag)  ? -1.0 : 1.0;

							if(INPUT_PARAM[ninputs]==2 & INPUT_PARAM[ninputs+1]==0 & !edgeflag_y0 ) //y0edge==2
								CHANGE_STAT[j] += (edgeflag)  ? -1.0 : 1.0;

							if(INPUT_PARAM[ninputs]==2 & INPUT_PARAM[ninputs+1]==1 & !edgeflag_y0 & nodalflag_y0) //y0edge==2,y0nodal==1
								CHANGE_STAT[j] += (edgeflag)  ? -1.0 : 1.0;

							if(INPUT_PARAM[ninputs]==2 & INPUT_PARAM[ninputs+1]==2 & !edgeflag_y0 & !nodalflag_y0) //y0edge==2,y0nodal==2
								CHANGE_STAT[j] += (edgeflag) ? -1.0 : 1.0;
						}
					}
				}
			}
		}
		else{ /*changestat for nodal status change*/
			statusi = NW_STATUS[node-1];

			if (ninputs==0) {
				STEP_THROUGH_OUTEDGES(node,e,v){
					matchnomatch += (NW_STATUS[v-1] == statusi) ? 1.0 : -1.0;
				}

				STEP_THROUGH_INEDGES(node,e,v){
					matchnomatch += (NW_STATUS[v-1] == statusi) ? 1.0 : -1.0;
				}

				CHANGE_STAT[0] += -matchnomatch;
			} else{

				STEP_THROUGH_OUTEDGES(node,e,v){
					ismatch += (NW_STATUS[v-1] == statusi) ? 1.0 : 0.0;
					isnomatch += (NW_STATUS[v-1] == statusi) ? 0.0 : 1.0;

					edgeflag_y0 = IS_OUTEDGE_Y0(node,v);
					nodalflag_y0 = Y0_STATUS[node-1]==Y0_STATUS[v-1];

					ismatchy0edge += (NW_STATUS[v-1] == statusi & edgeflag_y0) ? 1.0 : 0.0;
					ismatchy0noedge += (NW_STATUS[v-1] == statusi & !edgeflag_y0) ? 1.0 : 0.0;
					ismatchy0nodal += (NW_STATUS[v-1] == statusi & nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0nonodal += (NW_STATUS[v-1] == statusi & !nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0edgey0nodal += (NW_STATUS[v-1] == statusi & edgeflag_y0 & nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0edgey0nonodal += (NW_STATUS[v-1] == statusi & edgeflag_y0 & !nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0noedgey0nodal += (NW_STATUS[v-1] == statusi & !edgeflag_y0 & nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0noedgey0nonodal += (NW_STATUS[v-1] == statusi & !edgeflag_y0 & !nodalflag_y0) ? 1.0 : 0.0;

					isnomatchy0edge += (NW_STATUS[v-1] != statusi & edgeflag_y0) ? 1.0 : 0.0;
					isnomatchy0noedge += (NW_STATUS[v-1] != statusi & !edgeflag_y0) ? 1.0 : 0.0;
					isnomatchy0nodal += (NW_STATUS[v-1] != statusi & nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0nonodal += (NW_STATUS[v-1] != statusi & !nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0edgey0nodal += (NW_STATUS[v-1] != statusi & edgeflag_y0 & nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0edgey0nonodal += (NW_STATUS[v-1] != statusi & edgeflag_y0 & !nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0noedgey0nodal += (NW_STATUS[v-1] != statusi & !edgeflag_y0 & nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0noedgey0nonodal += (NW_STATUS[v-1] != statusi & !edgeflag_y0 & !nodalflag_y0) ? 1.0 : 0.0;

				}
				STEP_THROUGH_INEDGES(node,e,v){
					ismatch += (NW_STATUS[v-1] == statusi) ? 1.0 : 0.0;
					isnomatch += (NW_STATUS[v-1] == statusi) ? 0.0 : 1.0;

					edgeflag_y0 = IS_INEDGE_Y0(node,v);
					nodalflag_y0 = Y0_STATUS[node-1]==Y0_STATUS[v-1];


					ismatchy0edge += (NW_STATUS[v-1] == statusi & edgeflag_y0) ? 1.0 : 0.0;
					ismatchy0noedge += (NW_STATUS[v-1] == statusi & !edgeflag_y0) ? 1.0 : 0.0;
					ismatchy0nodal += (NW_STATUS[v-1] == statusi & nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0nonodal += (NW_STATUS[v-1] == statusi & !nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0edgey0nodal += (NW_STATUS[v-1] == statusi & edgeflag_y0 & nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0edgey0nonodal += (NW_STATUS[v-1] == statusi & edgeflag_y0 & !nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0noedgey0nodal += (NW_STATUS[v-1] == statusi & !edgeflag_y0 & nodalflag_y0) ? 1.0 : 0.0;
					ismatchy0noedgey0nonodal += (NW_STATUS[v-1] == statusi & !edgeflag_y0 & !nodalflag_y0) ? 1.0 : 0.0;


					isnomatchy0edge += (NW_STATUS[v-1] != statusi & edgeflag_y0) ? 1.0 : 0.0;
					isnomatchy0noedge += (NW_STATUS[v-1] != statusi & !edgeflag_y0) ? 1.0 : 0.0;
					isnomatchy0nodal += (NW_STATUS[v-1] != statusi & nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0nonodal += (NW_STATUS[v-1] != statusi & !nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0edgey0nodal += (NW_STATUS[v-1] != statusi & edgeflag_y0 & nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0edgey0nonodal += (NW_STATUS[v-1] != statusi & edgeflag_y0 & !nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0noedgey0nodal += (NW_STATUS[v-1] != statusi & !edgeflag_y0 & nodalflag_y0) ? 1.0 : 0.0;
					isnomatchy0noedgey0nonodal += (NW_STATUS[v-1] != statusi & !edgeflag_y0 & !nodalflag_y0) ? 1.0 : 0.0;

				}

				if(ninputs==1){
					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 0 & INPUT_PARAM[2] == 0){
						CHANGE_STAT[0] += -ismatch;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 0 & INPUT_PARAM[2] == 0){
						CHANGE_STAT[0] += isnomatch;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 0 & INPUT_PARAM[2] == 1){
						CHANGE_STAT[0] += -ismatchy0nodal;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 0 & INPUT_PARAM[2] == 1){
						CHANGE_STAT[0] += isnomatchy0nodal;
					}

					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 0 & INPUT_PARAM[2] == 2){
						CHANGE_STAT[0] += -ismatchy0nonodal;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 0 & INPUT_PARAM[2] == 2){
						CHANGE_STAT[0] += isnomatchy0nonodal;
					}

					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 1 & INPUT_PARAM[2] == 0){
						CHANGE_STAT[0] += -ismatchy0edge;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 1 & INPUT_PARAM[2] == 0){
						CHANGE_STAT[0] += isnomatchy0edge;
					}

					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 1 & INPUT_PARAM[2] == 1){
						CHANGE_STAT[0] += -ismatchy0edgey0nodal;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 1 & INPUT_PARAM[2] == 1){
						CHANGE_STAT[0] += isnomatchy0edgey0nodal;
					}

					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 1 & INPUT_PARAM[2] == 2){
						CHANGE_STAT[0] += -ismatchy0edgey0nonodal;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 1 & INPUT_PARAM[2] == 2){
						CHANGE_STAT[0] += isnomatchy0edgey0nonodal;
					}

					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 2 & INPUT_PARAM[2] == 0){
						CHANGE_STAT[0] += -ismatchy0noedge;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 2 & INPUT_PARAM[2] == 0){
						CHANGE_STAT[0] += isnomatchy0noedge;
					}

					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 2 & INPUT_PARAM[2] == 1){
						CHANGE_STAT[0] += -ismatchy0noedgey0nodal;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 2 & INPUT_PARAM[2] == 1){
						CHANGE_STAT[0] += isnomatchy0noedgey0nodal;
					}

					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[1] == 2 & INPUT_PARAM[2] == 2){
						CHANGE_STAT[0] += -ismatchy0noedgey0nonodal;
					}

					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[1] == 2 & INPUT_PARAM[2] == 2){
						CHANGE_STAT[0] += isnomatchy0noedgey0nonodal;
					}

				}
				if(ninputs==2){

					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 0 & INPUT_PARAM[3] == 0){  /*statusi = minus*/
						CHANGE_STAT[0] += -ismatch;
						CHANGE_STAT[1] += isnomatch;
					}
					if(statusi == INPUT_PARAM[1] & INPUT_PARAM[2] == 0 & INPUT_PARAM[3] == 0){ /*statusi = plus*/
						CHANGE_STAT[1] += -ismatch;
						CHANGE_STAT[0] += isnomatch;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 0 & INPUT_PARAM[3] == 1){
						CHANGE_STAT[0] += -ismatchy0nodal;
						CHANGE_STAT[1] += isnomatchy0nodal;
					}
					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[2] == 0 & INPUT_PARAM[3] == 1){
						CHANGE_STAT[1] += -ismatchy0nodal;
						CHANGE_STAT[0] += isnomatchy0nodal;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 0 & INPUT_PARAM[3] == 2){
						CHANGE_STAT[0] += -ismatchy0nonodal;
						CHANGE_STAT[1] +=  isnomatchy0nonodal;
					}
					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[2] == 0 & INPUT_PARAM[3] == 2){
						CHANGE_STAT[1] += -ismatchy0nonodal;
						CHANGE_STAT[0] +=  isnomatchy0nonodal;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 1 & INPUT_PARAM[3] == 0){
						CHANGE_STAT[0] += -ismatchy0edge;
						CHANGE_STAT[1] += isnomatchy0edge;
					}
					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[2] == 1 & INPUT_PARAM[3] == 0){
						CHANGE_STAT[1] += -ismatchy0edge;
						CHANGE_STAT[0] += isnomatchy0edge;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 1 & INPUT_PARAM[3] == 1){
						CHANGE_STAT[0] += -ismatchy0edgey0nodal;
						CHANGE_STAT[1] += isnomatchy0edgey0nodal;
					}
					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[2] == 1 & INPUT_PARAM[3] == 1){
						CHANGE_STAT[1] += -ismatchy0edgey0nodal;
						CHANGE_STAT[0] += isnomatchy0edgey0nodal;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 1 & INPUT_PARAM[3] == 2){
						CHANGE_STAT[0] += -ismatchy0edgey0nonodal;
						CHANGE_STAT[1] += isnomatchy0edgey0nonodal;
					}
					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[2] == 1 & INPUT_PARAM[3] == 2){
						CHANGE_STAT[1] += -ismatchy0edgey0nonodal;
						CHANGE_STAT[0] += isnomatchy0edgey0nonodal;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 2 & INPUT_PARAM[3] == 0){
						CHANGE_STAT[0] += -ismatchy0noedge;
						CHANGE_STAT[1] += isnomatchy0noedge;
					}
					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[2] == 2 & INPUT_PARAM[3] == 0){
						CHANGE_STAT[1] += -ismatchy0noedge;
						CHANGE_STAT[0] += isnomatchy0noedge;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 2 & INPUT_PARAM[3] == 1){
						CHANGE_STAT[0] += -ismatchy0noedgey0nodal;
						CHANGE_STAT[1] += isnomatchy0noedgey0nodal;
					}
					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[2] == 2 & INPUT_PARAM[3] == 1){
						CHANGE_STAT[1] += -ismatchy0noedgey0nodal;
						CHANGE_STAT[0] += isnomatchy0noedgey0nodal;
					}


					if(statusi == INPUT_PARAM[0] & INPUT_PARAM[2] == 2 & INPUT_PARAM[3] == 2){
						CHANGE_STAT[0] += -ismatchy0noedgey0nonodal;
						CHANGE_STAT[1] += isnomatchy0noedgey0nonodal;
					}
					if(statusi != INPUT_PARAM[0] & INPUT_PARAM[2] == 2 & INPUT_PARAM[3] == 2){
						CHANGE_STAT[1] += -ismatchy0noedgey0nonodal;
						CHANGE_STAT[0] += isnomatchy0noedgey0nonodal;					}
				}
			}

		}
		TOGGLE_IF_MORE_TO_COME(i);
	}
	UNDO_PREVIOUS_TOGGLES(i);
}






/*****************
 changestat: d_gwesp
 *****************/
D_CHANGESTAT_FN(d_gwesp_cotergm) {
	Edge e, f;
	int i, echange, ochange;
	int L2th, L2tu, L2uh;
	Vertex tail, head, node, u, v;
	double alpha, oneexpa, cumchange;

	CHANGE_STAT[0] = 0.0;
	alpha = INPUT_PARAM[0];
	oneexpa = 1.0-exp(-alpha);

	/* *** don't forget tail -> head */
	FOR_EACH_TOGGLE(i){
		node=NODE(i);
		if(node!=0) break;
		cumchange=0.0;
		L2th=0;
		ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
		echange = 2*ochange + 1;
		/* step through outedges of head  */
		for(e = EdgetreeMinimum(nwp->outedges, head);
				(u = nwp->outedges[e].value) != 0;
				e = EdgetreeSuccessor(nwp->outedges, e)){
			if (EdgetreeSearch(MIN(u,tail), MAX(u,tail), nwp->outedges) != 0){
				L2th++;
				L2tu=ochange;
				L2uh=ochange;
				/* step through outedges of u */
				for(f = EdgetreeMinimum(nwp->outedges, u);
						(v = nwp->outedges[f].value) != 0;
						f = EdgetreeSuccessor(nwp->outedges, f)){
					if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
					if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
				}
				/* step through inedges of u */
				for(f = EdgetreeMinimum(nwp->inedges, u);
						(v = nwp->inedges[f].value) != 0;
						f = EdgetreeSuccessor(nwp->inedges, f)){
					if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
					if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
				}
				cumchange += pow(oneexpa,(double)L2tu) +
						pow(oneexpa,(double)L2uh) ;
			}
		}
		/* step through inedges of head */

		for(e = EdgetreeMinimum(nwp->inedges, head);
				(u = nwp->inedges[e].value) != 0;
				e = EdgetreeSuccessor(nwp->inedges, e)){
			if (EdgetreeSearch(MIN(u,tail), MAX(u,tail), nwp->outedges) != 0){
				L2th++;
				L2tu=ochange;
				L2uh=ochange;
				/* step through outedges of u */
				for(f = EdgetreeMinimum(nwp->outedges, u);
						(v = nwp->outedges[f].value) != 0;
						f = EdgetreeSuccessor(nwp->outedges, f)){
					if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
					if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
				}
				/* step through inedges of u */
				for(f = EdgetreeMinimum(nwp->inedges, u);
						(v = nwp->inedges[f].value) != 0;
						f = EdgetreeSuccessor(nwp->inedges, f)){
					if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
					if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
				}
				cumchange += pow(oneexpa,(double)L2tu) +
						pow(oneexpa,(double)L2uh) ;
			}
		}

		if(alpha < 100.0){
			cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2th)) ;
		}else{
			cumchange += (double)L2th;
		}
		cumchange  = echange*cumchange;
		(CHANGE_STAT[0]) += cumchange;
		TOGGLE_IF_MORE_TO_COME(i);
	}

	UNDO_PREVIOUS_TOGGLES(i);
}

