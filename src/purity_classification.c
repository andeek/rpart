/* SCCS @(#)anova.c	1.8 08/13/01  */
/*
 * ALG 9/19/2012
** The four routines for splitting based on one-sided
** purity in a classification setting.  Follows Buja & Lee,
** with extension beyond 2-class problems.
**
** A modification of gini.c, but we only use Gini, following
** Buja & Lee.
*/

#include <math.h>
#include "rpart.h"
#include "rpartproto.h"

static int    numclass;
static double *left,     /*left branch n (weighted)*/
	            *right,
	            **ccnt;
static double *prior,
	            *aprior,   /*altered priors */
	            *freq,     //alg.  freq is supposed to hold the weighted number of observ. for each class
              *loss;      /* loss matrix */
static int    *tsplit,
	            *countn;
static double *awt,
	            *rate;

static double (*impurity)();


//ALG 4/12/2012: Gini: squared error loss on +/-1 - fitted prob
//9/19/2012: Now we use p(1-p)
static double gini_purity_impure(p) double p; {  return(p - p*p); }


//9/19/2012
//Init function.  Mostly the same as the regular Gini's.
int purity_classification_init(int n,        double **y, int maxcat, char **error,
	     double *parm, int *size,  int who,    double *wt)
{
  
    rp.collapse_is_possible = 0;

    int i, j, k;
    double temp;

    /* allocate memory  and setup losses */
    if (who == 1) {
        numclass = 0;   /*number of classes */
        for (i=0; i<n; i++){
        	if (*y[i] > numclass)  numclass = *y[i];
        }

        //For now we follow B&L and only use Gini.
        impurity = gini_purity_impure;

        left = (double *) ALLOC(numclass*2, sizeof(double));
        right = left+numclass;

        tsplit= (int *) ALLOC(maxcat*2, sizeof(int));
        countn= tsplit + maxcat;

		awt = (double *) ALLOC(maxcat*2, sizeof(double));
		rate= awt +maxcat;

		if (maxcat>0) {
			graycode_init0(maxcat);
			ccnt    = (double **) ALLOC(numclass, sizeof(double *));
			if (ccnt==0) {*error=_("Out of memory"); return(1);}
			ccnt[0] = (double *) ALLOC(numclass*maxcat, sizeof(double));
			if (ccnt[0]==0) {*error=_("Out of memory"); return(1);}
			for (i=1; i<numclass; i++)
			ccnt[i] = ccnt[i-1] + maxcat;
		}

		i = 3*numclass + numclass*numclass;
		prior = (double *) ALLOC(i, sizeof (double));
		if (prior==0) {*error=_("Out of memory"); return(1);}
		aprior = prior + numclass;
		freq   = aprior+ numclass;
		loss   = freq  + numclass;

		for (i=0; i<numclass; i++)  freq[i] =0;
		temp =0;
		for (i=0; i<n; i++) {
			j = *y[i] -1;
			freq[j] += wt[i];
			temp += wt[i];   /*sum total of weights */
	    }
		for (i=0; i<numclass; i++)  freq[i] /=temp;   /*relative frequency */

		temp =0;
		for (i=0; i<numclass; i++) {
			prior[i] = parm[i];
			aprior[i] =0;
			for (j=0; j<numclass; j++) {
				k = numclass*i + j;
				loss[k] = parm[numclass+k];
				temp += loss[k] * prior[i];
				aprior[i] += loss[k] * prior[i];
			}
	    }

		for (i=0; i<numclass; i++) {
			if (freq[i]>0) {  /* watch out for a missing class */
				prior[i] /= freq[i];
				aprior[i] /= (temp * freq[i]);  /* pi_i / n_i */
			}
		}
	}//end if(who==1)

	*size = 2 + numclass;
	return(0);
}

/*
** Prediction. straight copy of ginipred function
** in gini.c
*/
double purity_classification_pred(double *y, double *pred)
    {
    int i, j;
    double temp;
    i = y[0] -1;
    j = *pred -1;
    temp = prior[i]*loss[i*numclass +j];
    return(temp);
}


void purity_classification_eval(int n, double **y, double *value, double *risk,
	     double *wt){

    int i, j, max = 0;
    double total_ss, temp, dev = 0;
    double prob;
    double rwt = 0;
    int  rtot = 0;

   /*
    * count up number in each class,
    *   and P(T), the probability of reaching this branch of the tree
    */
    for (i = 0; i < numclass; i++) {
      freq[i] = 0;
      right[i] = 0;
    }
      
    temp = 0;
    for (i = 0; i < n; i++) {
	    j = (int) y[i][0] - 1;
	    freq[j] += wt[i];
	    temp += wt[i] * prior[j];
      right[j] += wt[i];
      rwt += aprior[j] * wt[i];
    }
    prob = temp;                /* this is actually P(T)*n; R code will fix
				 * it up */

   /*
    * Now compute best class and its error
    */
    total_ss = 0;
    for (i = 0; i < numclass; i++) {    /* assume class i were the prediction */
    	temp = aprior[i] * right[i]/ rwt; 
      total_ss += rwt * (*impurity)(temp);
      
      temp = 0;
    	for (j = 0; j < numclass; j++)
    	    temp += freq[j] * loss[i * numclass + j] * prior[j];
    	if (i == 0 || temp < dev) {
    	    max = i;
    	    dev = temp;
    	}
    }

    value[0] = max + 1;         /* remember: external groups start at 1 */
    for (i = 0; i < numclass; i++)
	    value[i + 1] = freq[i];
    value[numclass + 1] = prob;
    


    //9/19/2012
    //Here we differ from Gini's parent_objective function in that
    //divide by total weight to make it a 'per observation' value.
    *risk = total_ss/rwt;
}


/*
 * ALG 9/19/2012
** The gini_purity splitting function.  Find that split point in x such that
**  we maximize
**  node(gini) - min(gini_left, gini_right) where each gini is on a per-observation
**  (ie indep of nodesize) level.
**  as possible.
*/
void purity_classification(int n,    double *y[],     double *x,     int numcat,
	  int edge, double *improve, double *split, int *csplit, double my_risk,
	  double *wt, double *parent_objective)
    {
    int i,j,k;
    double lwt, rwt;
    int rtot, ltot;
    int direction = LEFT, where = 0;
    double total_ss, best, temp, p, parent_gini;
    double left_impure, right_impure;  //ALG 9/19/2012 need these so we can take min.

    double lmean, rmean;    /* used to decide direction */

    //reset right & left counts of y classes
    for (i=0; i<numclass; i++) {
    	left[i] =0;
    	right[i]=0;
	  }
    lwt =0;  rwt=0;
    rtot=0;  ltot=0;

    //put everything to the right to start
    for (i=0; i<n; i++) {
    	j = *y[i] -1;   //actual value
    	rwt += aprior[j] * wt[i];    /*altered weight = prior * case_weight */
    	right[j] += wt[i];
  		rtot++;
  	}

    total_ss =0;
    for (i=0; i<numclass; i++){
    	temp = aprior[i] * right[i]/ rwt;      /* p(class=i, given node A) */
    	total_ss += rwt * (*impurity)(temp);    /* p(A) * I(A). this is total node weight * impur fcn */
	  }
    best = total_ss/rwt; /* ALG: DIVDE BY rwt!*/
    parent_gini = best; /* ALG: set parent Gini*/
    /*
    ** at this point we split into 2 disjoint paths
    */
    if (numcat >0) goto categorical;

    //cts predictor.
    //iterate through observations.
    for (i=0;  rtot >edge; i++) {

    	//bookkeeping with weights.
  		j = (int)*y[i] -1;  //the ith observation's class.
  		rwt -= aprior[j] * wt[i];
  		lwt += aprior[j] * wt[i];
  		rtot--;
  		ltot++;
  		right[j] -= wt[i];  //right[j] is weight of jth class in right node.
  		left[j]  += wt[i];

		//check to see if this is a new class and
		//we meet the min bucket size.
		if (x[i+1] != x[i] &&  (ltot>=edge)) {
			temp =0;
			lmean =0; rmean =0;
			left_impure = 0;
			right_impure = 0;

			//iterate through the classes
			// ALG: unless set otherwise, all aprior=1.
			// unless set otherwise, all wt=1, in which case left[j]=#(left node == class j)
			for (j=0; j<numclass; j++) {
				//LEFT
				p = aprior[j]*left[j]/lwt;    /* p(j | left) */
				left_impure += lwt * (*impurity)(p);      /* p(left) * I(left) */
				lmean += p*j;  //for figuring class proportions, for RIGHT or LEFT and visual displays

				p =  aprior[j]*right[j]/rwt;   /* p(j | right) */
				right_impure += rwt * (*impurity)(p);      /*p(right) * I(right) */
				rmean += p*j; //for figuring class proportions, for RIGHT or LEFT and visual displays
			}

			// ALG 9/19/2012: divide the impurities by their node-weights.
			//Set 'temp' to the minimum such value.
			left_impure = left_impure/lwt;
			right_impure = right_impure/rwt;
			temp = fmin(left_impure, right_impure);
			

			//Proceed as in gini split
			if (temp < best) {
				best = temp;
				where = i;
				if (lmean < rmean) direction = LEFT;
					else           direction = RIGHT;
			}
		}//end of if
	}//end for loop through cts predictor

    //this is the impurity of parent - min(weight(right)*Impurity(right)+weight(left)*Impurity(left)
    //The weights are btwn (0,n) not (0,1).
    *improve =  (parent_gini - best);
    if ( *improve > 0 ) {   /* found something */
    	csplit[0] = direction;
    	*split = (x[where] + x[where+1]) /2;
	  }
    return;


    categorical:;
    /*
    ** First collapse the data into a numclass x numcat array
    **  ccnt[i][j] = number of class i obs, category j of the predictor
    */
    for (j=0; j<numcat; j++) {
    	awt[j] =0;
    	countn[j]=0;
    	for (i=0; i<numclass; i++)
    		ccnt[i][j] =0;
	}

    for (i=0; i<n; i++) {
    	j = (int)*y[i] -1;
    	k = x[i] -1;
    	awt[k] += aprior[j] * wt[i];
    	countn[k]++;
    	ccnt[j][k] += wt[i];
	}

    for (i=0; i<numcat; i++){
    	if (awt[i]==0) tsplit[i] =0;
    	else {
    		rate[i] = ccnt[0][i] / awt[i];   /* a scratch array */
    		tsplit[i]=RIGHT;
	    }
    }

    if (numclass==2) graycode_init2(numcat, countn, rate);
                else graycode_init1(numcat, countn);

    //main loop for categorical
    while((i=graycode()) < numcat) {
	/* item i changes groups */
		if (tsplit[i]==LEFT) {
			tsplit[i]=RIGHT;
			rwt  += awt[i];
			lwt -= awt[i];
			rtot += countn[i];
			ltot -= countn[i];
			for (j=0; j<numclass; j++) {
				right[j] += ccnt[j][i];
				left[j]  -= ccnt[j][i];
			}
		}
		else {
			tsplit[i]=LEFT;
			rwt -= awt[i];
			lwt += awt[i];
			rtot -= countn[i];
			ltot += countn[i];
			for (j=0; j<numclass; j++) {
				right[j] -= ccnt[j][i];
				left[j]  += ccnt[j][i];
			}
		}

		if (ltot>=edge  &&  rtot>=edge) {
			temp =0;
			lmean=0; rmean =0;

			right_impure = 0;  //ALG 9/19/2012
			left_impure = 0;   //ALG  9/19/2012

			for (j=0; j<numclass; j++) {
				p = aprior[j]*left[j] /lwt;
				left_impure +=  lwt * (*impurity)(p);
				lmean += p*j;

				p =  aprior[j]*right[j]/rwt;       /* p(j | right) */
				right_impure += rwt * (*impurity)(p);      /*p(right) * I(right) */
				rmean += p*j;
		    }

			//ALG 9/19/2012
			//just like above, divide by weights, and set min = temp.
			left_impure = left_impure/lwt;
			right_impure = right_impure/rwt;
			temp = fmin(left_impure, right_impure);

			//proceed as before...
			if (temp < best) {
				best=temp;
				if (lmean < rmean)
					for (j=0; j<numcat; j++) csplit[j] = tsplit[j];
				else
					for (j=0; j<numcat; j++) csplit[j] = -tsplit[j];
			}
		}
    }//end of while loop
    *improve = (parent_gini - best);
}
