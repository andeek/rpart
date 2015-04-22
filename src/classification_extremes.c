/* SCCS @(#)anova.c	1.8 08/13/01  */
/*
 *  THIS IS A PLACEHOLDER
 *
** ALG 3/3/2012: One sided Classification extremes criteria from
** "Data Mining Criteria for Tree-Based Regression and Classification",
** Andreas Buja & Yung-Seop Lee, 2001.
**
*/
#include <math.h>  //for max/min
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
static int *tsplit,
	         *countn;
static double *awt,
	            *rate;

static int class_of_interest;  //Ranges from 0,...,num_class-1, corresponding to 1,...,num_class

int classification_extremes_init(int n, double **y, int maxcat, char **error,
	      double *parm, int *size, int who, double *wt)
{
    int i, j, k;
    double temp;

    /* allocate memory  and setup losses */
    if (who==1) {
        numclass =0;   /*number of classes */
        for (i=0; i<n; i++){
        	if (*y[i] > numclass)  numclass = *y[i];
        }

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
				loss[k] = parm[numclass+k]; //should all be =1 but leave it as a placeholder
				temp += loss[k] * prior[i];
				aprior[i] += loss[k] * prior[i]; //hence aprior=prior.
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

    // ALG 10/22/2012:  set class of interest from the last double in 'parms'.
	// external classes are from 1,...k  internal classes
	// are from 0,...,(k-1).  Hence we subtract 1 from whatever was passed.
    /* ALG 4/12/2012: Choose which impurity function to use
     */
	  class_of_interest = parm[numclass + numclass*numclass + 1]-1;
    return(0);
}


/*
 * 10/22
 */
/*
** Prediction. straight copy of ginipred function
** in gini.c
*/
double classification_extremes_pred(double *y, double *pred)
    {
    int i, j;
    double temp;
    i = y[0] -1;
    j = *pred -1;
    temp = prior[i]*loss[i*numclass +j];
    return(temp);
}

/*
** 10/18/2012: Classification extremes eval: phat for class of interest.
*/
void classification_extremes_eval(int n, double **y, double *value, double *risk,
	     double *wt){
    int i, j, max = 0;
    double temp, dev = 0;
    double prob;
    double rwt = 0;
    int  rtot = 0;

   /*
    * count up number in each class,
    *   and P(T), the probability of reaching this branch of the tree
    */
    for (i = 0; i < numclass; i++)
      freq[i] = 0;
    temp = 0;
    for (i = 0; i < n; i++) {
	    j = (int) y[i][0] - 1;
	    freq[j] += wt[i];
	    temp += wt[i] * prior[j];
      
      rwt += aprior[j] * wt[i];    /*altered weight = class prior * case_weight */
    }
    prob = temp;                /* this is actually P(T)*n; R code will fix
				 * it up */

   /*
    * Now compute best class and its error
    */
    for (i = 0; i < numclass; i++) {    /* assume class i were the prediction */
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
    
    
    //10/18/2012
    //the risk is just the 1-phat(class_of_interest) at this node, adjusted for priors and
    //case weights. ....higher phat --> lower risk.
    *risk = (1 - freq[class_of_interest]/rwt);
}


/*
 * 10/18/2012
** The splitting function.  Find that split point in x such that
**  we have the maximum value over
**      max(phat(class of interest)_left, phat(class of interest)_right),
**
**     classification_extremes assumes no loss matrix (all =1), but
**     we leave it in the code in case I think of a way to modify this later.
**
*/
void classification_extremes(int n, double *y[], double *x, int numcat,
	   int edge, double *improve, double *split, int *csplit,
	   double my_risk, double *wt)
{
    int i, j, k;
    double lwt, rwt;
    int rtot, ltot;
    int direction = LEFT, where = 0;
    double best, temp, p;
    double lmean, rmean;    /* used to decide direction */

    //ALG: 10/18/2012
    double class_of_interest_prop; //weighted proportion of class of interest in the parent node.
    double left_prop, right_prop;
    left_prop = 0;
    right_prop = 0;


    //reset right & left counts of y classes
    for (i=0; i<numclass; i++) {
    	left[i] = 0;
    	right[i]= 0;
	  }
    lwt = 0;  
    rwt = 0;
    rtot= 0;  
    ltot = 0;

    //put everything to the right to start
    for (i = 0; i < n; i++) {
    	j = (int) *y[i] - 1;   //actual value
    	rwt += aprior[j] * wt[i];    /*altered weight = prior * case_weight */
    	right[j] += wt[i];
		  rtot++;
	  }
    //now rwt has total node weight including priors and observation weights.

    //ALG 10/18/2012: weighted proportion of class of interest in the current node.
    class_of_interest_prop = right[class_of_interest]*aprior[class_of_interest]/rwt;
    best = class_of_interest_prop; //what we need to beat in order to split.
    
    /*
    ** at this point we split into 2 disjoint paths
    */
    if (numcat >0) goto categorical;

    //cts predictor
    for (i = 0; rtot > edge; i++) { //as we increment i, we take from right and put in left.
    	j = (int) *y[i] - 1; //class of this observation
    	rwt -= aprior[j] * wt[i];
    	lwt += aprior[j] * wt[i];
    	rtot--;
    	ltot++;
    	right[j] -= wt[i];
    	left[j]  += wt[i];

		  if (x[i + 1] != x[i] &&  (ltot>=edge)) { //are we allowed to split here?
  			temp = 0; 
        lmean = 0; 
        rmean = 0;
  
  			//class of interest left and right
  			left_prop  = left[class_of_interest]*aprior[class_of_interest]/lwt;
  			right_prop = right[class_of_interest]*aprior[class_of_interest]/rwt;
  
  			//we want the max value..
  			temp = fmax(left_prop, right_prop);
  
  			//book-keeping for left and right decision
  			for (j=0; j<numclass; j++) { //now interate through the classes
  				p = aprior[j]*left[j]/lwt;    /* p(j | left) */
  				lmean += p*j;
  				p =  aprior[j]*right[j]/rwt;   /* p(j | right) */
  				rmean += p*j;
  			}
  			if (temp > best) {
  					best = temp;
  					where = i;
  					if (lmean < rmean){
  						direction = LEFT;
  					}
  					else{
  						direction = RIGHT;
  					}
  				}
			}
		}//end of for loop

    //To redefine as impurity where higher is worse, we set it to 1-p_hat, and so
    //to fit in the previous context: (1-parent_prop) - (1- best_prop)
    *improve =  (1 - class_of_interest_prop) - (1 - best);
    if (*improve > 0 ) {   /* found something */
    	csplit[0] = direction;
    	*split = (x[where] + x[where+1]) /2;
	  }
    return; //end of cts predictor.

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
    	j = *y[i] -1;
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

    	if (ltot>=edge  &&  rtot>=edge) { //ALG: adjustments can go here.
    	    temp =0;
    	    lmean=0; rmean =0;

	        //figure out the proportion of class of interest left and right
      		left_prop  = left[class_of_interest]*aprior[class_of_interest]/lwt;
      		right_prop = right[class_of_interest]*aprior[class_of_interest]/rwt;

      		//we want the max value..
      		temp = fmax(left_prop, right_prop);

    	    //right left book-keeping
    	    for (j=0; j<numclass; j++) {
    			  p = aprior[j]*left[j] /lwt;
    			  lmean += p*j;
    			  p =  aprior[j]*right[j]/rwt;       /* p(j | right) */
    			  rmean += p*j;
    	    }
    	    if (temp > best) {
    	    	best=temp;
		        if (lmean < rmean)
			        for (j=0; j<numcat; j++) csplit[j] = tsplit[j];
	        	else
			        for (j=0; j<numcat; j++) csplit[j] = -tsplit[j];
	        }
	    }//end of loop through observations.
    }//end of while loop

    //set the improve just as before.
    *improve = (1- class_of_interest_prop) - (1 - best);
}
