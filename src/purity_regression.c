/* SCCS @(#)anova.c	1.8 08/13/01  */
/*
 * ALG 3/21/2012
** The four routines for splitting based on one-sided
** purity in a regression setting.  Follows Buja & Lee.
*/

#include <math.h>  //for max/min, infinity
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

static double *mean, *sums, *square_sums; //added square_sums
static double *wts;
static int *countn;
static int *tsplit;

int purity_regression_init(int n,        double *y[],  int maxcat, char **error,
	      double *parm, int *size,    int who,    double *wt)
 {
   
   rp.collapse_is_possible = 0;

    if (who==1 && maxcat >0) {
	graycode_init0(maxcat);
	countn  = (int *)ALLOC(2*maxcat, sizeof(int));
	tsplit  = countn + maxcat;
	mean   = (double *)ALLOC(3*maxcat, sizeof(double));
	wts    = mean + maxcat;
	sums   = wts + maxcat;
	square_sums = sums + maxcat; //ALG 3/21/2012. memory space for squared observations
	}
    *size =1;

    //ALG 3/22/2012
    double ss, temp, twt, ymean;
    int i;

    temp =0;
    twt  =0;   		/* sum of the weights */
    for (i=0; i<n; i++) {
	temp += *y[i] * wt[i];
	twt  += wt[i];
	}
    ymean = temp/twt;

    /*
    ss =0;
    for (i=0; i<n; i++) {
	temp = *y[i] - ymean;
	ss += temp*temp* wt[i];
	}
    rp.root_risk = ss/twt; //in this case we use sample variance instead of SS
	*/

    return(0);
}

/*
** Evaluation function.
*  made this active as of 3/18/2013.
*/
void purity_regression_eval(int n, double *y[], double *value, double *risk,
	     double *wt) {

    int i;
    double temp, twt;
    double nodemean, ss;

    temp =0;
    twt  =0;   		/* sum of the weights */

    for (i=0; i<n; i++) {
    	temp += *y[i] * wt[i];
    	twt  += wt[i];
	}
    nodemean = temp/twt;

    ss =0;
    for (i=0; i<n; i++) {
    	temp = *y[i] - nodemean;
    	ss += temp*temp* wt[i];
	}

    *value = nodemean;
    *risk = ss/twt;
}

/*
** The splitting function.  Find that split point in x such that
**   we minimize min(ssq_left , ssq_right)
**   where ssq_left = sample variance of left node.
**
**  ALG 4/30/2012: added parent_objective
*/
void purity_regression(int n,    double *y[],     double *x,     int nclass,
	   int edge, double *improve, double *split, int *csplit,
	   double myrisk, double *wt, double *parent_objective)
    {
    int i,j;
    double temp;
    double left_sum, right_sum;
    double left_wt, right_wt;
    int    left_n,  right_n;
    double parent_var, best;
    int direction = LEFT;
    int where = 0;
    double right_ssq, left_ssq;  //sum of squared observations left and right
    double right_var, left_var;  //right and left sample variances
    /*
    ** The improvement of a node is SSQ_parent - min(ssq_left, ssq_right)
    ** for the minimum min(left, right) value.
    */
    right_ssq = 0;
    right_wt =0;
    right_n  = n;
    right_sum =0;
    for (i=0; i<n; i++) {
    	right_sum += *y[i] * wt[i];
    	right_wt  += wt[i];
    	right_ssq += (*y[i]) * (*y[i]) * wt[i];
	}
    //sample variance of parent
    parent_var = (right_ssq/right_wt) - (right_sum/right_wt)*(right_sum/right_wt);

    //ALG 3/18/2013: if parent_var = 0, then NOTHING IS ALLOWED.
    if(parent_var==0){
		*improve=0;
    }
    //only try things if parent_var>0
    else{

	    if (nclass==0) {   /* continuous predictor */
		left_sum=0;   /* No data in left branch, to start */
		left_wt =0;   left_n =0;
		best  = parent_var;
		left_ssq = 0;
		for (i=0; right_n>edge; i++) {
		    left_wt += wt[i];
		    right_wt -= wt[i];
		    left_n++;
		    right_n--;
		    temp = (*y[i]) * wt[i];
		    left_sum  +=temp;
		    right_sum -=temp;
		    left_ssq  += temp * (*y[i]);
		    right_ssq -= temp * (*y[i]);

		    //alg 3/3/2012: check if we're allowed to do this split
		    if (x[i+1] !=x[i] &&  left_n>=edge) {
		    	right_var =(1/right_wt) *( right_ssq - (right_sum*right_sum / right_wt));
		    	left_var = (1/left_wt) * (left_ssq - (left_sum*left_sum / left_wt));
				temp = fmin(right_var, left_var); //minimum sample variance
				if (temp < best) {  //are we lower than the min so far?
					best=temp;
					where =i;
					if ( (left_sum/left_wt) < (right_sum/right_wt)) direction = LEFT;
							  else    direction = RIGHT;
					}
		    }
		 }//end loop through each split
		*improve =  (parent_var - best);
		if (best<parent_var) {   /* found something */
		    csplit[0] = direction;
		    *split = (x[where] + x[where+1]) /2;
		    }
		} /* end of continuous predictor */

	    /*
	    ** Categorical predictor
	    */
	    else {

	    //initialize
		for (i=0; i<nclass; i++) {
		    sums[i] =0;
		    countn[i]=0;
		    wts[i] =0;
		    square_sums[i]=0;
		}

		/* rank the classes by their mean y value */
		for (i=0; i<n; i++) {
		    j = (int)x[i] -1;
		    countn[j]++;
		    wts[j] += wt[i];
		    sums[j] += *y[i] * wt[i];
		    square_sums[j] += (*y[i]) * (*y[i]) * wt[i];
		}
		for (i=0; i<nclass; i++)  {
		    if (countn[i] >0) {
			tsplit[i] = RIGHT;
			mean[i] = sums[i]/ wts[i];
			}
		    else tsplit[i] = 0;
		    }

		//might not be the most useful ordering since criteria
		//is sample variance instead of mean.
		graycode_init2(nclass, countn, mean);

		/*
		** Now find the split that we want
		*/
	    best  = parent_var;
		left_wt =0;
		left_sum=0;
		left_n = 0;
	    left_ssq = 0;
		where =0;
		while((j=graycode()) < nclass) {
		    tsplit[j] = LEFT;
		    left_n += countn[j];
		    right_n-= countn[j];
		    left_wt += wts[j];
		    right_wt-= wts[j];
		    left_sum += sums[j];
		    right_sum -= sums[j];
		    left_ssq  += square_sums[j];
		    right_ssq -= square_sums[j];

		    if (left_n>=edge  &&  right_n>=edge) {  //we can do the split

		    	//compute sample var left and right, take the minimum
		    	right_var = (right_ssq/right_wt) - (right_sum/right_wt)*(right_sum/right_wt);
		    	left_var = (left_ssq/left_wt) - (left_sum/left_wt)*(left_sum/left_wt);
				temp = fmin(right_var, left_var); //minimum sample variance

		    	if (temp < best) {
		    		best=temp;
		    		if ( (left_sum/left_wt) > (right_sum/right_wt)) {
		    			for (i=0; i<nclass; i++) csplit[i] = -tsplit[i];
		    		}
		    		else {
		    			for (i=0; i<nclass; i++) csplit[i] = tsplit[i];
		    		}
			    }
			}
		  }
			*improve = (parent_var - best);
			//added 3/18/2013
			if (best<parent_var) {   /* found something */
		    	csplit[0] = direction;
			    *split = (x[where] + x[where+1]) /2;
		    }
		} /* end of categorical predictor */
	}
}

/*direct copy of anovapred for consistency */
double purity_regression_pred(double *y, double *yhat)
{
    double temp = y[0] - *yhat;
    return temp * temp;
}
