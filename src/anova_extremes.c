/* SCCS @(#)anova.c	1.8 08/13/01  */
/*
** ALG 3/3/2012: One sided extremes from criteria
** "Data Mining Criteria for Tree-Based Regression and Classification",
** Andreas Buja & Yung-Seop Lee, 2001.
*/
#include <math.h>  //for max/min
#include "rpart.h"
#include "rpartproto.h"

static double *mean, *sums;
static double *wts;
static int *countn;
static int *tsplit;
static int high_or_low;  //ALG 3/21/2012 looking for high means (1) or low means (-1)

int anova_extremes_init(int n,        double *y[],  int maxcat, char **error,
	      double *parm, int *size,    int who,    double *wt)
{

	//collapsing isn't possible for this method
	rp.collapse_is_possible = 0;


    if (who==1 && maxcat >0) {
	graycode_init0(maxcat);
	countn  = (int *)ALLOC(2*maxcat, sizeof(int));
	tsplit  = countn + maxcat;
	mean   = (double *)ALLOC(3*maxcat, sizeof(double));
	wts    = mean + maxcat;
	sums   = wts + maxcat;
	}
    *size =1;

    //ALG 3/21/2012.  set high or low mean
    // 4/11/2012: adjusted to read high or low from parm
    if( *parm == -1){
    	high_or_low = -1;
    }
    else{
    	high_or_low = 1;
    }
    return(0);

    //4/2/2012
    //This ensures splitting without regard to complexity.
    //Complexity does not make sense for criteria that are min/max over left
    //and right nodes becasue the complexity number is defined to be path-independent.
    //rp.split_check_offset = -HUGE_VAL;

}

/*
** The anova_extremes evaluation function.  Direct copy of
** since now we use anovass eval function.
*/
void anova_extremes_eval(int n, double *y[], double *value, double *risk, double *wt)
			{
			    int i;
			    double temp = 0., twt = 0.; /* sum of the weights */
			    double mean, ss;

			    for (i = 0; i < n; i++) {
				temp += *y[i] * wt[i];
				twt += wt[i];
			    }
			    mean = temp / twt;

			    ss = 0;
			    for (i = 0; i < n; i++) {
				temp = *y[i] - mean;
				ss += temp * temp * wt[i];
			    }

			    *value = mean;
			    *risk = ss;
			}

/*
** The anova splitting function.  Find that split point in x such that
**  we have the maximum value over
**      max(mean_left, mean_right),
**  or alternatively the minimum over minimums.
**
**  ALG 4/11/2012: adjusted to take just the pre-computed penalty
**  ALG 4/19/2012: took away penalty--- done in bsplit
*/
void anova_extremes(int n,    double *y[],     double *x,     int nclass,
	   int edge, double *improve, double *split, int *csplit,
	   double myrisk, double *wt)
    {
    int i,j;
    double temp;
    double left_sum, right_sum;
    double left_wt, right_wt;
    int    left_n,  right_n;
    double parent_mean, most_extreme_so_far;
    int direction = LEFT;
    int where = 0;
    double left_mean, right_mean;


    /*
     * ALG 3/3/2012
    ** The improvement of a node is max(mean_L,mean_R) - mean_parent, when
    ** we look for max means. Alternatively
    **  -1*( min(mean_L, mean_R) - mean_parent)
    **  when looking for small means.
    **
    **  Accomodates weights.
    */

    right_wt =0;
    right_n  = n;
    right_sum =0;
    for (i=0; i<n; i++) {
    	right_sum += *y[i] * wt[i];
    	right_wt  += wt[i];
	}

    parent_mean = right_sum/right_wt;  //3/14/2012 for this function risk=node mean...
    right_mean = parent_mean;  //right_mean = node mean to start


    if (nclass==0) {   /* continuous predictor */
		left_sum=0;   /* No data in left branch, to start */
		left_wt =0;   left_n =0;

		//at the beginning, most extreme mean is parent mean
		//multiply by -1 (if looking for low means so we preserve checking for > most_extreme_so_far)
		most_extreme_so_far = parent_mean * high_or_low;

		//loop through each split point
		for (i=0; right_n>edge; i++) {
			left_wt += wt[i];
			right_wt -= wt[i];
			left_n++;
			right_n--;
			temp = (*y[i]) * wt[i];
			left_sum  +=temp;
			right_sum -=temp;

			//is this split allowed?  make sure node isn't too small
			if (x[i+1] != x[i] && left_n>=edge) {

				//if yes, calculate new means, then switch
				left_mean = left_sum/left_wt;
				right_mean = right_sum/right_wt;

				//pick maximum (or minimum)
				//3/21/2012: multiply by high or low to accomodate low means
				temp = fmax(high_or_low*left_mean, high_or_low*right_mean);

				//are we better than what we've seen so far??
				if (temp > most_extreme_so_far) {
					most_extreme_so_far=temp;
					where =i;

					//ALG 3/19/2012. remove due to left/right adjustment for high(low) means
					/*
					if (left_sum < right_sum) direction = LEFT;
							  else    direction = RIGHT;
					*/
				}
			}
		}//end of loop through all split points

		/*  ALG 3/21/2012
		 * If we're looking for small means, most_extreme_so_far is -1*smallest value,
		 * and so to get improvement we need to subtract -1*parent_mean. ie. looking for mins,
		 * lowest seen is -2, parent mean is 10 ==> improvement is 12.
		 */
		*improve =  (most_extreme_so_far - parent_mean*high_or_low);

		//ALG 3/3/2012
		if (most_extreme_so_far > parent_mean*high_or_low) {   /* found something */
			csplit[0] = direction;
			*split = (x[where] + x[where+1]) /2;
		}
	} /* end of continuous predictor */


    /*
    ** Categorical predictor
    */
    else {

		for (i=0; i<nclass; i++) {
			sums[i] =0;
			countn[i]=0;
			wts[i] =0;
		}

		/* rank the classes by their mean y value */
		for (i=0; i<n; i++) {
			j = (int)x[i] -1;
			countn[j]++;
			wts[j] += wt[i];
			sums[j] += (*y[i]) * wt[i]; //alg 3/5/2012, removed subtracting parent_mean
		}

		for (i=0; i<nclass; i++)  {
			if (countn[i] >0) {
				tsplit[i] = RIGHT;
				mean[i] = sums[i]/ wts[i];
			}
			else tsplit[i] = 0;
		}
		graycode_init2(nclass, countn, mean);

		/*
		** Now find the split that we want
		*/
		left_wt =0;
		left_sum=0;
		left_n = 0;
		right_sum = parent_mean*n;  //start with everything on the right side

		//3/21/2012 adjustment to deal with low means
		most_extreme_so_far = parent_mean * high_or_low;
		where =0;
		while((j=graycode()) < nclass) {
			tsplit[j] = LEFT;
			left_n += countn[j];
			right_n-= countn[j];
			left_wt += wts[j];
			right_wt-= wts[j];
			left_sum += sums[j];
			right_sum-= sums[j];
			if (left_n>=edge  &&  right_n>=edge) {
				//ALG 3/5/2012: changed from SS to means...

				//if yes, calculate new means, then switch
				left_mean = left_sum/left_wt;
				right_mean = right_sum/right_wt;

				//pick maximum
				//3/21/2012 adjustment for low means
				temp = fmax(left_mean*high_or_low, right_mean*high_or_low);

				if (temp > most_extreme_so_far) {
					most_extreme_so_far=temp;

					if ( left_mean > right_mean) {
						for (i=0; i<nclass; i++) csplit[i] = -tsplit[i];
					}
					else {
						for (i=0; i<nclass; i++) csplit[i] = tsplit[i];
					}
				}//end if(temp > best)
			}
		}//end of finding best split...

		//3/21/2012: adjusted for low means
		*improve = (most_extreme_so_far - parent_mean*high_or_low);      /* % improvement */
	} /* end of categorical predictor */
}

/*direct copy of anovapred for consistency */
double anova_extremes_pred(double *y, double *yhat)
{
    double temp = y[0] - *yhat;
    return temp * temp;
}
