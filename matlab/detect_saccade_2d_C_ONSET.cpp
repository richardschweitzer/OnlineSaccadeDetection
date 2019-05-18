/* Velocity-based online saccade detection, by Richard Schweitzer.
 * Inspired by the approach of Engbert & Mergenthaler (2006).
 * Input:
 * - x, y gaze positions and corresponding timestamps (must be sorted!) in milliseconds
 * - threshold factor: is multiplied with the median-based standard deviation of velocity samples
 * - above thresholds needed: how many samples that pass the test criteria do we need to accept a saccade?
 * - restrict direction min and max: We can limit the direction of detected saccades to certain directions specified by two values.
 *   For example, 250 and 290 limit the detection to saccades that go upwards (given that 0/0 is in the upper left).
 *   If you don't want to use the direction criterion, set both values to zero. 
 * - sampling_freq: specify a sampling rate in Hz to resample to.
 * - anchor_vel_thres: for detection of actual saccade onset, how many SDs are the lower limit?
 * - print_results: If 0, then results are not printed and interpolation vectors are not returned
 *
 * Output:
 * - sac_detected: was saccade detected, i.e., 0 or 1
 * - sac_t, sac_vx, sac_vy: time and velocity when detection happened
 * - threshold_vx, threshold_vy: velocity thresholds applied (sd * threshold factor)
 * - sac_t_onset: actual saccade onset.
 * - t_inter_out, x_inter_out, y_inter_out: interpolated/resampled position vectors 
 *   (only if print_results was set to 1)
 *
 * Remarks: This version applied additional smoothing to the first and last velocity samples.
 */ 

#include "mex.h"
#include <matrix.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


/* QUICK SELECT FUNCTION TO COMPUTE MEDIAN: implementation by N. Devillard. */
/* retrieved from: http://ndevilla.free.fr/median/median/index.html */
/* Algorithm described in:
 * Numerical recipes in C, second edition, Cambridge University Press, 1992, section 8.5, ISBN 0-521-43108-5 */
#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }
double quick_select(double arr[], int n) {
    int low, high ;
    int median;
    int middle, ll, hh;
    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;
        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }
    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;
    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;
        if (hh < ll)
        break;
        ELEM_SWAP(arr[ll], arr[hh]) ;
    }
    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;
    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

/* VECTOR TO DEGREE FUNCTION */
double vec_to_degree(double x_here, double y_here) {
    double deg_here;
    deg_here = atan2(y_here, x_here)*180*1.0/M_PI;
    if (deg_here < 0) {
        deg_here = deg_here + 360;
    }
    return deg_here;
}

/* FIND ADJACENT VALUES IN VECTOR FUNCTION */
unsigned int find_adjacent_values(double val, unsigned int start_here, unsigned int end_here, double *vector) {
    unsigned int i = start_here, found = 0, found_at = -1;
    while (i < end_here-1 && found == 0) {
        if (val >= vector[i] && val < vector[i+1]) { // check whether we found adjacent values
            found = 1;
            found_at = i;
        } else {
            i += 1;
        }
    }
    return found_at;
}

/* MEXFUNCTION HERE */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // output
	double *sac_detected, *sac_t, *sac_vx, *sac_vy, // results of saccade detection
           *threshold_vx, *threshold_vy, // velocity threshold criteria
           *sac_t_onset, // actual saccade onset
           *t_inter_out, *x_inter_out, *y_inter_out; // interpolated values
    // input
	double *x, *y, *t, // input vectors   
            *thres_fac, *above_thres_needed, *restrict_dir_min, *restrict_dir_max, // input options
            *sampling_freq, *anchor_vel_thres, *print_results; 
    // others
	unsigned int ncols_x, ncols_y, ncols_t, ncols_inter, ncols_v, ncols_v_sq;     // size of vectors
    unsigned int it; // the iterator for our loops
		
	// check whether we have enough input arguments
	if (nrhs != 10) {
		mexErrMsgTxt("Input must be 10 arguments: x, y, t, thres_fac, above_thres_needed, restrict_dir_min, restrict_dir_max, sampling_freq, anchor_vel_thres, print_results!\n");
	}

	// make sure the second input argument is type double
    if ( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) ||
            !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]) || !mxIsDouble(prhs[5]) || !mxIsDouble(prhs[6]) ||
            !mxIsDouble(prhs[7]) || !mxIsDouble(prhs[8]) || !mxIsDouble(prhs[9]) ) {
        mexErrMsgTxt("All inputs must be type double.");
    }
    
	// get dimensions of the input vectors and check that they all have the same length
    ncols_x = mxGetN(prhs[0]);
	ncols_y = mxGetN(prhs[1]);
    ncols_t = mxGetN(prhs[2]);
    if ( ncols_x != ncols_t || ncols_y != ncols_t ) {
        mexErrMsgTxt("All input vectors (i.e., arguments 1-3) must have the same length.");
    }
    
    // check that we have multiple elements in row vectors
    if( mxGetM(prhs[0])!=1 || mxGetM(prhs[1])!=1 || mxGetM(prhs[2])!=1 || 
            ncols_x<=1 || ncols_y<=1 || ncols_t<=1) {
        mexErrMsgTxt("Input vectors (i.e., arguments 1-3) must be a row vectors with multiple elements.");
    }
    
    // check whether all input options have length 1
    if ( mxGetN(prhs[3])!=1 || mxGetN(prhs[4])!=1 || mxGetN(prhs[5])!=1 || 
            mxGetN(prhs[6])!=1 || mxGetN(prhs[7])!=1 || mxGetN(prhs[8])!=1 || 
            mxGetN(prhs[9])!=1 ||
            mxGetM(prhs[3])!=1 || mxGetM(prhs[4])!=1 || mxGetM(prhs[5])!=1 || 
            mxGetM(prhs[6])!=1 || mxGetM(prhs[7])!=1 || mxGetM(prhs[8])!=1 ||
            mxGetM(prhs[9])!=1 ) {
        mexErrMsgTxt("All option inputs (i.e., arguments 4-10) must be of length = 1.");
    }
    
	// now, read the input ...
	x = mxGetPr(prhs[0]);
	y = mxGetPr(prhs[1]);
    t = mxGetPr(prhs[2]);
    thres_fac = mxGetPr(prhs[3]);
    above_thres_needed = mxGetPr(prhs[4]); 
    above_thres_needed[0] = round(above_thres_needed[0]);
    restrict_dir_min = mxGetPr(prhs[5]); 
    restrict_dir_max = mxGetPr(prhs[6]);
    sampling_freq = mxGetPr(prhs[7]);
    if (sampling_freq[0] < 0) { // make sure sampling_freq is equal-larger than zero
        sampling_freq[0] = 0;
    }
    anchor_vel_thres = mxGetPr(prhs[8]);
    print_results = mxGetPr(prhs[9]);
    print_results[0] = round(print_results[0]);
    if (print_results[0] != 0) {
        fprintf(stderr, "\n\n\nRunning Velocity-based Sac Detection by Richard Schweitzer with:\n");
        fprintf(stderr, "thres_fac=%g above_thres_needed=%g restrict_dir_min=%g restrict_dir_max=%g sampling_freq=%g anchor_vel_thres=%g print_results=%g\n",
                thres_fac[0], above_thres_needed[0], restrict_dir_min[0], restrict_dir_max[0], sampling_freq[0], anchor_vel_thres[0], print_results[0]);
    }
    
    // make sure that we have a long enough input vector
    if (ncols_t < (2*above_thres_needed[0])) {
        mexErrMsgTxt("You must supply a vector of a length that is at least twice of the number of samples needed above threshold!");
    }
    
    // create the output variables and get a pointer to the real data in the output matrix
    if (print_results[0] != 0) {
        nlhs = 10;
    } else {
        nlhs = 7;
    }
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); // sac_detected
    sac_detected = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); // sac_t
    sac_t = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); // sac_vx
    sac_vx = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); // sac_vy
    sac_vy = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); // threshold_vx
    threshold_vx = mxGetPr(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL); // threshold_vy
    threshold_vy = mxGetPr(plhs[5]);
    plhs[6] = mxCreateDoubleMatrix(1,1,mxREAL); // sac_t_onset
    sac_t_onset = mxGetPr(plhs[6]);
    
    // transform the input vectors, in case they're nonuniformly sampled
    double time_passed = t[ncols_t-1] - t[0];
    double fd_uniform = time_passed / ncols_t; // frame duration in milliseconds
    ncols_inter = ncols_t + 1; // length for the interpolated samples
    if (sampling_freq[0] > 0) { // we interpolate based on the given sampling rate
        fd_uniform = 1000 / sampling_freq[0]; // frame duration in milliseconds
        ncols_inter = round(time_passed / fd_uniform) + 1; // length of the interpolated samples
    }
    if (print_results[0] != 0) {
        fprintf(stderr, "Last t = %g, total time = %g\n", t[ncols_t-1], time_passed);
        fprintf(stderr, "Last x = %g, last y = %g\n", x[ncols_x-1], y[ncols_y-1]);
        fprintf(stderr, "Number of samples = %i\n", ncols_t);
    }
    
    // make output vectors based on ncols_inter
    if (print_results[0] != 0) {
        plhs[7] = mxCreateDoubleMatrix(1,(mwSize)ncols_inter,mxREAL); // interpolated values T
        t_inter_out = mxGetPr(plhs[7]);
        plhs[8] = mxCreateDoubleMatrix(1,(mwSize)ncols_inter,mxREAL); // interpolated values X
        x_inter_out = mxGetPr(plhs[8]);
        plhs[9] = mxCreateDoubleMatrix(1,(mwSize)ncols_inter,mxREAL); // interpolated values Y
        y_inter_out = mxGetPr(plhs[9]);
    }
    
    // run through time and interpolate X and Y for all uniform timestamps
    double cum_t = t[0];
    unsigned int kt = 0;
    double t_inter[ncols_inter], x_inter[ncols_inter], y_inter[ncols_inter];
    for (it = 0; it < ncols_inter-1; it++) { 
        // make a uniform time series
        t_inter[it] = cum_t;
        cum_t += fd_uniform;
        // find adjacent values for interpolation
        kt = find_adjacent_values(t_inter[it], kt, ncols_t, t);
//        if (print_results[0] != 0) {
//            fprintf(stderr, "it = %i, t_inter[it] = %g, t[kt] = %g, t[kt+1] = %g\n", it, t_inter[it], t[kt], t[kt+1]);
//        }
        // linear interpolation here for X
        x_inter[it] = ((x[kt+1]-x[kt]) * (1.0/(t[kt+1]-t[kt]))) * (t_inter[it]-t[kt]) + x[kt];
        // linear interpolation here for Y
        y_inter[it] = ((y[kt+1]-y[kt]) * (1.0/(t[kt+1]-t[kt]))) * (t_inter[it]-t[kt]) + y[kt];
    }
    // the last interpolated sample is always the raw sample
    t_inter[ncols_inter-1] = t[ncols_t-1];
    x_inter[ncols_inter-1] = x[ncols_x-1];
    y_inter[ncols_inter-1] = y[ncols_y-1];
    // write out interpolated values, if needed and make interpolation results available to matlab
    if (print_results[0] != 0) {
        for (it = 0; it < ncols_inter; it++) {
            if (it == 0) {
                fprintf(stderr, "it, t_inter, x_inter, y_inter\n");
            }
            fprintf(stderr, "%i, %g, %g, %g\n", it, t_inter[it], x_inter[it], y_inter[it]);
            t_inter_out[it] = t_inter[it]; 
            x_inter_out[it] = x_inter[it]; 
            y_inter_out[it] = y_inter[it]; 
        }
    }
    
	// compute the velocity from uniformly sampled data
    ncols_v = ncols_inter-1;
    double vx_temp[ncols_v], vy_temp[ncols_v], vt[ncols_v]; // temporary arrays for velocity
	for (it = 0; it < ncols_v; it++) { // compute velocity here
        vx_temp[it] = (x_inter[it+1] - x_inter[it]) * (1.0/(t_inter[it+1] - t_inter[it]));
        vy_temp[it] = (y_inter[it+1] - y_inter[it]) * (1.0/(t_inter[it+1] - t_inter[it]));
		vt[it] = t_inter[it+1];
    }  
    
    // since data is now uniformly sampled, use a 5-point moving window here, just like Engbert & Mergenthaler
    double vx[ncols_v], vy[ncols_v]; // arrays for velocity (smooth)
    for (it = 0; it < ncols_v; it++) { 
        if (it==ncols_v-1) { // last value
            vy[it] = (vy_temp[it-1] + vy_temp[it] + vy_temp[it]) * (1.0/3.0);
            vx[it] = (vx_temp[it-1] + vx_temp[it] + vx_temp[it]) * (1.0/3.0);
        } else if (it == 0) { // first value
            vy[it] = (vy_temp[it] + vy_temp[it] + vy_temp[it+1]) * (1.0/3.0);
            vx[it] = (vx_temp[it] + vx_temp[it] + vx_temp[it+1]) * (1.0/3.0);
        } else if (it==ncols_v-2) { // second last value
            vy[it] = (vy_temp[it-2] + vy_temp[it-1] + vy_temp[it] + vy_temp[it+1] + vy_temp[it+1]) * (0.2);
            vx[it] = (vx_temp[it-2] + vx_temp[it-1] + vx_temp[it] + vx_temp[it+1] + vx_temp[it+1]) * (0.2);
        } else if (it == 1) { // second value
            vy[it] = (vy_temp[it-1] + vy_temp[it-1] + vy_temp[it] + vy_temp[it+1] + vy_temp[it+2]) * (0.2);
            vx[it] = (vx_temp[it-1] + vx_temp[it-1] + vx_temp[it] + vx_temp[it+1] + vx_temp[it+2]) * (0.2);
        } else { // all others
            vy[it] = (vy_temp[it-2] + vy_temp[it-1] + vy_temp[it] + vy_temp[it+1] + vy_temp[it+2]) * (0.2);
            vx[it] = (vx_temp[it-2] + vx_temp[it-1] + vx_temp[it] + vx_temp[it+1] + vx_temp[it+2]) * (0.2);
        }
    }
    
    // compute the median velocity and standard deviation for each x and y
    ncols_v_sq = ncols_v - above_thres_needed[0]; // number of elements without the ones that will be tested
    double median_vx, median_vy, sd_vx, sd_vy; // values to compute the thresholds
    double vx_sq[ncols_v_sq], vy_sq[ncols_v_sq]; // median-cleaned and squared vectors
    median_vx = quick_select(vx, ncols_v_sq);
    median_vy = quick_select(vy, ncols_v_sq);
    for (it = 0; it < ncols_v_sq; it++) {
        vx_sq[it] = pow(vx[it] - median_vx, 2);
        vy_sq[it] = pow(vy[it] - median_vy, 2);
    }
    sd_vx = sqrt(quick_select(vx_sq, ncols_v_sq));
    sd_vy = sqrt(quick_select(vy_sq, ncols_v_sq));
    if (print_results[0] != 0) {
        fprintf(stderr, "Median: x=%g y=%g\n",
                median_vx, median_vy);
        fprintf(stderr, "SD: x=%g y=%g\n",
                sd_vx, sd_vy);
    }
    
    // compute the thresholds
    threshold_vx[0] = sd_vx * thres_fac[0];
    threshold_vy[0] = sd_vy * thres_fac[0];
    if (print_results[0] != 0) {
        fprintf(stderr, "Threshold: x=%g y=%g\n", threshold_vx[0], threshold_vy[0]);
    }
    
    // compute the test criteria for the last few samples
    short int test_criterion, direction_criterion; 
    short int above_thresholds = 0;
    unsigned int detected_at = 0; // position at which we have the needed number of samples
    
    double deg;
    if (isnan(threshold_vx[0]) || isnan(threshold_vy[0]) || 
            threshold_vx[0] == 0 || threshold_vy[0] == 0) { 
        // one of the thresholds is NaN, meaning we have recorded a blink
        sac_detected[0] = NAN;
        sac_t[0] = NAN;
        sac_vx[0] = NAN;
        sac_vy[0] = NAN;
        sac_t_onset[0] = NAN;
        if (print_results[0] != 0) {
            fprintf(stderr, "->STOP. There are NaNs or zeros in computed thresholds!\n");
        }
    } else {
        // let's try to detect a saccade in the last n samples specified by above_thres_needed
        for (it=ncols_v-above_thres_needed[0]; it<ncols_v; it++) { 
            if (print_results[0] != 0) {
                fprintf(stderr, "->TEST sample %i with ncols_v=%i\n", it+1, ncols_v);
            }
            // compute the Engbert & Mergenthaler test criterion
            if ( (pow(vx[it] * (1.0/threshold_vx[0]), 2) + pow(vy[it] * (1.0/threshold_vy[0]), 2)) > 1 ) {
                test_criterion = 1;
            } else {
                test_criterion = 0;
            }
            if (print_results[0] != 0) {
                fprintf(stderr, "vx=%g vy=%g\n", vx[it], vy[it]);
                fprintf(stderr, "test_criterion=%i\n", test_criterion);
            }
            // compute the new direction criterion, if specified
            if (test_criterion == 1) {
                if (restrict_dir_min[0]==0 || restrict_dir_max[0]==0) {
                    direction_criterion = 1;
                } else {
                    deg = vec_to_degree(vx[it], vy[it]);
                    if (print_results[0] != 0) {
                        fprintf(stderr, "deg=%g\n", deg);
                    }
                    // is direction still in our acceptable range?
                    if (restrict_dir_max[0] > restrict_dir_min[0]) { // any saccade
                        if ( deg < restrict_dir_max[0] && 
                                deg > restrict_dir_min[0] ) {
                            direction_criterion = 1;
                        } else {
                            direction_criterion = 0;
                        }
                    } else { // rightward saccade
                        if ( (deg >= 0 && deg < restrict_dir_max[0]) ||
                                (deg <= 360 && deg > restrict_dir_min[0]) ) {
                            direction_criterion = 1;
                        } else {
                            direction_criterion = 0;
                        }
                    }
                }
                if (print_results[0] != 0) {
                    fprintf(stderr, "direction_criterion=%i\n", direction_criterion);
                }
                // has this sample succeeded at both criteria?
                if (direction_criterion == 1) {
                    above_thresholds += 1;
                }
                if (above_thresholds == above_thres_needed[0]) {
                    detected_at = it;
                }
                if (print_results[0] != 0) {
                    fprintf(stderr, "above_thresholds=%i\n", above_thresholds);
                }
            }
        }
        // did all samples succeed the test criteria? If so, it's a saccade!
        if (detected_at > 0) {
            if (print_results[0] != 0) {
                fprintf(stderr, "SACCADE DETECTED!\n");
            }
            sac_detected[0] = 1;
            sac_t[0] = vt[detected_at];
            sac_vx[0] = vx[detected_at];
            sac_vy[0] = vy[detected_at];
        } else {
            if (print_results[0] != 0) {
                fprintf(stderr, "NO SACCADE DETECTED.\n");
            }
            sac_detected[0] = 0;
            sac_t[0] = NAN;
            sac_vx[0] = NAN;
            sac_vy[0] = NAN;
        }
        
        // NEW: here we walk back in time to find the actual onset of the saccade.
        // i.e., we consider the beginning of the saccade as the first sample that
        // falls below the specified velocity threshold (formerly lambda=5).
        double anchor_thres_vx = sd_vx * anchor_vel_thres[0];
        double anchor_thres_vy = sd_vy * anchor_vel_thres[0];
        sac_t_onset[0] = NAN;       
        if (sac_detected[0] == 1) {
            if (print_results[0] != 0) {
                fprintf(stderr, "Now checking for actual saccade onset with:\n");
                fprintf(stderr, "anchor_thres_vx=%g anchor_thres_vy=%g\n", anchor_thres_vx, anchor_thres_vy);
            }
            it = detected_at; // start searching backwards from detection point
            // as long as we haven't decided about the real saccade onset and
            // we still have samples to look at:
            while (isnan(sac_t_onset[0]) && it > 0) { 
                // we decide for the onset of the saccade when we are below the regular E&M threshold
                if (print_results[0] != 0) {
                    fprintf(stderr, "it=%i vx=%g vy=%g\n", it, vx[it], vy[it]);
                }
                if ( (pow(vx[it] * (1.0/anchor_thres_vx), 2) + pow(vy[it] * (1.0/anchor_thres_vy), 2)) < 1 
                        ) {
                    sac_t_onset[0] = vt[it]; // we have found a sample below threshold
                    if (print_results[0] != 0) {
                        fprintf(stderr, "Actual saccade onset detected.\n");
                    }
                } else {
                    it -= 1; // next sample backwards
                }
            }
        }
        
    } // end of saccade detection
};
