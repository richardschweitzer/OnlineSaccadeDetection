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
 *
 * Remarks: This version applied additional smoothing to the first and last velocity samples.
 */ 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include "detect_saccade_pure_C_ONSET.h"

/* QUICK SELECT FUNCTION TO COMPUTE MEDIAN: implementation by N. Devillard. */
/* retrieved from: http://ndevilla.free.fr/median/median/index.html */
/* Algorithm described in:
 * Numerical recipes in C, second edition, Cambridge University Press, 1992, section 8.5, ISBN 0-521-43108-5 */
#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }
double quick_select(double *arr, int n) { // previously arr[]
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
//    printf("x_here=%g y_here=%g deg_here=%g \n", x_here, y_here, deg_here);
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


/* here we wrap the detection results into a struct, so that we can return it */
struct detection_results {
    double sac_detected, sac_t, sac_vx, sac_vy; // results of saccade detection
    double threshold_vx, threshold_vy; // velocity threshold criteria
    double sac_t_onset; // actual saccade onset
};


/* main function here */
struct detection_results run_detection(double *x, double *y, double *t, 
                   unsigned int thres_fac, unsigned int above_thres_needed, 
                   unsigned int restrict_dir_min, unsigned int restrict_dir_max, 
                   unsigned int sampling_freq, unsigned int anchor_vel_thres, unsigned int print_results, 
                   unsigned int n_elements
                  ) {
    // output
    double sac_detected = 0, sac_t, sac_vx, sac_vy; // results of saccade detection
    double threshold_vx, threshold_vy; // velocity threshold criteria
    double sac_t_onset; // actual saccade onset
    struct detection_results results_struct;
    // others
	unsigned int ncols_x, ncols_y, ncols_t, ncols_inter, ncols_v, ncols_v_sq;     // size of vectors
    unsigned int it; // the iterator for our loops
    unsigned int do_run_detection; // is set to zero, if any error should occur
    

	// get dimensions of the input vectors and check that they all have the same length
    ncols_x = ncols_y = ncols_t = n_elements;
    
    // print some intro
    if (print_results != 0) {
        fprintf(stdout, "\n\n\nRunning Velocity-based Sac Detection by Richard Schweitzer with following parameters:\n");
        fprintf(stdout, "thres_fac=%i above_thres_needed=%i restrict_dir_min=%i restrict_dir_max=%i sampling_freq=%i anchor_vel_thres=%i print_results=%i n_samples=%i\n",
                thres_fac, above_thres_needed, restrict_dir_min, restrict_dir_max, sampling_freq, anchor_vel_thres, print_results, ncols_x);
    }
    
    // make sure that we have a long enough input vector
    if (ncols_t < (2*above_thres_needed)) {
        do_run_detection = 0;
        fprintf(stdout, "ERROR: You must supply a vector of a length that is at least twice of the number of samples needed above threshold!");
    } else {
        do_run_detection = 1; // ALL GOOD!
    }
    
    // transform the input vectors, in case they're nonuniformly sampled
    double time_passed = t[ncols_t-1] - t[0];
    double fd_uniform = time_passed / (double) ncols_t; // frame duration in milliseconds
    ncols_inter = ncols_t + 1; // length for the interpolated samples
    if (sampling_freq > 0) { // we interpolate based on the given sampling rate
        fd_uniform = 1000 / (double) sampling_freq; // frame duration in milliseconds
        ncols_inter = round(time_passed / fd_uniform) + 1; // length of the interpolated samples
    }
    if (print_results != 0) {
        fprintf(stdout, "Last t = %g, total time = %g\n", t[ncols_t-1], time_passed);
        fprintf(stdout, "Last x = %g, last y = %g\n", x[ncols_x-1], y[ncols_y-1]);
        fprintf(stdout, "Number of samples = %i\n", ncols_t);
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
//            fprintf(stdout, "it = %i, t_inter[it] = %g, t[kt] = %g, t[kt+1] = %g\n", it, t_inter[it], t[kt], t[kt+1]);
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
    if (print_results != 0) {
        for (it = 0; it < ncols_inter; it++) {
            if (it == 0) {
                fprintf(stdout, "it, t_inter, x_inter, y_inter\n");
            }
            fprintf(stdout, "%i, %g, %g, %g\n", it, t_inter[it], x_inter[it], y_inter[it]);
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
    double median_vx, median_vy, sd_vx, sd_vy; // values to compute the thresholds
    if (do_run_detection == 1) {
        ncols_v_sq = ncols_v - above_thres_needed; // number of elements without the ones that will be tested
        double vx_sq[ncols_v_sq], vy_sq[ncols_v_sq]; // median-cleaned and squared vectors
        median_vx = quick_select(vx, ncols_v_sq);
        median_vy = quick_select(vy, ncols_v_sq);
        for (it = 0; it < ncols_v_sq; it++) {
            vx_sq[it] = pow(vx[it] - median_vx, 2);
            vy_sq[it] = pow(vy[it] - median_vy, 2);
        }
        sd_vx = sqrt(quick_select(vx_sq, ncols_v_sq));
        sd_vy = sqrt(quick_select(vy_sq, ncols_v_sq));
        if (print_results != 0) {
            fprintf(stdout, "Median: x=%g y=%g\n",
                    median_vx, median_vy);
            fprintf(stdout, "SD: x=%g y=%g\n",
                    sd_vx, sd_vy);
        }
    } else { // we cannot compute thresholds, because we don't have enough samples!
        median_vx = NAN;
        median_vy = NAN;
        sd_vx = NAN;
        sd_vy = NAN;
    }
    
    // compute the thresholds
    threshold_vx = sd_vx * (double) thres_fac;
    threshold_vy = sd_vy * (double) thres_fac;
    if (print_results != 0) {
        fprintf(stdout, "Threshold: x=%g y=%g\n", threshold_vx, threshold_vy);
    }
    
    // compute the test criteria for the last few samples
    short int test_criterion, direction_criterion; 
    unsigned int above_thresholds = 0;
    unsigned int detected_at = 0; // position at which we have the needed number of samples
    double deg;
    if (do_run_detection == 0 || 
            isnan(threshold_vx) || isnan(threshold_vy) || 
            threshold_vx == 0 || threshold_vy == 0) { 
        // one of the thresholds is NaN, meaning we have recorded a blink
        sac_detected = NAN;
        sac_t = NAN;
        sac_vx = NAN;
        sac_vy = NAN;
        sac_t_onset = NAN;
//        if (print_results != 0) {
            fprintf(stdout, "->STOP. There are NaNs or zeros in computed thresholds!\n");
//        }
    } else {
        // let's try to detect a saccade in the last n samples specified by above_thres_needed
        for (it=ncols_v-above_thres_needed; it<ncols_v; it++) { 
            if (print_results != 0) {
                fprintf(stdout, "->TEST sample %i with ncols_v=%i\n", it+1, ncols_v);
            }
            // compute the Engbert & Mergenthaler test criterion
            if ( (pow(vx[it] * (1.0/threshold_vx), 2) + pow(vy[it] * (1.0/threshold_vy), 2)) > 1 ) {
                test_criterion = 1;
            } else {
                test_criterion = 0;
            }
            if (print_results != 0) {
                fprintf(stdout, "vx=%g vy=%g\n", vx[it], vy[it]);
                fprintf(stdout, "test_criterion=%i\n", test_criterion);
            }
            // compute the new direction criterion, if specified
            if (test_criterion == 1) {
                if (restrict_dir_min==0 || restrict_dir_max==0) {
                    direction_criterion = 1;
                } else {
                    deg = vec_to_degree(vx[it], vy[it]);
                    if (print_results != 0) {
                        fprintf(stdout, "deg=%g\n", deg);
                    }
                    // is direction still in our acceptable range?
                    if (restrict_dir_max > restrict_dir_min) { // any saccade
                        if ( deg < (double) restrict_dir_max && 
                                deg > (double) restrict_dir_min ) {
                            direction_criterion = 1;
                        } else {
                            direction_criterion = 0;
                        }
                    } else { // rightward saccade
                        if ( (deg >= 0.0 && deg < (double) restrict_dir_max) ||
                                (deg <= 360.0 && deg > (double) restrict_dir_min) ) {
                            direction_criterion = 1;
                        } else {
                            direction_criterion = 0;
                        }
                    }
                }
                if (print_results != 0) {
                    fprintf(stdout, "direction_criterion=%i\n", direction_criterion);
                }
                // has this sample succeeded at both criteria?
                if (direction_criterion == 1) {
                    above_thresholds += 1;
                }
                if (above_thresholds == above_thres_needed) {
                    detected_at = it;
                }
                if (print_results != 0) {
                    fprintf(stdout, "above_thresholds=%i\n", above_thresholds);
                }
            }
        }
        // did all samples succeed the test criteria? If so, it's a saccade!
        if (detected_at > 0) {
            if (print_results != 0) {
                fprintf(stdout, "SACCADE DETECTED!\n");
            }
            sac_detected = 1;
            sac_t = vt[detected_at];
            sac_vx = vx[detected_at];
            sac_vy = vy[detected_at];
        } else {
            if (print_results != 0) {
                fprintf(stdout, "NO SACCADE DETECTED.\n");
            }
            sac_detected = 0;
            sac_t = NAN;
            sac_vx = NAN;
            sac_vy = NAN;
        }
        
        // NEW: here we walk back in time to find the actual onset of the saccade.
        // i.e., we consider the beginning of the saccade as the first sample that
        // falls below the specified velocity threshold (formerly lambda=5).
        double anchor_thres_vx = sd_vx * (double) anchor_vel_thres;
        double anchor_thres_vy = sd_vy * (double) anchor_vel_thres;
        sac_t_onset = NAN;       
        if (sac_detected == 1) {
            if (print_results != 0) {
                fprintf(stdout, "Now checking for actual saccade onset with:\n");
                fprintf(stdout, "anchor_thres_vx=%g anchor_thres_vy=%g\n", anchor_thres_vx, anchor_thres_vy);
            }
            it = detected_at; // start searching backwards from detection point
            // as long as we haven't decided about the real saccade onset and
            // we still have samples to look at:
            while (isnan(sac_t_onset) && it > 0) { 
                // we decide for the onset of the saccade when we are below the regular E&M threshold
                if (print_results != 0) {
                    fprintf(stdout, "it=%i vx=%g vy=%g\n", it, vx[it], vy[it]);
                }
                if ( (pow(vx[it] * (1.0/anchor_thres_vx), 2) + pow(vy[it] * (1.0/anchor_thres_vy), 2)) < 1 
                        ) {
                    sac_t_onset = vt[it]; // we have found a sample below threshold
                    if (print_results != 0) {
                        fprintf(stdout, "Actual saccade onset detected.\n");
                    }
                } else {
                    it -= 1; // next sample backwards
                }
            }
        }
    } // end of saccade detection
    
    // RETURN THE RESULTS STRUCT:
    /* struct detection_results {
        double sac_detected, sac_t, sac_vx, sac_vy; // results of saccade detection
        double threshold_vx, threshold_vy; // velocity threshold criteria
        double sac_t_onset; // actual saccade onset }; */
    results_struct.sac_detected = sac_detected;
    results_struct.sac_t = sac_t;
    results_struct.sac_vx = sac_vx;
    results_struct.sac_vy = sac_vy;
    results_struct.threshold_vx = threshold_vx;
    results_struct.threshold_vy = threshold_vy;
    results_struct.sac_t_onset = sac_t_onset;
    return results_struct;
};
