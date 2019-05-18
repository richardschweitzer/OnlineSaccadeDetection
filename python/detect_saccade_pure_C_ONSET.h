/* header file for detect_saccade_pure_C_ONSET.cpp */

double quick_select(double arr[], int n);

double vec_to_degree(double x_here, double y_here);

unsigned int find_adjacent_values(double val, unsigned int start_here, unsigned int end_here, double *vector);

struct detection_results run_detection(double *x, double *y, double *t, 
                   unsigned int thres_fac, unsigned int above_thres_needed, 
                   unsigned int restrict_dir_min, unsigned int restrict_dir_max, 
                   unsigned int sampling_freq, unsigned int anchor_vel_thres, unsigned int print_results, 
                   unsigned int n_elements
                  );

