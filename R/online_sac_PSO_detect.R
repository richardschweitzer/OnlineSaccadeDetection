dyn_load_detect_online <- function(path_to_online_sac_detect = "/home/richard/Dropbox/PROMOTION WORKING FOLDER/Projects/Online Saccade Detection/GITHUB/R", 
                                   fun_name = "detect_saccade_for_R.so") {
  # This function loads the compiled C function. 
  # To compile, run in console: R CMD SHLIB detect_saccade_for_R.c
  # Source: http://users.stat.umn.edu/~geyer/rc/
  path_to <- file.path(path_to_online_sac_detect, fun_name)
  print(paste("dyn.loading:", path_to))
  dyn.load(path_to)
  return(is.loaded("run_detection"))
}

detect_online_saccade <- function(thres_fac=5, above_thres_needed=1, restrict_dir_min=0, restrict_dir_max=0, 
                                  sampling_freq=0, anchor_vel_thres=5, print_results=FALSE, 
                                  x, y, t) {
  # This function runs the online saccade detection algorithm (Schweitzer & Rolfs, 2019)
  # and returns a data.frame of results.
  require(assertthat)
  assert_that(is.loaded("run_detection"), msg = "'run_detection' is not loaded!") # check if library is loaded
  # check whether all vectors have same length and are otherwise right
  assert_that(length(x)==length(y), length(x)==length(t), 
              msg = "Vectors x,y,t must have same length!")
  assert_that(thres_fac>=0, above_thres_needed>=0, sampling_freq>=0, anchor_vel_thres>=0, 
              msg = "All parameters must be unsigned integers!")
  assert_that(restrict_dir_min>=0 & restrict_dir_min<=360, 
              restrict_dir_max>=0 & restrict_dir_max<=360, 
              msg = "restrict_dir_min and _max must be values between 0 and 360 degrees!")
  assert_that(is.logical(print_results), msg = "print_results must be logical!")
  # run the function
  number_of_output_values = 8
  res <- .C("run_detection", 
            thres_fac = as.integer(thres_fac), above_thres_needed = as.integer(above_thres_needed), 
            restrict_dir_min = as.integer(restrict_dir_min), restrict_dir_max = as.integer(restrict_dir_max), 
            sampling_freq = as.integer(sampling_freq), anchor_vel_thres = as.integer(anchor_vel_thres), 
            print_results = as.integer(print_results), n_elements = as.integer(length(x)), 
            x = as.double(x), 
            y = as.double(y), 
            t = as.double(t), 
            detection_results_array = as.double(rep(0, number_of_output_values))
  )
  # return results
  df <- data.frame(sac_detected = res$detection_results_array[1], 
                   sac_t = res$detection_results_array[2], 
                   sac_vx = res$detection_results_array[3], 
                   sac_vy = res$detection_results_array[4], 
                   threshold_vx = res$detection_results_array[5], 
                   threshold_vy = res$detection_results_array[6], 
                   sac_t_onset = res$detection_results_array[7],
                   sac_direction = res$detection_results_array[8],
                   thres_fac = res$thres_fac, above_thres_needed = res$above_thres_needed, 
                   restrict_dir_min = res$restrict_dir_min, restrict_dir_max = res$restrict_dir_max, 
                   sampling_freq = res$sampling_freq, anchor_vel_thres = res$anchor_vel_thres, 
                   print_results = res$print_results, n_elements = res$n_elements, 
                   stringsAsFactors = FALSE)
  return(df)
}

vec_to_degree_C <- function(x, y) {
  # This function converts X and Y to a direction. 
  # Can be used to estimate saccade direction from the saccade's horizontal and vertical components.
  require(assertthat)
  assert_that(is.loaded("vec_to_degree_for_R"), msg = "'vec_to_degree_for_R' is not loaded!")
  res_deg <- .C("vec_to_degree_for_R", 
                x_here = as.double(x), y_here = as.double(y), deg_here = as.double(0))
  return(res_deg$deg_here)
}

check_for_PSO <- function(eye_x, eye_y, eye_t, em_table, 
                          direction_range = 50, velocity_thres = 5,
                          add_direction_too = FALSE) {
  # This function performs the detection of the post-saccadic oscillation (PSO) based on direction.
  # For each saccade in the EM-table, we estimate overall direction and detect the point 
  # where the direction inversion occurs and add it to the EM-table. This is the peak of the PSO.
  # by Richard Schweitzer

  # how many saccades were detected using the Engbert-kliegl algorithm?
  n_saccades_in_table <- nrow(em_table)
  # if the table is not empty
  if (!is.null(em_table) && n_saccades_in_table>=1) {
    # new column for our more conservative offset estimate
    em_table <- cbind(em_table, rep(NaN, n_saccades_in_table))
    if (add_direction_too) {
      em_table <- cbind(em_table, rep(NaN, n_saccades_in_table))
    }
    n_cols_in_table <- ncol(em_table)
    # check for every saccade
    for (current_saccade in 1:n_saccades_in_table) {
      # extract the offline detection results
      offline_onset <- em_table[current_saccade,1]
      offline_offset <- em_table[current_saccade,2]
      # what is the overall direction of the saccade?
      overall_sac_direction <- vec_to_degree_C(eye_x[offline_offset]-eye_x[offline_onset], 
                                               eye_y[offline_offset]-eye_y[offline_onset])
      # what is the direction criterion?
      opposite_direction_limits <- c(overall_sac_direction+direction_range, overall_sac_direction-direction_range) # max, min
      opposite_direction_limits[opposite_direction_limits>360] <- opposite_direction_limits[opposite_direction_limits>360] - 360
      opposite_direction_limits[opposite_direction_limits<0] <- opposite_direction_limits[opposite_direction_limits<0] + 360
      opposite_direction_limits[opposite_direction_limits==0] <- opposite_direction_limits[opposite_direction_limits==0] + 0.00001
      # where should we start looking for the inversion point?
      start_looking <- floor(offline_onset+(offline_offset-offline_onset)/2)
      # pre-allocate
      degrees <- rep(NaN, length(eye_x))
      detected_velocity <- rep(NaN, length(eye_x))
      detected_direction <- rep(NaN, length(eye_x))
      # extract the point, where direction changes by 180 degrees, i.e., likely the PSO
      for (it in start_looking:offline_offset) {
        # run once without direction criterion (purely velocity based)
        detect_results <- detect_online_saccade(x = eye_x[1:it], y = eye_y[1:it], t = eye_t[1:it], 
                                                thres_fac = velocity_thres,
                                                restrict_dir_min = 0, 
                                                restrict_dir_max = 0)
        detected_velocity[it] <- detect_results$sac_detected
        degrees[it] <- detect_results$sac_direction
        # run once with direction criterion (velocity + direction)
        detect_results <- detect_online_saccade(x = eye_x[1:it], y = eye_y[1:it], t = eye_t[1:it], 
                                                thres_fac = velocity_thres,
                                                restrict_dir_min = opposite_direction_limits[2], 
                                                restrict_dir_max = opposite_direction_limits[1])
        detected_direction[it] <- detect_results$sac_detected
      }
      # what's the first mismatch? (if there's any)
      if (any(!is.nan(detected_direction) & !is.nan(detected_velocity) & detected_velocity!=detected_direction)) {
        first_mismatch <- min(which(!is.nan(detected_direction) & !is.nan(detected_velocity) & detected_velocity!=detected_direction))
        if (add_direction_too) {
          em_table[current_saccade, n_cols_in_table-1] <- first_mismatch
          em_table[current_saccade, n_cols_in_table] <- overall_sac_direction
        } else {
          em_table[current_saccade, n_cols_in_table] <- first_mismatch
        } # add direction
      } # mismatch found
    } # across saccades
  } else {
    print(paste("A valid em_table must be provided. check_for_PSO received:"))
    print(em_table)
  }
  return(em_table)
}

### test this
# source the Engbert-Kliegl microsaccade detection algorithm
microsacc_where <- "~/Dropbox/PROMOTION WORKING FOLDER/Projects/Online Saccade Detection/GITHUB/R"
source(file.path(microsacc_where, "vecvel.R"))
source(file.path(microsacc_where, "microsacc.R"))
# get some arbitrary saccade data
load(file.path(microsacc_where, "saccade_data.rda"))
# detect the saccades using the Engbert-Kliegl algorithm
(em_table <- microsacc(x = as.matrix(saccade_data[ ,c("sac_x_raw", "sac_y_raw")]), 
                       VFAC = 5, MINDUR = 15, SAMPLING = 500) ) 
# plot X-dimension over time, as well as saccade onsets and offsets
plot(saccade_data$sac_t_raw, saccade_data$sac_x_raw) 
abline(v = saccade_data$sac_t_raw[em_table$table[,1]], col = "blue", lty = 2) 
abline(v = saccade_data$sac_t_raw[em_table$table[,2]], col = "blue", lty = 2) 
# now load the PSO detection
dyn_load_detect_online()
# test whether the online saccade detection algorithm works in R
(some_random_result <- detect_online_saccade(x = rnorm(30), y = rnorm(30), t = 1:30, 
                                             above_thres_needed = 2))
# run the PSO detection algorithm
(new_table <- check_for_PSO(eye_x = saccade_data$sac_x_raw, eye_y = saccade_data$sac_y_raw, 
                            eye_t = saccade_data$sac_t_raw, em_table = em_table$table, 
                            velocity_thres = 5, direction_range = 45, 
                            add_direction_too = TRUE))
# plot saccade offsets corrected for PSO
abline(v = saccade_data$sac_t_raw[new_table[,8]], col = "red", lty = 1) 

