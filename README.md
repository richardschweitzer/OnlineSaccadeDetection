# Online Saccade Detection

Code for the online saccade detection algorithm, as described in:

An adaptive algorithm for fast and reliable online saccade detection\
Richard Schweitzer & Martin Rolfs\
Behavior Research Methods; doi: https://doi.org/10.3758/s13428-019-01304-3 \
bioRxiv 693309; doi: https://doi.org/10.1101/693309 

Two versions are available, one mex-function for usage in Matlab, one C-function along with a Python module (using ctypes). Compiling the mex-function has been successfully tested with Matlab2015b and Matlab2016b on K/Ubuntu 18.04 using g++ and with Matlab2019b on Windows 10 using the MinGW64 compiler installed via Matlab's Add-on manager. 

![real_saccade_simulation_1.png](https://raw.githubusercontent.com/richardschweitzer/OnlineSaccadeDetection/master/python/real_saccade_simulation_1.png)
Example figure created from online_sac_detect_module.py. Samples of a real saccade (saccade onset marked by dashed line) retrieved sequentially. Saccade detection was run after each retrieval (violet: not detected, yellow: detected). 

**Update:** The function is now also available in R and along with it a solution for the direction-based detection of post-saccadic oscillation (PSO). If you are interested in more conservative estimates of saccade duration that are not confounded by PSOs, then check out *online_sac_PSO_detect.R*, as the code is a convenient add-on for the Engbert-Kliegl saccade detection algorithm. The figure below shows a saccade sequence with saccade onsets and offsets determined by the former (blue dashed lines). Red solid lines indicate the more conservative estimates for saccade offset that take into account PSOs based on a transient direction inversion (of course only if the latter was detected).
![saccade_data_result.png](https://raw.githubusercontent.com/richardschweitzer/OnlineSaccadeDetection/master/R/saccade_data_result.png)
