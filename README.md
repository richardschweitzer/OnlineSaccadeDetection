# Online Saccade Detection

Code for the online saccade detection algorithm, as described in:

An adaptive algorithm for fast and reliable online saccade detection\
Richard Schweitzer & Martin Rolfs\
Behavior Research Methods; doi: https://doi.org/10.3758/s13428-019-01304-3 \
bioRxiv 693309; doi: https://doi.org/10.1101/693309 

Two versions are available, one mex-function for usage in Matlab, one C-function along with a Python module (using ctypes).

![real_saccade_simulation_1.png](https://raw.githubusercontent.com/richardschweitzer/OnlineSaccadeDetection/master/python/real_saccade_simulation_1.png)
Example figure created from online_sac_detect_module.py. Samples of a real saccade (saccade onset marked by dashed line) retrieved sequentially. Saccade detection was run after each retrieval (violet: not detected, yellow: detected). 

