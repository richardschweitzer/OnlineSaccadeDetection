%% short script to elucidate the usage the online saccade detection algorithm in matlab.
% first, compile to mex function: mex detect_saccade_2d_C_ONSET.cpp
% by Richard Schweitzer

% parameters to play around with
how_many_millisec_after_sac_onset = 4; % 0 is saccade onset
thres_fac = 10;
above_thres_needed = 3;
restrict_dir_min = 0;
restrict_dir_max = 0;
anchor_thres_fac = 5;
print_results = 1;

% get some test data and 
% make a subset up to some samples after offline detected saccade onset
load('example_saccade_2.mat')
subset_here = 20:min(find(sac.t_sac_on>=how_many_millisec_after_sac_onset));
t_eye = sac.t_sac_on(subset_here);
t_eye = t_eye';
x = sac.sac_x(subset_here);
x = x';
y = sac.sac_y(subset_here);
y = y';

% sampling rate?
effective_fd = (t_eye(end)-t_eye(1)) / length(t_eye);
effective_samp_freq = 1000 / effective_fd;
resample_to_freq = 1000;

%% Run the version with saccade onset detection
[detected, d_t, d_vx, d_vy, th_vx, th_vy, d_sac_on, inter_t, inter_x, inter_y] = ...
    detect_saccade_2d_C_ONSET( x, y, t_eye,... % x, y, and time stamps [ms]
    thres_fac, above_thres_needed, restrict_dir_min, restrict_dir_max,... 
    resample_to_freq, anchor_thres_fac, print_results);     % re/sampling frequency, anchor velocity threshold, print results
disp(['detected=', num2str(detected), ' at t=', num2str(d_t) ' with thresholds of: ', num2str([th_vx, th_vy])])


%% make a figure of the resampling results
figure(1000);
subplot(1,2,1) % position x
plot(t_eye, x, '-o', inter_t, inter_x, '-x');
line([d_sac_on, d_sac_on], [min(x), max(x)], 'LineStyle','-');
subplot(1,2,2) % position y
plot(t_eye, y, '-o', inter_t, inter_y, '-x');
line([d_sac_on, d_sac_on], [min(y), max(y)], 'LineStyle','-');
