% test speed of saccade detection

max_seconds = 5;
sampling_rate = 500;
resampling_rate = 500;
n_samples = (1:max_seconds)*sampling_rate;
n_iter = 30000;
timings = NaN(numel(n_samples), n_iter);

for sample_i = 1:numel(n_samples)
    samples = n_samples(sample_i);
    for iter = 1:n_iter
        % data
        x = randn(1, samples);
        y = randn(1, samples);
        t = cumsum(repmat(1000/sampling_rate, 1, samples)); % simulate the time stamps
        % timing and algorithm
        t0 = GetSecs;
        [detected, t_detected, v_detected_x, v_detected_y, ...
            sac_thres_x, sac_thres_y, t_sac_onset] = ...
            detect_saccade_2d_C_ONSET( x, y, ...
            t, ...
            10, 5, ...
            0, 0, ... % no direction criterion
            resampling_rate, 10, 0);
        timings(sample_i, iter) = (GetSecs-t0)*1000;
    end
end

% analyze timing
n_samples' 
[median(timings, 2), std(timings, 0, 2)] % mean and sd run times across 1:max_seconds
median(timings, 2)*1000 ./ n_samples' % run time per sample in microseconds


