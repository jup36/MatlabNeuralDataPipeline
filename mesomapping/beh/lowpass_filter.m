function filtered_ts = lowpass_filter(time_series, samp_rate, cutoff_freq)
    % Apply a low-pass filter to the input time series
    % Inputs:
    %   time_series  - the input time series (vector)
    %   samp_rate    - the sampling rate (in Hz)
    %   cutoff_freq  - the cutoff frequency for the low-pass filter (in Hz)
    % Output:
    %   filtered_ts  - the filtered time series (vector)
    
    % Calculate the normalized cutoff frequency (cutoff_freq / Nyquist frequency)
    nyquist_freq = samp_rate / 2;
    normalized_cutoff = cutoff_freq / nyquist_freq;
    
    % Design a Butterworth low-pass filter of order 4
    [b, a] = butter(4, normalized_cutoff, 'low');
    
    % Apply the filter to the time series using filtfilt for zero-phase distortion
    filtered_ts = filtfilt(b, a, time_series);
end
