function [out_timeseries] = digit_decompose(digit_norm, nSamp)
    % Decompose a time series with sums of 2^0, 2^1, ..., 2^7 into 8 time series
    % and filter out values that last shorter than nSamp/10.
    % Inputs:
    %   digit_norm - input time series (vector)
    %   nSamp - sampling rate
    % Outputs:
    %   out_timeseries - cell array of 8 time series vectors (binary), where each 
    %                    indicates the presence of each respective power of 2.
    
    % Define the minimum duration (in samples) that a value must last
    min_duration = round(nSamp / 200); % to cutoff short noises that last less than 5ms
    
    % Initialize a filtered version of digit_norm
    digit_filtered = digit_norm;
    
    % Find changes in the input time series
    changes = [1, find(diff(digit_norm) ~= 0) + 1, length(digit_norm) + 1];
    
    % Filter out any segment that lasts shorter than min_duration
    for ii = 1:(length(changes) - 1)
        segment_start = changes(ii);
        segment_end = changes(ii+1) - 1;
        if (segment_end - segment_start + 1) < min_duration
            % Set the values of this short segment to 0 (treated as noise)
            digit_filtered(segment_start:segment_end) = 0;
        end
    end
    
    % Define the 8 powers of 2 (2^0, 2^1, ..., 2^7)
    powers_of_2 = 2.^(0:7);
    
    % Initialize output timeseries as a matrix with 8 rows (one for each power of 2)
    out_timeseries = zeros(8, length(digit_filtered));
    
    % Initialize a variable to store the previous decomposition state
    previous_decomposition = zeros(8, 1);
    
    % Iterate only at change points
    for ii = 1:(length(changes) - 1)
        segment_start = changes(ii);
        segment_end = changes(ii+1) - 1;
        value = digit_filtered(segment_start);
        
        if value == 0
            % If value is zero, all powers of 2 are zero
            current_decomposition = zeros(8, 1);
        elseif all(ismember(dec2bin(value), ['0', '1'])) && value >= 0 && value < 256
            % Decompose the value only if it's valid
            current_decomposition = bitand(value, powers_of_2) ~= 0;
        else
            % Treat it as noise, skip decomposition
            current_decomposition = zeros(8, 1);
        end
        
        % Assign the current decomposition to the timeseries for the current segment
        for t = segment_start:segment_end
            if t == segment_start
                % For the first point in the segment, use current decomposition
                out_timeseries(:, t) = current_decomposition;
            else
                % For subsequent points in the segment, inherit the previous values
                out_timeseries(:, t) = out_timeseries(:, t-1);
            end
        end
    end
    
    % Convert the output matrix to a cell array with 8 separate timeseries
    out_timeseries = mat2cell(out_timeseries, ones(1,8), length(digit_filtered));
end


% function [out_timeseries] = digit_decompose(digit_norm, nSamp)
%     % Decompose a time series with sums of 2^0, 2^1, ..., 2^7 into 8 time series
%     % and filter out values that last shorter than nSamp/10.
%     % Inputs:
%     %   digit_norm - input time series (vector)
%     %   nSamp - sampling rate
%     % Outputs:
%     %   out_timeseries - cell array of 8 time series vectors (binary), where each 
%     %                    indicates the presence of each respective power of 2.
%     
%     % Define the minimum duration (in samples) that a value must last
%     min_duration = round(nSamp / 200); % to cutoff short noises that last less than 5ms
%     
%     % Initialize a filtered version of digit_norm
%     digit_filtered = digit_norm;
%     
%     % Find changes in the input time series
%     changes = [1, find(diff(digit_norm) ~= 0) + 1, length(digit_norm) + 1];
%     
%     % Filter out any segment that lasts shorter than min_duration
%     for ii = 1:(length(changes) - 1)
%         segment_start = changes(ii);
%         segment_end = changes(ii+1) - 1;
%         if (segment_end - segment_start + 1) < min_duration
%             % Set the values of this short segment to 0 (treated as noise)
%             digit_filtered(segment_start:segment_end) = 0;
%         end
%     end
%     
%     % Define the 8 powers of 2 (2^0, 2^1, ..., 2^7)
%     powers_of_2 = 2.^(0:7);
%     
%     % Initialize output timeseries as a matrix with 8 rows (one for each power of 2)
%     out_timeseries = zeros(8, length(digit_filtered));
%     
%     % Iterate through each time point in the filtered timeseries
%     for t = 1:length(digit_filtered)
%         % Check if the value is a valid sum of powers of 2
%         value = digit_filtered(t);
%         if all(ismember(dec2bin(value), ['0', '1'])) && value >= 0 && value < 256
%             % Decompose the value using bitwise comparison
%             for ii = 1:8
%                 if bitand(value, powers_of_2(ii)) ~= 0
%                     out_timeseries(ii, t) = 1;
%                 end
%             end
%         else
%             % Treat it as noise, skip decomposition
%             % (output remains 0 by default)
%         end
%     end
%     
%     % Convert the output matrix to a cell array with 8 separate timeseries
%     out_timeseries = mat2cell(out_timeseries, ones(1,8), length(digit_filtered));
% end
