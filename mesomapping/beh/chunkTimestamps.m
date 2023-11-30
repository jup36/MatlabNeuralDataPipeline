function chunks = chunkTimestamps(timestamps, cutoff)
% chunk timestamps whose intervals are shorter than the cutoff
    timestamps = timestamps(:); 
    
    % Calculate the difference between consecutive timestamps
    intervals = [Inf; diff(timestamps)];

    % Initialize variables
    chunks = {};
    startIdx = 1; % Start index of the current chunk

    % Iterate through the intervals
    for i = 2:length(timestamps)
        if intervals(i) > cutoff
            % If interval exceeds cutoff, end current chunk and start a new one
            chunks{end + 1} = timestamps(startIdx:i-1);
            startIdx = i;
        end
    end
    
    % Add the last chunk
    chunks{end + 1} = timestamps(startIdx:end);
end