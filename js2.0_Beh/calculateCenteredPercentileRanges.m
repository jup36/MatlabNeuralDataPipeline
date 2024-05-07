function ranges = calculateCenteredSignedPercentileRanges(data, step)
    % Calculate symmetric percentile ranges centered at zero, extending outward, with signs.
    %
    % Parameters:
    %   data - An array containing the data points.
    %   step - The step size for the percentiles in terms of percentage (e.g., 10 for 10%, 20% ...).
    %
    % Returns:
    %   ranges - A matrix where each row represents the [lower_bound, upper_bound] of the percentile range.
    %            Ranges are given with explicit signs to indicate their position relative to zero.
    
    if nargin < 2
        step = 10; % Default step size is 10 if not specified
    end

    % Separate positive and negative data
    positiveData = data(data >= 0);
    negativeData = data(data < 0);

    % Calculate the percentiles for positive and negative data
    positivePercentiles = prctile(positiveData, 0:step:100);
    negativePercentiles = prctile(negativeData, 100:-step:0);

    % Number of percentile ranges to calculate
    numRanges = length(positivePercentiles) - 1;
    ranges = zeros(numRanges * 2, 2);

    % Fill in the negative ranges and flip them
    for i = 1:numRanges
        ranges(i, 1) = -negativePercentiles(numRanges - i + 2);
        ranges(i, 2) = -negativePercentiles(numRanges - i + 1);
    end

    % Fill in the positive ranges
    for i = 1:numRanges
        ranges(numRanges + i, 1) = positivePercentiles(i);
        ranges(numRanges + i, 2) = positivePercentiles(i + 1);
    end

    % Sort the ranges to have a continuous sequence from negative to positive
    ranges = sortrows(ranges);
end
