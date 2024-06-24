function SI = calculateXcorrSymmetryIndex(xTimeLag, xcorrMatrix)
    % Ensure the inputs are column vectors
    if isrow(xTimeLag)
        xTimeLag = xTimeLag';
    end

    % Find the center index (zero lag)
    zeroLagIndex = find(xTimeLag == 0);

    % Initialize the output Symmetry Index vector
    numCorrelograms = size(xcorrMatrix, 1);
    SI = zeros(numCorrelograms, 1);

    for i = 1:numCorrelograms
        % Extract the current cross-correlogram
        xcorr = xcorrMatrix(i, :)';

        % Split the cross-correlation into negative and positive lags
        negLags = xcorr(1:zeroLagIndex-1);
        posLags = xcorr(zeroLagIndex+1:end);

        % Reverse the negative lags to align with positive lags
        negLagsReversed = flip(negLags);

        % Ensure both halves are of the same length
        minLength = min(length(negLagsReversed), length(posLags));
        negLagsReversed = negLagsReversed(1:minLength);
        posLags = posLags(1:minLength);

        % Calculate the Symmetry Index (SI) for the current cross-correlogram
        %SI(i) = sum(abs(negLagsReversed - posLags)) / sum(abs(negLagsReversed) + abs(posLags));
        SI(i) = sum(negLagsReversed - posLags) / sum(abs(negLagsReversed) + abs(posLags));
    end
end
