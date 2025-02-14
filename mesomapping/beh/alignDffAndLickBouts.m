function [licksOnTf, lickBoxOnTf, dffsOnTfItp, tF] = alignDffAndLickBouts(lickBouts, dffs, frameT, timeW)
% lickBouts = tbytDat(t).LickBoutRel;
% dffs = dffM1;
% frameT = tbytDat(t).frameTrel;
% timeW = [-0.9 5];

% get licks on the specified time frame
tF = timeW(1):0.001:timeW(end); % time frame
licksOnTf = zeros(1, length(tF));
lickBoxOnTf = zeros(1, length(tF));

for ts = 1:length(lickBouts)
    lb = lickBouts{ts};

    vlicks = lb(cell2mat(arrayfun(@(a) a>=min(tF) & a<=max(tF), lb, 'UniformOutput', false)));

    if ~isempty(vlicks)
        licksOnTfI = locateOnArray(tF, vlicks);
        licksOnTf(licksOnTfI) = 1;
        lickBoxOnTf(min(licksOnTfI):max(licksOnTfI)) = 1;
    end
end

% get dff on the specified time frame with interpolation
ftI = frameT >= min(tF) & frameT <= max(tF);
dffsOnTf = dffs(ftI);
frameTOnTf = frameT(ftI);

% Define the valid range based on frameTOnTf
validIdx = tF >= min(frameTOnTf) & tF <= max(frameTOnTf);

% Curtail tF to only valid indices
tF = tF(validIdx);
licksOnTf = licksOnTf(validIdx); 
lickBoxOnTf = lickBoxOnTf(validIdx); 

% Perform strict interpolation (no extrapolation)
dffsOnTfItp = interp1(frameTOnTf, dffsOnTf, tF, 'linear');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function loc = locateOnArray(array, target)
        % Convert the array to a column vector for broadcasting
        array = array(:);
        target = target(:);

        % Subtract each target value from the entire array and find the minimum absolute value
        [~, loc] = min(abs(array - target.'), [], 1);
    end


end