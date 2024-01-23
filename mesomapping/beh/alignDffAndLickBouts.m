function [licksOnTf, lickBoxOnTf, dffsOnTfItp, tF] = alignDffAndLickBouts(lickBouts, dffs, frameT, timeW)
% lickBouts = tbytDat(1).LickChunkRel;
% dffs = dffM1;
% frameT = tbytDat(1).frameTrel;
% timeW = [-5 5];

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
% Interpolate
dffsOnTfItp = interp1(frameTOnTf, dffsOnTf, tF, 'linear', 'extrap');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function loc = locateOnArray(array, target)
        % Convert the array to a column vector for broadcasting
        array = array(:);
        target = target(:);

        % Subtract each target value from the entire array and find the minimum absolute value
        [~, loc] = min(abs(array - target.'), [], 1);
    end


end