function drez = signalDetectionLickAnalysis(tbytDat, trRwdI, trPnsI)
%trRwdI = var.stimopts.rewarded_stim; %
%trPnsI = var.stimopts.punished_stim;

% classify trials
rewarded = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.water}, 'un', 0));
punished = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.airpuff}, 'un', 0));

% compute overall hit, miss, fa, cr rates
rawHitRate = sum(rewarded(:) & trRwdI(:))/sum(trRwdI); % raw hit rate
drez.hitRate = adjustRate(rawHitRate, sum(trRwdI));

drez.missRate = 1 - drez.hitRate; % miss rate

if sum(trPnsI) > 5 % If there were No-Go trials
    rawFaRate = sum(punished(:) & trPnsI(:))/sum(trPnsI); % raw false alarm rate
    drez.FaRate = adjustRate(rawFaRate, sum(trPnsI));

    drez.CrRate = 1 - drez.FaRate; % correct rejection rate

    % d' (discriminability): The metric d' is a measure of how well an
    % observer can distinguish between two different stimuli or
    % conditions (signal vs. noise).
    drez.dprime = norminv(drez.hitRate) - norminv(drez.FaRate);
end

    function adjustedRate = adjustRate(rawRate, N)
        % Adjusts hit or false alarm rate if it's 1 or 0
        if rawRate == 1
            adjustedRate = (N-0.5) / N;
        elseif rawRate == 0
            adjustedRate = 0.5 / N;
        else
            adjustedRate = rawRate;
        end
    end

end