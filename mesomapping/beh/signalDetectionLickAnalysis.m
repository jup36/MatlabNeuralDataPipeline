function drez = signalDetectionLickAnalysis(tbytDat, trRwdI, trPnsI)
%trRwdI = var.stimopts.rewarded_stim; %
%trPnsI = var.stimopts.punished_stim;

% classify trials
rewarded = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.water}, 'un', 0));
punished = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.airpuff}, 'un', 0));

% compute overall hit, miss, fa, cr rates
drez.hitRate = sum(rewarded(:) & trRwdI(:))/sum(trRwdI); % hit rate
drez.missRate = 1-drez.hitRate; % miss rate sum(~rewarded(:) & trRwdI(:))/sum(trRwdI);

if sum(trPnsI)>5 % If there were No-Go trials

    drez.FaRate = sum(punished(:) & trPnsI(:))/sum(trPnsI); % false alarm rate
    drez.CrRate = 1-drez.FaRate; % correct rejection rate sum(~punished(:) & trPnsI(:))/sum(trPnsI);

    % d' (discriminability): The metric d' is a measure of how well an
    % observer can distinguish between two different stimuli or
    % conditions (signal vs. noise).
    drez.dprime = norminv(drez.hitRate) - norminv(drez.FaRate);

end

end