function dPrmRez = dPrime_trialType_vector_norm_stat(spike_time_trial,ttI,twoDconvMat)
%spike_time_trial = spike_rStart_conv; 
%ttI = [tV.leloI(pStartI), tV.lehiI(pStartI), tV.riloI(pStartI), tV.rihiI(pStartI)]; 
%twoDconvMat = [-1 -1 1 1; -1 1 -1 1]; % column(1-4): left-low, left-high, right-low, right-high

nTr = size(ttI,1);  
nTr_sample = round(nTr/4); 

distf = @(a, b) sqrt(a^2 + b^2); 

%% trial-shuffled distribution of d' 
for tt = 1:size(ttI,2)
    for rs = 1:1000 % 1000 resampling
        tmpTr = randsample(nTr, nTr_sample); % randsample(nTr,sum(ttI(:,tt))); % random sample trials 
        % random sampled spike_time_trial mat
        tmpMat = spike_time_trial(:,tmpTr); % spike (random sample trials)
        tmpNoMat = spike_time_trial(:,~ismember(1:nTr,tmpTr)); % spike (rest)
        rs_dPrm_1d{tt}(:,rs) = (nanmean(tmpMat,2)-nanmean(tmpNoMat,2))./(nanstd(tmpMat,[],2)+nanstd(tmpNoMat,[],2)); 
        rs_dPrm_2d{tt}(:,:,rs) = twoDconvMat(:,tt)*rs_dPrm_1d{tt}(:,rs)'; % projection onto the 2d space
    end
end

rs_dPrm_2d_X = nansum(cell2mat(reshape(cellfun(@(a) squeeze(a(1,:,:)), rs_dPrm_2d, 'un', 0),1,1,[])),3); % X-dim (left-right)
rs_dPrm_2d_Y = nansum(cell2mat(reshape(cellfun(@(a) squeeze(a(2,:,:)), rs_dPrm_2d, 'un', 0),1,1,[])),3); % Y-dim (low-high)

max_dist = nan(1000, 1);  
for j = 1:1000 % 1000 resampling
    max_dist(j) = max(cell2mat(arrayfun(@(a, b) distf(a, b), rs_dPrm_2d_X(:,j), rs_dPrm_2d_Y(:,j), 'un', 0))); 
end

%% actual d' and p-values 
for tt = 1:size(ttI,2) % trial Types (e.g. lelo, lehi, rilo, rihi)
    tmpMat = spike_time_trial(:,ttI(:,tt)); 
    tmpNoMat = spike_time_trial(:,~ttI(:,tt)); 
    dPrm_1d{tt} = (nanmean(tmpMat,2)-nanmean(tmpNoMat,2))./(nanstd(tmpMat,[],2)+nanstd(tmpNoMat,[],2)+eps); 
    dPrm_2d{tt} = twoDconvMat(:,tt)*dPrm_1d{tt}'; 
end

dPrmRez.trj2d = sum(cell2mat(reshape(dPrm_2d,1,1,[])),3)'; % sum across all trial-type scores

dist2d = @(a) sqrt(a(:,1).^2+a(:,2).^2); 
dPrmRez.trjDistZero = dist2d(dPrmRez.trj2d); % distance from zero
dPrmRez.trjDistZero_sigI = dPrmRez.trjDistZero > mean(max_dist) + 2*std(max_dist); 

end