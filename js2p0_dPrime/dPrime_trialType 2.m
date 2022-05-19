function dPrmRez = dPrime_trialType(spike_time_trial,ttI,twoDconvMat)
%spike_time_trial = spike_rStart_conv; 
%ttI = [tV.leloI(pStartI), tV.lehiI(pStartI), tV.riloI(pStartI), tV.rihiI(pStartI)]; 
%twoDconvMat = [-1 -1 1 1; -1 1 -1 1]; % column(1-4): left-low, left-high, right-low, right-high

nTr = size(ttI,1);  

%% trial-shuffled distribution of d' 
for tt = 1:size(ttI,2)
    for rs = 1:1000 % 1000 resampling
        tmpTr = randsample(nTr,sum(ttI(:,tt))); % random sample trials 
        % random sampled spike_time_trial mat
        tmpMat = spike_time_trial(:,tmpTr); % spike (random sample trials)
        tmpNoMat = spike_time_trial(:,~ismember(1:nTr,tmpTr)); % spike (rest)
        rs_dPrm_1d{tt}(:,rs) = (nanmean(tmpMat,2)-nanmean(tmpNoMat,2))./(nanstd(tmpMat,[],2)+nanstd(tmpNoMat,[],2)); 
        rs_dPrm_2d{tt}(:,:,rs) = twoDconvMat(:,tt)*rs_dPrm_1d{tt}(:,rs)'; % projection onto the 2d space
    end
end

rs_dPrm_2d_X = nansum(cell2mat(reshape(cellfun(@(a) squeeze(a(1,:,:)), rs_dPrm_2d, 'un', 0),1,1,[])),3); % X-dim (left-right)
rs_dPrm_2d_Y = nansum(cell2mat(reshape(cellfun(@(a) squeeze(a(2,:,:)), rs_dPrm_2d, 'un', 0),1,1,[])),3); % Y-dim (low-high)

rs_dPrm_2d_XY = [reshape(rs_dPrm_2d_X,[],1),reshape(rs_dPrm_2d_Y,[],1)]; 

% resampled/shuffled bivariate distribution mean, cov to get the pdf 
dPrmRez.mu_shuff_2d_XY = nanmean(rs_dPrm_2d_XY,1); 
dPrmRez.cov_shuff_2d_XY = cov(rs_dPrm_2d_XY); 

% cdf
[Xs,Ys] = meshgrid(linspace(-3,3,600)',linspace(-3,3,600)'); 
XYs = [Xs(:),Ys(:)]; 
p_cdf0 = mvncdf(XYs, dPrmRez.mu_shuff_2d_XY, dPrmRez.cov_shuff_2d_XY); 
%p_cdf = reshape(p_cdf0,600,600); 





p_cdf_img = nan(600, 600); 

for j = 1:length(p_cdf0) 
    xy = XYs(j,:); 
    p_cdf_img(Xs==XYs(j,1) & Ys==XYs(j,2)) = p_cdf0(j);  
end







%% actual d' and p-values based on the cdf of the suffled distribution
for tt = 1:size(ttI,2) % trial Types (e.g. lelo, lehi, rilo, rihi)
    tmpMat = spike_time_trial(:,ttI(:,tt)); 
    tmpNoMat = spike_time_trial(:,~ttI(:,tt)); 
    dPrm_1d{tt} = (nanmean(tmpMat,2)-nanmean(tmpNoMat,2))./(nanstd(tmpMat,[],2)+nanstd(tmpNoMat,[],2)+eps); 
    dPrm_2d{tt} = twoDconvMat(:,tt)*dPrm_1d{tt}'; 
end

dPrmRez.trj2d = sum(cell2mat(reshape(dPrm_2d,1,1,[])),3)'; % sum across all trial-type scores

dist2d = @(a) sqrt(a(:,1).^2+a(:,2).^2); 
dPrmRez.trjDistZero = dist2d(dPrmRez.trj2d); % distance from zero

for j = 1:size(dPrmRez.trj2d,1)
    [~,minDistI] = min(dist2d(dPrmRez.trj2d(j,:)-XYs)); % locate each point of the 2d trajectory within the grid
    dPrmRez.dist2d_cdf_p(j,1) = min(p_cdf0(minDistI),1-p_cdf0(minDistI));    
end

%dPrm_2d_distZero(dist2d_cdf_p>0.01,1)=NaN;

%[~,dPrm_sig_maxDistI] = max(dPrm_2d_distZero); 

end