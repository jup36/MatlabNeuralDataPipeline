function [corrRez,r2Rez] = trjCorrByTrialType(state, stateEst)
%Computes pearson correlation and r-squared metrics between actual and
% estimated states, each input should be trial-by-trialType cell, whose cell
% element should be variables-by-timeBin. 
% state = s.dat.stateP; % actual state
% stateEst = s.dat.statePCtx; % estimated state
s1 = cell2mat(reshape(state,[],1)')'; % state mat
s2 = cell2mat(reshape(stateEst,[],1)')'; % stateEst mat

r2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
r2Rez.all = r2(s1,s2); % r^2 across all trials

corrRez.all = corr(s1,s2,'Rows','complete'); % corr without smoothing
corrRez.all_sm = corr(smooth2a(s1,3,0),smooth2a(s2,3,0),'Rows','complete'); % corr with smoothing across time bins (smoothing seems to hurt usually)

for tt = 1:size(state,2) % trial types
    s1tt = cell2mat(reshape(state(:,tt),[],1)')';
    s2tt = cell2mat(reshape(stateEst(:,tt),[],1)')';
    tempCorr = corr(s1tt,s2tt,'Rows','complete');
    tempCorr_sm = corr(smooth2a(s1tt,3,0),smooth2a(s2tt,3,0),'Rows','complete');
    tempR2 = r2(s1tt,s2tt);
    
    switch tt
        case 1
            corrRez.lelt = tempCorr;
            corrRez.lelt_sm = tempCorr_sm;
            r2Rez.lelt = tempR2;
        case 2
            corrRez.leht = tempCorr;
            corrRez.leht_sm = tempCorr_sm;
            r2Rez.leht = tempR2;
        case 3
            corrRez.rilt = tempCorr;
            corrRez.rilt_sm = tempCorr_sm;
            r2Rez.rilt = tempR2;
        case 4
            corrRez.riht = tempCorr;
            corrRez.riht_sm = tempCorr_sm;
            r2Rez.riht = tempR2;
    end
end
end