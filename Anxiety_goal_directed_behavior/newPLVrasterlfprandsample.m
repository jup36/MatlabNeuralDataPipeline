function [ plv, plvRay, pfPhase, spkphasemat, lfpPhase ] = newPLVrasterlfprandsample( lfp, raster, uv, blockorder )
%This function computes the PLVs (Phase-locking values) 
% Updated from the 'PLVrasterlfprandsample' to compute the spkphasemat and lfpPhase, 
% that provide each spike's phase angles as well as the entire LFP signal's
% phase angles. 

phsedge = 0:pi/10:2*pi; % phase edge to be used for histc to determine the preferred phase of a unit
spkphasemat = zeros(size(raster,1),size(raster,2)); % zeromat to contain instantaneous phases for spikes

% Detrend lfp
lfpdt = locdetrend(lfp, uv.Fs, uv.detrend_win); % detrend lfp; signal, sampling frequency, moving window
t = (1:size(lfpdt, 1))'/uv.Fs; % time axis
% Filter lfp
filter_lfp = FilterLFP([t, lfpdt], 'passband', uv.passband);
% Instantaneous phase
lfpPhase = nan(length(uv.edges),size(lfp,2));    % matrix to contain the instantaneous phases
lfpPhase = angle(hilbert(filter_lfp(:, 2:end))); % instantaneous phase of oscillation
lfpPhase(lfpPhase<0) = lfpPhase(lfpPhase<0)+2*pi; % convert the radians to span the 0 to 2 pi range

if strcmp(blockorder,'A')==1     % Ascending order block
    tempidx = [1:size(lfp,2)];
    B1idx = tempidx<=50;           % B1 index
    B2idx = ~B1idx & tempidx<=100; % B2 index
    B3idx = ~B1idx & ~B2idx;       % B3 index
    
elseif strcmp(blockorder,'D')==1 % Descending order block
    tempidx = [1:size(lfp,2)];
    B3idx = tempidx<=50;           % B3 index
    B2idx = ~B3idx & tempidx<=100; % B2 index
    B1idx = ~B2idx & ~B3idx;       % B1 index
end


% B1
spkidx = zeros(size(full(raster)));
spkidx(:,B1idx) = full(raster(:,B1idx));
spkidx = logical(spkidx);   % Spike index for B1


if sum(spkidx(:))>=100 % if there are more than 100 spikes in B1
    B1Phase = lfpPhase(spkidx);
    B1Phase = B1Phase(~isnan(B1Phase));
    spkphasemat(spkidx) = lfpPhase(spkidx);
    B1spkphasemat = spkphasemat(:,B1idx);
    
    for i = 1:1000 % resample 1000 times
        k(:,i) = randsample(length(B1Phase),100);     
    end
    
    % get B1 PLV
    rsB1Phase = B1Phase(k); % randomly resampled B1 phase       
    % rsB1PhaseVec = reshape(rsB1Phase,[100000 1]);
    % rsB1PhaseVec(rsB1PhaseVec<0) = rsB1PhaseVec(rsB1PhaseVec<0)+2*pi; % convert the range of radian to 0 to 2 pi 
    % rsB1phasedist = histc(rsB1PhaseVec,phsedge)./(length(rsB1PhaseVec)/100); % rsB1 phase distribution ranging 0 to 2 pi in Percentage
    % rsB1phasedist = rsB1phasedist(1:end-1,:); % cut the last element which is always zero
    % bar([rsB1phasedist;rsB1phasedist],1) % histogram for phase distribution
    % axis([0.5 40.5 0 15])
    
    plv(1,1) = mean(abs(sum(exp(1i*(rsB1Phase))))./100); % phase-locking value (mean resultant vector)
    [plvRay(1,1),~] = circ_rtest(B1Phase);    % Rayleigh's test for nonuniformity
    
    % get B1 preferred phase
    B1PhaseVec = B1Phase; % B1 phase vector
    [~,B1phsmaxidx] = max(histc(B1PhaseVec, phsedge)); % find the preferred phase
    %bar(histc(B1PhaseVec, phsedge)./(sum(spkidx(:))/100),1) % histogram for phase distribution 
    pfPhase(1,1) = phsedge(B1phsmaxidx+1); % get the Block1 preferred phase
    
else % if there are fewer than 100 spikes in B1
    plv(1,1) = NaN;
    plvRay(1,1) = NaN;
    pfPhase(1,1) = NaN;
    B1spkphasemat = nan(length(raster),sum(B1idx));
end

% B2
spkidx = zeros(size(full(raster)));
spkidx(:,B2idx) = full(raster(:,B2idx));
spkidx = logical(spkidx);   % Spike index for B2
if sum(spkidx(:))>=100 % if there are more than 100 spikes in B2
    B2Phase = lfpPhase(spkidx);
    B2Phase = B2Phase(~isnan(B2Phase));
    spkphasemat(spkidx) = lfpPhase(spkidx);
    B2spkphasemat = spkphasemat(:,B2idx);
    
    for i = 1:1000 % resample 1000 times
        k(:,i) = randsample(length(B2Phase),100);     
    end
    
    rsB2Phase = B2Phase(k); % randomly resampled B2 phase
    plv(1,2) = mean(abs(sum(exp(1i*(rsB2Phase))))./100); % phase-locking value (mean resultant vector)
    [plvRay(1,2),~] = circ_rtest(B2Phase);    % Rayleigh's test for nonuniformity
    
    % get B2 preferred phase
    B2PhaseVec = B2Phase; % B2 phase vector
    [~,B2phsmaxidx] = max(histc(B2PhaseVec, phsedge)); % find the preferred phase
    %bar(histc(B2PhaseVec, phsedge)) % histogram for phase distribution 
    pfPhase(1,2) = phsedge(B2phsmaxidx+1); % get the Block1 preferred phase

else % if there are fewer than 100 spikes in B2
    plv(1,2) = NaN;
    plvRay(1,2) = NaN;
    pfPhase(1,2) = NaN;
    B2spkphasemat = nan(length(raster),sum(B2idx));
end

% B3
spkidx = zeros(size(full(raster)));
spkidx(:,B3idx) = full(raster(:,B3idx));
spkidx = logical(spkidx);   % Spike index for B3
if sum(spkidx(:))>=100 % if there are more than 100 spikes in B3
    B3Phase = lfpPhase(spkidx);
    B3Phase = B3Phase(~isnan(B3Phase));
    spkphasemat(spkidx) = lfpPhase(spkidx);
    B3spkphasemat = spkphasemat(:,B3idx);
    
    for i = 1:1000 % resample 1000 times
        k(:,i) = randsample(length(B3Phase),100);     
    end
    
    rsB3Phase = B3Phase(k); % randomly resampled B3 phase   
    plv(1,3) = mean(abs(sum(exp(1i*(rsB3Phase))))./100); % phase-locking value (mean resultant vector)
    [plvRay(1,3),~] = circ_rtest(B3Phase);    % Rayleigh's test for nonuniformity
    
    % get B3 preferred phase
    B3PhaseVec = B3Phase; % B3 phase vector
    [~,B3phsmaxidx] = max(histc(B3PhaseVec, phsedge)); % find the preferred phase
    %bar(histc(B3PhaseVec, phsedge)) % histogram for phase distribution 
    pfPhase(1,3) = phsedge(B3phsmaxidx+1); % get the Block1 preferred phase

else % if there are fewer than 100 spikes in B3
    plv(1,3) = NaN;
    plvRay(1,3) = NaN;
    pfPhase(1,3) = NaN;
    B3spkphasemat = nan(length(raster),sum(B3idx));
end

spkphasemat = [B1spkphasemat,B2spkphasemat,B3spkphasemat]; % concatenate the spkphasemat

end

