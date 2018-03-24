function [pairlist,ccr,cc,ccs,ccHz] = jpsth_psthjitter(S,tmax,jbins,winsize)
% [pairlist,ccr,cc,ccs,ccHz] = jpsth_psthjitter(S,tmax,jbins,winsize)
%
% Uses the PSTH jitter method to calculate a JPSTH
%
% S is a cell array of cells x trials (spike times in milliseconds)
% tmax = maximum spike time (stimulus duration in milliseconds)
% jbins = number of bins in the jitter window (default is 50 ms)
% winsize = size of CCG (+/- this value, default is 250 ms)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sizes
%[N,T] = size(sp);
tmin = 0;
%tstep = .001; % size of bins (if spike times are in seconds)
tstep = 1; % size of bins (if spike times are in milliseconds)

% tmax should have the same units as spike times
%tmax = 1280; % stimulus duration (max spike time)
%tmax = 1.28; % stimulus duration (max spike time)

% Number of cells and number of trials
[N,T] = size(S);

% pair mappings
%mapping = reshape(1:N*N,N,N);
mapping = zeros(N,N);
count = 1;
for I=1:N
    for J=I+1:N
        mapping(I,J) = count;
        count = count + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cstr = {};
%cstr = 'rowvalid';  % 'rowvalid' leaves out spikes in train 1 to prevent ccg edge effects
ccalgK = [];  % algorithm control, [] uses default
ccscalgK = [];  % algorithm control, [] uses default
ccralgK = []; % algorithm control, [] uses default
ccjitalgK = []; % algorithm control, [] uses default

if (nargin < 4)
    w = 250;  % number of bins (+/-) to keep in ccg
else
    w = winsize;
    if (w > tmax)  
        error('CCG window must be no larger than tmax');
    end
end

if (nargin < 3)
    L = 50;  % number of bins in jitter window
else
    L = jbins;
    if (L > w)
        error('Jitter bin size must be less than ccg bin size');
    end
end
Lstart = 1;  % partition start, [] uses default, 1 is faster but less elegant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert to indices and count spikes
nspike = zeros(N,1);

for j = 1:N
    [S(j,:), M] = spikes2ndx(S(j,:),tmin+tstep/2,tstep,tmax+tstep/2);

    % count the contributing spikes in the jth cell
    % leave out spikes that are omitted by 'rowvalid'
    if strcmp(cstr,'rowvalid')
        rmin = w;
        rmax = M-w+1;
    else
        rmin = 1;
        rmax = M;
    end
    nspikej = 0;
    for t = 1:T
        nspikej = nspikej + sum(S{j,t} >= rmin & S{j,t} <= rmax);
    end
    nspike(j) = nspikej;
end

% initialize the ccg storage
cc = zeros(2*w-1,nchoosek(N,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0; tic % used by progress report

% loop over cell 1
for j = 1:N
    
    % get the spikes for cell 1
    ndx = S(j,:);
    nspikej = nspike(j);
    
    % PSTH jitter preprocessing
    ss = ndx2spikeshuffle(ndx,[1 M+1]);
    p = mean4spikeshuffle2ndx2bins(ss,M,L,Lstart);

    % loop over cell 2
    for k = j+1:N

        % get the index for this cell pair into the ccg data storage
        col = mapping(j,k);
    
        % get the spikes for cell 2
        ndy = S(k,:);
        nspikek = nspike(k);

        % PSTH jitter preprocessing
        ss2 = ndx2spikeshuffle(ndy,[1 M+1]);
        p2 = mean4spikeshuffle2ndx2bins(ss2,M,L,Lstart);

        % store the pair indices for this CCG
        pairlist(:,col) = [j,k];
        
        % jpsth
        jp = bins2jpsth(p,p2,'raw');
        imagesc(jp); pause;
        
        % progress report
        count = count + 1; mytimer(count,nchoosek(N,2),true);
    end
end
