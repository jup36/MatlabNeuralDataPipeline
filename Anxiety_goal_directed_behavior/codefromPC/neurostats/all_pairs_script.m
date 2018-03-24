%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load jitterstuff.mat

% sizes
%[N,T] = size(sp);
tmin = 0;
tstep = .001;
tmax = 20;

% reduce number for testing
N = 5;

% pair mappings
mapping = reshape(1:N*N,N,N);

% get the distances
%D = rand(N,N);
D = D(1:N,1:N);
D(diag(true(N,1))) = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = 250;  % number of bins (+/-) to keep in ccg
cstr = 'rowvalid';  % 'rowvalid' leaves out spikes in train 1 to prevent ccg edge effects
ccalgK = [];  % algorithm control, [] uses default
ccscalgK = [];  % algorithm control, [] uses default
ccssalgK = []; % algorithm control, [] uses default
ccjitalgK = []; % algorithm control, [] uses default

L = 50;  % number of bins in jitter window
Lstart = 1;  % partition start, [] uses default, 1 is faster but less elegant

dplot = [0 0.5 2.5:2.5:12.5]; % partition for distances (used for plotting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert to indices and count spikes
nspike = zeros(N,1);

for j = 1:N
    [S(j,:), M] = spikes2ndx(S(j,:),tmin,tstep,tmax);

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
cc = zeros(2*w-1,N*N);
ccsc = zeros(2*w-1,N*N);
ccss = zeros(2*w-1,N*N);
ccjit = zeros(2*w-1,N*N);

ccgHz = zeros(2*w-1,N*N);
ccHz = zeros(2*w-1,N*N);
ccFrac = zeros(2*w-1,N*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0; tic % used by progress report

% loop over cell 1
for j = 1:N 
    
    % get the spikes for cell 1
    ndx = S(j,:);
    nspikej = nspike(j);
    
    % spike shuffle jitter preprocessing
    ss = ndx2spikeshuffle(ndx,[1 M+1]);
    p = mean4spikeshuffle2ndx2bins(ss,M,L,Lstart);
    
    % jitter preprocessing
    jit = ndx2jitter(ndx,[1 M+1]);
    
    % loop over cell 2
    for k = j+1:N 
        
        % get the index for this cell pair into the ccg data storage
        col = mapping(j,k);
    
        % get the spikes for cell 2
        ndy = S(k,:);
        nspikek = nspike(k);
        
        % compute the ccg
        [cc(:,col),ccdenom] = ndx2jpsth2ccg(ndx,M,ndy,M,'raw',w,cstr,ccalgK);
        
        % record the ccHz and ccFrac rescaling
        gmean(col) = ((nspikej+eps)*(nspikek+eps))^0.5;
        ccdenomkeep(:,col) = ccdenom;
        ccgHz(:,col) = ccdenom .* (T./(gmean(col)));
        ccHz(:,col) = ccdenom .* (1000.*T./(nspikej+eps));
        ccFrac(:,col) = ccdenom .* (1000.*T./(nspikej+eps)./(nspikej./(T.*M./1000)+eps));
        
        % compute the shuffle corrected ccg
        ccsc(:,col) = ndx2jpsth2ccg(ndx,M,ndy,M,'cov',w,cstr,ccalgK);
                
        %% compute the spike shuffle corrected ccg
        ccss(:,col) = cc(:,col) - mean4spikeshuffle2ndx2jpsth2ccg(ss,M,L,Lstart,ndy,M,w,cstr,ccssalgK,p);
        
        % compute the jitter corrected ccg
        ccjit(:,col) = cc(:,col) - mean4jitter2ndx2jpsth2ccg(jit,M,L,Lstart,ndy,M,w,cstr,ccssalgK);
        
        % progress report
        count = count + 1; mytimer(count,N*N,true);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time lags for x-axis plotting
x = 1-w:w-1;

% distance plotting

dplotn = numel(dplot)-1;

% raw ccgs
figure('name','raw ccgs')

for k = 1:dplotn
    
    colndx = mapping(D >= dplot(k) & D < dplot(k+1));
    
    subplot(dplotn,4,(k-1)*4+1)
    plot(x,mean(cc(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+2)
    plot(x,mean(ccsc(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+3)
    plot(x,mean(ccss(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+4)
    plot(x,mean(ccjit(:,colndx),2))    
    
end

% Hz ccgs
figure('name','Hz ccgs')

for k = 1:dplotn
    colndx = mapping(D >= dplot(k) & D < dplot(k+1));

    subplot(dplotn,4,(k-1)*4+1)
    plot(x,mean(cc(:,colndx).*ccgHz(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+2)
    plot(x,mean(ccsc(:,colndx).*ccgHz(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+3)
    plot(x,mean(ccss(:,colndx).*ccgHz(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+4)
    plot(x,mean(ccjit(:,colndx).*ccgHz(:,colndx),2))    
    
end

% frac ccgs
figure('name','fractional Hz ccgs')

for k = 1:dplotn
    
    colndx = mapping(D >= dplot(k) & D < dplot(k+1));
    
    subplot(dplotn,4,(k-1)*4+1)
    plot(x,mean(cc(:,colndx).*ccFrac(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+2)
    plot(x,mean(ccsc(:,colndx).*ccFrac(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+3)
    plot(x,mean(ccss(:,colndx).*ccFrac(:,colndx),2))
    
    subplot(dplotn,4,(k-1)*4+4)
    plot(x,mean(ccjit(:,colndx).*ccFrac(:,colndx),2))    
    
end
