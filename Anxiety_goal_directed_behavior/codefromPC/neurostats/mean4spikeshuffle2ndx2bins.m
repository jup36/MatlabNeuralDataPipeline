function p = mean4spikeshuffle2ndx2bins(ss,M,L,Lstart)
% function p = mean4spikeshuffle2ndx2bins(ss,M,L,Lstart)
%
% p is a MxT matrix 
% p is the expected value of ndx2bins(spikeshuffle2ndx(ss,L,[],Lstart),M),
%  i.e., p(m,t) is the expected number of spikes indices that are equal to m 
%  in the tth trial using a call to spikeshuffle2ndx(ss,L,[],Lstart)
%
% (if Lstart is not provided, then the expected value includes averaging
% over Lstart, which is randomly chosen by spikeshuffle2ndx)
%
% spikes outisde of [1,M] are ignored by p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract data
event = ss{1};
minfix = ss{2};
walls = ss{3};
T = ss{4};

% M check
if numel(M) ~= 1 || M <= 0 || M ~= round(M)
    error('M must be a positive integer scalar')
end

% L check
if numel(L) ~= 1 || L <= 0 || L ~= round(L)
    error('L must be a positive integer scalar')
end
Lminus1 = L-1;

% Lstart check
LstartF = true;
if nargin < 4 || isempty(Lstart)
    LstartF = false;
elseif numel(Lstart) ~= 1 || Lstart <= 0 || Lstart > L || Lstart ~= round(Lstart)
    error('Lstart must be a scalar integer between 1 and L')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE OUT OF BOUNDS SPIKES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isinf(minfix)
    event = event(event(:,1) >= 2-L & event(:,1) <= M+Lminus1,:);
    N = size(event,1);
    Nmax = N;
else
    ndx = event(:,1) >= 2-L & event(:,1) <= M+Lminus1;
    minfix = sum(ndx(1:minfix-1))+1;
    event = event(ndx,:);
    N = size(event,1);
    Nmax = minfix-1;
end

if N == 0, p = zeros(M,T); return, end

% remove out of bounds walls
% shift walls to index into b
% shift walls back 1 for use with cumsum
walls = walls(walls >= 2-L & walls <= M+Lminus1) + Lminus1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% different behavior based on L
if LstartF
    Lst = Lstart;
    Len = Lstart;
    Ln = 1;
else
    Lst = 1;
    Len = L;
    Ln = L;
end

% generate the cumulative bins
% sort for fast indexing into b
ev = sortrows(event(1:Nmax,:),[2 1]);
% shift for indexing into b
ev(:,1) = ev(:,1) + L;
b = zeros(M+2*L-1,T);
% populate
for k = 1:Nmax
    rk = ev(k,1);
    ck = ev(k,2);
    
    b(rk,ck) = b(rk,ck) + 1;
end
% cumsum
b = cumsum(b);

% generate the cumulative psth (normalized)
d = sum(b,2);

% get the output
p = zeros(M,T);

% loop over Lstarts
for Ls = Lst:Len
    
    % get the jitter boundaries
    % shift for indexing into b
    % shift back one for use with cumsum
    bnds = ((floor((1-Ls)./L)*L+Ls+Lminus1):L:(floor((M-Ls)./L)*L+Ls+L+Lminus1)).';
    % add walls to this partition
    bnds = unique([min(max(walls,bnds(1)),bnds(end)); bnds]);
    bndsN = numel(bnds)-1;

    % loop over boundaries
    for k = 1:bndsN
        
        js = bnds(k); % one bin before the starting bin of this piece
        je = bnds(k+1); % ending bin of this piece
        
        % sum over all trials
        % deal with divide by zero without if statement
        dk = (d(je)-d(js)) + eps;
        
        % get the indices into p
        jsp = max(js-Lminus1,1);
        jep = min(je-L,M);
        
        % update p
        p(jsp:jep,:) = p(jsp:jep,:) + ((ones(jep-jsp+1,1)./dk) * (b(je,:)-b(js,:)));
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCALE p BY PSTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = p.*((d(L+1:M+L)-d(L:M+L-1))*(ones(1,T)./Ln));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEAL WITH FIXED SPIKES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if N > Nmax
    
    
    % sort fixed spikes according to trial to speed up moving through p
    efix = sortrows(event(Nmax+1:N,:),[2 1]);
    % remove out of bounds spikes
    efix = efix(efix(:,1) >= 1 & efix(:,1) <= M,:);
    
    % add the spikes to p
    for k = 1:size(efix,1)
        
        rk = efix(k,1);
        ck = efix(k,2);
        
        p(rk,ck) = p(rk,ck) + 1;
        
    end
    
end
