function p = mean4jitter2ndx2bins(jit,M,L,Lstart)
% function p = mean4jitter2ndx2bins(jit,M,L,Lstart)
%
% p is a MxT matrix 
% p is the expect value of ndx2bins(jitter2ndx(jit,L,[],Lstart),M), i.e.,
%  p(m,t) is the expected number of spikes indices that are equal to m 
%  in the tth trial using a call to jitter2ndx(jit,L,[],Lstart)
%
% (if Lstart is not provided, then the expected value includes averaging
% over Lstart, which is randomly chosen by jitter2ndx)
%
% spikes outisde of [1,M] are ignored by p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no cell input
if ~iscell(jit)
    jit = {jit};
end

T = numel(jit);

% M check
if numel(M) ~= 1 || M <= 0 || M ~= round(M)
    error('M must be a positive integer scalar')
end

% L check
if numel(L) ~= 1 || L <= 0 || L ~= round(L)
    error('L must be a positive integer scalar')
end

% Lstart check
LstartF = true;
if nargin < 4 || isempty(Lstart)
    LstartF = false;
elseif numel(Lstart) ~= 1 || Lstart <= 0 || Lstart > L || Lstart ~= round(Lstart)
    error('Lstart must be a scalar integer between 1 and L')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM FOR FIXED LSTART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = zeros(M,T);

if LstartF

    for t = 1:T
    
        jitt = jit{t};
        [nt,tmp3] = size(jitt);
        if nt == 0, continue, end
        if tmp3 ~= 3, error(['jit{' num2str(t) '} must be a matrix with 3 columns']), end
    
        % remove noncontributing spikes
        jitt = jitt(jitt(:,2) >= 2-L & jitt(:,2) <= M+L-1 & jitt(:,1) <= M & jitt(:,3) > 1,:);
        [nt,tmp3] = size(jitt);
        if nt == 0, continue, end
        
        % temporary storage
        pt = zeros(M+4*L-2,1);

        % get the jitter window for this Lstart
        st = floor((jitt(:,2)-Lstart)./L).*L + Lstart;
        a = max(st,jitt(:,1))+(2*L-2);
        b = min(st+L,jitt(:,3))+(2*L-2);

        % get the contribution of each jitter window
        Lt = 1./(b-a);

        % loop over spikes and create the step up and down for each jitter window
        for k = 1:nt

            ndx1 = a(k);
            ndx2 = b(k);
            Ltk = Lt(k);

            pt(ndx1) = pt(ndx1) + Ltk;
            pt(ndx2) = pt(ndx2) - Ltk;

        end % k

        % cumsum converts to jitter distribution
        pt = cumsum(pt);

        % remove the initial and final garbage collection entries
        p(:,t) = pt(2*L-1:M+2*L-2);

    end % t
    
    return
    
end % LstartF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM FOR RANDOM LSTART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ltri = 1./(L.*L);
Ltri2 = -2.*Ltri;

for t = 1:T
    
    jitt = jit{t};
    [nt,tmp3] = size(jitt);
    if nt == 0, continue, end
    if tmp3 ~= 3, error(['jit{' num2str(t) '} must be a matrix with 3 columns']), end

    % remove noncontributing spikes
    jitt = jitt(jitt(:,2) >= 2-L & jitt(:,2) <= M+L-1 & jitt(:,1) <= M & jitt(:,3) > 1,:);
    [nt,tmp3] = size(jitt);
    if nt == 0, continue, end

    %%%%%%%%%%%%%%%%%%%
    % triangle windows
    %%%%%%%%%%%%%%%%%%%
    
    % find indices that allow for full triangle 
    triind = jitt(:,1) <= jitt(:,2)-(L-1) & jitt(:,2)+L <= jitt(:,3); 
    
    % get the triangle update indices
    a = jitt(triind,2)+(L-1); % start (1 <= a <= M+2*L-2)
    ab = a + L; % mid
    b = ab + L; % end (1+2L <= b <= M+4*L-2)
    
    nt = numel(a);
    
    % temporary storage
    pt = zeros(M+4*L-2,1);

    % loop over spikes and create the triangular step up and down for each jitter window
    for k = 1:nt

        ak = a(k);
        abk = ab(k);
        bk = b(k);
        
        pt(ak) = pt(ak) + Ltri;
        pt(abk) = pt(abk) + Ltri2;
        pt(bk) = pt(bk) + Ltri;

    end % k

    % use cumsum to turn this into normal steps up and down for jitter
    % windows
    pt = cumsum(pt);
    
    %%%%%%%%%%%%%%%%%%%
    % edge windows
    %%%%%%%%%%%%%%%%%%%
    
    triind = ~triind;
    j1 = jitt(triind,1)+(2*L-2);
    j2 = jitt(triind,2)+(2*L-2);
    j3 = jitt(triind,3)+(2*L-2);
    
    nt = numel(j1);
    
    % loop over Lstarts
    for Ls = 1:L
        
        % get the jitter window for this Lstart
        st = floor((j2-Ls)./L).*L + Ls;
        a = max(st,j1);
        b = min(st+L,j3);
        
        % get the contribution of each jitter window
        Lt = 1./(L.*(b-a));
        
        % loop over spikes and create the step up and down for each jitter window 
        for k = 1:nt
        
            ndx1 = a(k);
            ndx2 = b(k);
            Ltk = Lt(k);
        
            pt(ndx1) = pt(ndx1) + Ltk;
            pt(ndx2) = pt(ndx2) - Ltk;
        
        end % k
    
    end % Ls

    % final cumsum converts to jitter distribution
    pt = cumsum(pt);

    % remove the initial and final garbage collection entries
    p(:,t) = pt(2*L-1:M+2*L-2);

end % t
