function [tsp, lambda1, pHist] = simIMI(lambda1, lambda2)
%% simuluate neuron with IMI model
% this simulation method is taken from Jonathan Pillow's code
% basic idea is from time-rescaling theorem:
% http://www.stat.columbia.edu/~liam/teaching/neurostat-fall13/papers/brown-et-al/time-rescaling.pdf

Nbin = length(lambda1);         %number of sample bins

if nargin<2
    loglambda2 = 0;
    t_post = 1;
else
    loglambda2 = log(lambda2);
    t_post = length(loglambda2);    %length of post-spike effect
end

dt = 0.001;         %time resolution

maxR = 50;          %estimated maximum firing rate
maxspk = ceil(Nbin*dt*maxR);    %maximum estimated number of spikes within whole interval
tsp = zeros(maxspk, 1);         %cell stucture data for saving spike trains
nbinsPerEval = ceil(1000/maxR); %test if a spike happens after this period

% tspnext = exprnd(1);    %time of next spike (in rescaled time)
% rprev = exprnd(1);      %integrated rescaled time up to current time
% if tspnext<=rprev
%     while tspnext<=rprev
%         rprev = exprnd(1);
%     end
%     nsp = 1;
%     tsp(nsp, 1) = dt;
%     rprev = 0;
%     jbin = 2;       %current time bin
% else
%     nsp = 0;     %number of spikes
%     jbin = 1;       %current time bin
% end
tspnext = exprnd(1);    %time of next spike (in rescaled time)
rprev = 0;      %integrated rescaled time up to current time
% while tspnext<=rprev
%     rprev = exprnd(1);
% end
nsp = 0;     %number of spikes
jbin = 1;       %current time bin

p = log(lambda1*dt);
pHist = 0*lambda1;
while jbin<=Nbin
    iinxt = jbin:min(jbin+nbinsPerEval-1, Nbin);    %indice of bins to be evaluated
    rrnxt = exp(p(iinxt));
    rrcum = cumsum(rrnxt)+rprev;        %compute rescaled time
    
    if tspnext>=rrcum(end)              %no spike happens
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else
        ispk = iinxt(find(rrcum>=tspnext,1)); %time bin where spike occurs
        nsp = nsp+1;
        tsp(nsp,1) = ispk*dt;           %spike times, unit: second
        mxi = min(Nbin, ispk+t_post);
        iiPostSpk = ispk+1:mxi;
        
        if ~isempty(iiPostSpk)      %updating history effect
            p(iiPostSpk) = p(iiPostSpk)+loglambda2(1:(mxi-ispk));
            pHist(iiPostSpk) = pHist(iiPostSpk)+loglambda2(1:(mxi-ispk));
        end
        
        tspnext = exprnd(1);    %draw next spike time
        rprev = 0;              %reset integrated intensity
        jbin = ispk+1;          %move to next bin
        
        muISI = jbin/nsp;       %re-estimate firing rate and modify nbinsPerEval
        nbinsPerEval = round(1.5*muISI);
    end
end

tsp = tsp(1:nsp);














