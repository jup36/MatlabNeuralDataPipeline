function tmpBasis = nsSrf(data, flag)
%natural spline basis for self-recovery function 
try 
    maxISI = data.maxISI; 
catch 
    maxISI = 100; 
end
Fs = 1000;  %sampling rate

% knots
if length(flag)>1
    %pre-assigned knots 
    knots = flag;
elseif flag<=0
    %draw knots automatically 
    knots = logspace(0, log10(maxISI), 10);
elseif flag==1
    %manually select knots 
    t0 = data.t0t1(1);
    t1 = data.t0t1(2);
    tsp = data.tsp; 
    srfHist = getHistSrf(tsp, t0, t1, Fs, maxISI);
    figure;
    plot(srfHist.x, srfHist.y);
    [knots, ~] = ginput();
    close; 
elseif flag>1
    %number of knots is pre-assigned. 
    knots = logspace(0, log10(maxISI), flag);
end
knots = unique([reshape(round(knots), 1, []), 1, maxISI]);
knots(knots<1) = [];
knots(knots>maxISI) = [];

nOrder = 4;   %cubic spline 
nBasis = length(knots)+nOrder-2; 
basisObj = create_bspline_basis([1, maxISI], nBasis, nOrder, knots);

t = 1:maxISI; 
X = full(eval_basis(t, basisObj));
tmpBasis.knots = knots;
tmpBasis.basisObj = basisObj;
tmpBasis.t = t;
%tmpBasis.basisX = X;
X = bsxfun(@minus, X, mean(X, 1)); 
tmpBasis.basisX = bsxfun(@times, X, 1./std(X)); 
