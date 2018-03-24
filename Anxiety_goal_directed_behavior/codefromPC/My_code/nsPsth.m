function tmpBasis = nsPsth(data, flag, wd, Fs)
%create natural spline basis for PSTH 
%var_range: time interval 
%flag: methods of choosing knots. 
%0-automatical; 1-manually; n-n knots; vecto-knots to be used
%Fs: sampling rate 
%wd: width of interval between two knots

if ~exist('wd', 'var')
    wd = 0.1;       %default width 
end
if ~exist('Fs', 'var')
    Fs = 1000; 
end
if ~exist('flag', 'var')
    flag = 0; 
end

T = ceil(diff(data.t0t1)*Fs);   %number of bins 

if length(flag)>1
    %pre-assigned knots 
    knots = flag;   
elseif flag==0
    %automatically generate knots 
    knots = linspace(1, T, min(ceil(T/(wd*Fs)), 30));
elseif flag==1
    %manually select knots
    tsp = data.tsp;
    
    %use histogram method to get a smoothed version of psth
    psthHist = getHistPsth(tsp, data.t0t1(1), data.t0t1(2), Fs);
    figure;
    plot(1:T, psthHist.y);
    [knots, ~] = ginput();
    close;
elseif flag>1
    %number of knots is pre-assigned
    knots = linspace(1, T, ceil(flag));
end
knots = unique([1, T, reshape(round(knots), 1, [])]);
knots(knots<1) = [];
knots(knots>T) = [];

nOrder = 4;

%psth
nBasis = length(knots)+nOrder-2;
basisObj = create_bspline_basis([1, T], nBasis, nOrder, knots);

t = 1:T;
X = full(eval_basis(t, basisObj));
tmpBasis.knots = knots;
tmpBasis.basisObj = basisObj;
tmpBasis.t = t;
%tmpBasis.basisX = X;
X = bsxfun(@minus, X, mean(X, 1)); 
tmpBasis.basisX = bsxfun(@times, X, 1./std(X)); 