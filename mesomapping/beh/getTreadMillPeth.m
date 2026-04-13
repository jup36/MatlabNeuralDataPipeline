function tbytDat = getTreadMillPeth(tbytDat, sigIn, tSig, tRel, varargin)
% getTreadMillPeth  Peri-event interpolation with pre-smoothing & zero-centering
%
%   [peth, tRel, inBounds] = getTreadMillPeth(tbytDat, sigIn, tSig, tRel, ...)
%
% Inputs
%   tbytDat : struct array with field .evtOn (event times, seconds)
%   sigIn   : T×1 (or T×D) raw voltage signal(s)
%   tSig    : T×1 timestamps (seconds), monotonic
%   tRel    : 1×M peri-event time grid (e.g., -1:0.02:6)
%
% Name-Value Options (all optional)
%   'Method'   : interp1 method, default 'linear'
%   'Extrap'   : logical, default false (NaN out-of-range)
%   'SgWinSec' : Savitzky–Golay window (seconds), default 0.10
%   'SgOrder'  : Savitzky–Golay polynomial order, default 3
%   'VoltToVel': [a b] linear map: vel = a*voltage + b (default [606.06 -1000])
%   'ZeroCenter':'median'|'mean', default 'median'
%
% Outputs
%   tbytDat     : Now contains 'trSpdLocom' as a field, which is the
%                treadmill velocity. 
%
% J. Park, 2025

% ---------- options
p = inputParser;
p.addParameter('Method','linear',@(s)ischar(s)||isstring(s));
p.addParameter('Extrap',false,@(b)islogical(b)||isscalar(b));
p.addParameter('SgWinSec',0.10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('SgOrder',3,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('VoltToVel',[606.06 -1000],@(v)isnumeric(v)&&numel(v)==2);
p.addParameter('ZeroCenter','median',@(s)ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;

% ---------- prepare inputs
sigIn = double(sigIn);
tSig  = tSig(:);
tRel  = tRel(:).';            % keep as row for N×M
[T, D] = size(sigIn);

if ~isvector(tSig) || numel(tSig) ~= T
    error('tSig must be T×1 and match size(sigIn,1).');
end
if any(diff(tSig) <= 0)
    error('tSig must be monotonically increasing.');
end

evtOn = [tbytDat.evtOn];
N = numel(evtOn);

% ---------- 1) convert voltage → velocity, smooth, zero-center
% fill NaNs gently to avoid filter artifacts
sigFilled = fillmissing(sigIn,'pchip');
sigFilled = fillmissing(sigFilled,'nearest','EndValues','nearest');

% linear conversion (commutes with smoothing, but we’ll do it first)
a = opt.VoltToVel(1); b = opt.VoltToVel(2);
sigVel = a*sigFilled + b;

% SG window in samples (odd, >= order+2 ideally)
fs  = 1/median(diff(tSig),'omitnan');
Lw  = max(5, round(opt.SgWinSec*fs));
if mod(Lw,2)==0, Lw = Lw+1; end
if Lw <= opt.SgOrder, Lw = opt.SgOrder+2 + mod(opt.SgOrder+2,2); end

sigSmooth = zeros(size(sigVel));
sigSmooth = sgolayfilt(sigVel, opt.SgOrder, Lw);

% zero-center (median by default)
switch lower(string(opt.ZeroCenter))
    case "median"
        ctr = median(sigSmooth, 1, 'omitnan');
    case "mean"
        ctr = mean(sigSmooth, 1, 'omitnan');
    otherwise
        error('ZeroCenter must be ''median'' or ''mean''.');
end
sigZero = sigSmooth - ctr;     % center around 0

% ---------- 2) build PETH by chunked interpolation
doExtrap = opt.Extrap;
extrapArg = 'extrap';
if ~doExtrap, extrapArg = NaN; end

tMin = tSig(1); tMax = tSig(end);

for n = 1:N
    tAbs = evtOn(n) + tRel;                
    inBounds(n,:) = (tAbs >= tMin) & (tAbs <= tMax);

    % minimal covering segment in tSig
    tLo = max(min(tAbs), tMin);
    tHi = min(max(tAbs), tMax);
    if tHi < tLo
        continue; % out of range → remains NaN
    end

    i1 = max(1, find(tSig >= tLo, 1, 'first') - 1);
    i2 = min(T, find(tSig <= tHi, 1, 'last')  + 1);
    tSeg = tSig(i1:i2);

    xSeg = sigZero(i1:i2);
    tbytDat(n).trSpdLocom = interp1(tSeg, xSeg, tAbs, opt.Method, extrapArg);
end

end