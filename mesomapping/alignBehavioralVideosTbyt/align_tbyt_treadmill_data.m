function tbytDat = align_tbyt_treadmill_data(filePath, fileKeyword)
% align_tbyt_treadmill_data  Align treadmill signal to trial events and append to tbytDat.
%
%   tbytDat = align_tbyt_treadmill_data(filePath, fileKeyword)
%
% Description
%   Loads the trial-by-trial structure (tbytDat) for a session and the
%   treadmill signal from evtInS.mat, builds a peri-event time base
%   (-1:0.02:6 s relative to cue), and calls getTreadMillPeth to interpolate
%   the treadmill signal around each trial’s event time (tbytDat(evtOn)).
%   The aligned trace for each trial is stored in:
%       tbytDat(tr).trSpdLocom = [aligned_signal; cueAlignedTs]
%
%   Finally, the updated tbytDat is saved back into the same Matfiles
%   location that contains the refit CNMF outputs.
%
% Inputs
%   filePath    : Full path to the session's 'task' folder
%                 e.g., 'Z:\Rodent Data\...\mXXXX_XXXXXX\task'
%   fileKeyword : Filename (or distinctive suffix) that locates the refit CNMF
%                 data file under 'Matfiles' (e.g., '_refitChunks_red_dff_combined.mat'
%                 or '_refitChunks_green_dff_combined.mat'). This file must contain tbytDat.
%
% Output
%   tbytDat     : Trial-by-trial struct with added field per trial:
%                 • trSpdLocom (2×N) → [aligned treadmill signal; cueAlignedTs]
%
% Data Expectations
%   • evtInS.mat exists in the session folder (or subfolders) and contains:
%       evtInS.treadMill  — T×2 matrix: [:,1] = raw treadmill signal
%                                            [:,2] = timestamps (seconds)
%   • The refit CNMF Matfile (found via fileKeyword) contains tbytDat and trial times
%     (tbytDat(tr).evtOn) in the same time base as evtInS.treadMill(:,2).
%
% Notes
%   • Peri-event time base is fixed here to cueAlignedTs = -1:0.02:6 (seconds).
%   • getTreadMillPeth performs the temporal chunking and interpolation; it should
%     populate each trial’s trSpdLocom. If your getTreadMillPeth returns a PETH
%     array instead of updating tbytDat in-place, adapt this caller accordingly.
%   • The function saves the updated tbytDat back into the Matfile located by fileKeyword.
%
% Author: Junchol Park (Buschman Lab)
% Date  : September 2025

%% whereabouts
%filePath = compatiblepath('Z:\Rodent Data\dualImaging_parkj\m1045_jRGECO_GRABda\m1045_121324\task');
%fileKeyword = '_refitChunks_red_dff_combined.mat';
filePath_mat = find_keyword_containing_folder(filePath, 'Matfiles', 'recursive', false);
filePath_dat = find_keyword_containing_files(filePath_mat{1}, fileKeyword, 'recursive', false);
if ~isempty(filePath_dat)
    filePath_evt = findFileInSubfolders(filePath, 'evtInS.mat');
    if ~isempty(filePath_evt)

        % load behavioral data
        if ~isempty(filePath_dat), load(filePath_dat{1}, 'tbytDat'); end
        if ~isempty(filePath_evt)
            S = load(filePath_evt, 'evtInS');
            treadMill = S.evtInS.treadMill;
        end

        cueAlignedTs = -1:0.02:6;
        %trBIdC = {tbytDat.limeLEDTrainI};

        tbytDat = getTreadMillPeth(tbytDat, treadMill(:, 1), treadMill(:, 2), cueAlignedTs);

        %% ---------- Save ----------
        save(fullfile(filePath_dat{1}), 'tbytDat');
    end
end

%% ===== Local helpers =====
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
end