function tbytDat = align_tbyt_pupil_orofacial_data(filePath, fileKeyword)
% align_tbyt_pupil_orofacial_data  Aligns behavioral video features to trial data.
%
%   align_tbyt_pupil_orofacial_data(filePath)
%
%   This function loads, preprocesses, and temporally aligns behavioral
%   features extracted from orofacial, pupil, and lick video recordings
%   to trial events in a behavioral task, and appends the aligned signals
%   into the trial-by-trial structure (tbytDat).
%
%   INPUT:
%     filePath   – Full path to the 'task' folder for a session
%                  (e.g., 'Z:\Rodent Data\...\mXXXX_XXXXXX\task')
%
%   OUTPUT:
%     Updates the in-memory tbytDat struct by adding the following fields
%     to each trial:
%       • trNoseTipEnv  – Envelope of nose-tip velocity (Hilbert)
%       • trWhiskerEnv  – Envelope of whisker velocity
%       • trPupilDia    – Filtered pupil area trace (low-pass, 10 Hz)
%       • trLickTrj     – 1D tongue trajectory
%     Each field is stored as a 2×N matrix: [signal; cueAlignedTs]
%
%   FUNCTIONAL OVERVIEW:
%     • Automatically detects analyzed .mat files for:
%         - Orofacial features (nose tip, whisker)
%         - Pupil area (ellipse)
%         - Lick trajectories
%     • Loads behavioral event timestamps (evtInS.faceCam) and trial struct (tbytDat)
%     • Preprocesses signals with:
%         - Interpolation for short length mismatches
%         - Velocity calculation (centered difference)
%         - Envelope extraction via Hilbert transform (for motion features)
%         - Low-pass filtering for pupil area
%     • Aligns each signal to cue-aligned timestamps (-1 to +6 s)
%     • Robust to missing files or incomplete data
%
%   INTERNAL HELPERS:
%     localMatchLength(sigIn, tgtTime, maxDiff)
%       → Harmonizes signal length to match target timestamps
%
%     localVelocity(sig, t)
%       → Computes velocity and returns midpoint timestamps
%
%     localEnvelope(sig)
%       → Computes analytic signal envelope using Hilbert transform
%
%   NOTES:
%     - Requires that the task folder contains Matfiles, evtInS.mat, and
%       analyzed behavioral video .mat files.
%     - Designed to tolerate small differences in data length and to run
%       independently on each signal type.
%     - Prints progress for each successfully processed block.
%
%   Author: Junchol Park (Buschman Lab)
%   Date:   August 2025

%% whereabouts
%filePath = compatiblepath('Z:\Rodent Data\dualImaging_parkj\m1045_jRGECO_GRABda\m1045_121324\task');
header = extract_date_animalID_header(filePath);
%mId = cell2mat(regexp(header, 'm\d{4}', 'match'));
%fileKeyword = '_refitChunks_red_dff_combined.mat';
filePath_mat = find_keyword_containing_folder(filePath, 'Matfiles', 'recursive', false);
filePath_dat = find_keyword_containing_files(filePath_mat{1}, fileKeyword, 'recursive', false);
if ~isempty(filePath_dat)
    filePath_evt = findFileInSubfolders(filePath, 'evtInS.mat');

    filePathC_orofacial = GrabFiles_sort_trials([header,'*','_vid_cropped_orofacial'], 0, {filePath});
    filePathC_pupil     = GrabFiles_sort_trials([header,'*','_vid_cropped_pupil'],     0, {filePath});

    if ~isempty(filePathC_pupil) && ~isempty(filePathC_orofacial)

        % keep analyzed MATs as cell arrays (no cell2mat!)
        filePathC_orofacialAnalyzed = GrabFiles_sort_trials('behvid*_orofacialAnalyzed', 0, filePathC_orofacial);
        filePathC_lickAnalyzed      = GrabFiles_sort_trials('behvid*_lickAnalyzed',      0, filePathC_orofacial);
        filePathC_pupilAnalyzed     = GrabFiles_sort_trials('behvid*_pupilAnalyzed',     0, filePathC_pupil);

        getBlockNum = @(fname) str2double(regexp(fname, '_m\d{4}_(\d{1,2})_', 'tokens', 'once'));

        % load behavioral data
        if ~isempty(filePath_dat), load(filePath_dat{1}, 'tbytDat'); end
        if ~isempty(filePath_evt)
            S = load(filePath_evt, 'evtInS');
            faceCamT = S.evtInS.faceCam;            % [time, blockNum]
        end

        cueAlignedTs = -1:0.02:6;
        trBIdC = {tbytDat.limeLEDTrainI};

        %% ---------- Process OROFACIAL (nose/whisker) ----------
        for i = 1:numel(filePathC_orofacialAnalyzed)
            fn = filePathC_orofacialAnalyzed{i};
            if isempty(fn) || ~exist(fn,'file'), continue; end

            [~, fileName] = fileparts(fn);
            blockNum = getBlockNum(fileName);

            load(fn, "ntip1d_smooth", "whisker1d_smooth");

            % camera time for this block
            currFaceCamT = faceCamT(faceCamT(:,2)==blockNum,1);
            if isempty(currFaceCamT), warning('No faceCamT for block %g', blockNum); continue; end

            % length harmonization
            ntip1d_smooth    = localMatchLength(ntip1d_smooth,    currFaceCamT, 20);
            whisker1d_smooth = localMatchLength(whisker1d_smooth, currFaceCamT, 20);

            % If either failed, skip this block cleanly
            if isempty(ntip1d_smooth) || isempty(whisker1d_smooth)
                warning('Skipping block %d due to length mismatch.', blockNum);
                continue;
            end

            % velocities + envelopes
            [ntipVel,   tVel] = localVelocity(ntip1d_smooth,    currFaceCamT);
            [whiskVel, ~    ] = localVelocity(whisker1d_smooth, currFaceCamT);

            ntipVelEnv   = localEnvelope(ntipVel);
            whiskerVelEnv= localEnvelope(whiskVel);

            % align to trials
            for t = 1:numel(tbytDat)
                if trBIdC{t} ~= blockNum, continue; end
                tempRelT  = tVel - tbytDat(t).evtOn;
                sel       = (cueAlignedTs(1) <= tempRelT) & (tempRelT <= cueAlignedTs(end));

                tbytDat(t).trNoseTipEnv = [interp1(tempRelT(sel), ntipVelEnv(sel),    cueAlignedTs, 'linear','extrap'); cueAlignedTs];
                tbytDat(t).trWhiskerEnv = [interp1(tempRelT(sel), whiskerVelEnv(sel), cueAlignedTs, 'linear','extrap'); cueAlignedTs];
            end
            fprintf("Processed %s orofacial data of block#%d\n", header, blockNum)
        end

        %% ---------- Process PUPIL ----------
        for i = 1:numel(filePathC_pupilAnalyzed)
            fn = filePathC_pupilAnalyzed{i};
            if isempty(fn) || ~exist(fn,'file'), continue; end

            [~, fileName] = fileparts(fn);
            blockNum = getBlockNum(fileName);

            S = load(fn, "ellipse_area");
            ellipse_area = S.ellipse_area;

            currFaceCamT = faceCamT(faceCamT(:,2)==blockNum,1);
            if isempty(currFaceCamT), warning('No faceCamT for block %g', blockNum); continue; end

            % length harmonization
            ellipse_area = localMatchLength(ellipse_area, currFaceCamT, 20);
            % If either failed, skip this block cleanly
            if isempty(ellipse_area)
                warning('Skipping block %d due to length mismatch.', blockNum);
                continue;
            end

            % filter + envelope (low-pass, then no hilbert for pupil unless needed)
            fs = 200; fc = 10; order = 4; % 10-Hz cutoff low-pass filter
            [b,a] = butter(order, fc/(fs/2), 'low');

            x = double(ellipse_area); x(~isfinite(x)) = NaN;
            x = fillmissing(x,'pchip');
            x = fillmissing(x,'nearest','EndValues','nearest');
            ellipse_area_filt = filtfilt(b,a,x);

            % align to trials
            for t = 1:numel(tbytDat)
                if trBIdC{t} ~= blockNum, continue; end
                tempRelT  = currFaceCamT - tbytDat(t).evtOn;
                sel       = (cueAlignedTs(1) <= tempRelT) & (tempRelT <= cueAlignedTs(end));
                tbytDat(t).trPupilDia = [interp1(tempRelT(sel), ellipse_area_filt(sel), cueAlignedTs, 'linear','extrap'); cueAlignedTs];
            end
            fprintf("Processed %s pupil data of block#%d\n", header, blockNum)
        end

        %% ---------- Process LICK ----------
        if ~isempty(filePathC_lickAnalyzed)
            for i = 1:numel(filePathC_lickAnalyzed)
                fn = filePathC_lickAnalyzed{i};
                if isempty(fn) || ~exist(fn,'file'), continue; end

                [~, fileName] = fileparts(fn);
                blockNum = getBlockNum(fileName);

                S = load(fn, "tongue1d_smooth");
                tongue1d_smooth = S.tongue1d_smooth;

                currFaceCamT = faceCamT(faceCamT(:,2)==blockNum,1);
                if isempty(currFaceCamT), warning('No faceCamT for block %g', blockNum); continue; end

                % length harmonization
                tongue1d_smooth = localMatchLength(tongue1d_smooth, currFaceCamT, 20);
                if isempty(tongue1d_smooth)
                    warning('Skipping block %d due to length mismatch.', blockNum);
                    continue;
                end

                % align to trials
                for t = 1:numel(tbytDat)
                    if trBIdC{t} ~= blockNum, continue; end
                    tempRelT  = currFaceCamT - tbytDat(t).evtOn;
                    sel       = (cueAlignedTs(1) <= tempRelT) & (tempRelT <= cueAlignedTs(end));
                    tbytDat(t).trLickTrj = [interp1(tempRelT(sel), tongue1d_smooth(sel), cueAlignedTs, 'linear','extrap'); cueAlignedTs];
                end
                fprintf("Processed %s lick data of block#%d\n", header, blockNum)
            end

            trLickTrjConcat = cell2mat({tbytDat.trLickTrj});
            [trLickTrjMean, trLickTrjStd] = meanstdsem(trLickTrjConcat(1, :)'); % for video lick thresholding
            lickTrjThres = trLickTrjMean+2*trLickTrjStd;

            % ---------- Count video-based LICKs ----------
            for tr = 1:numel(tbytDat)
                if ~isempty(tbytDat(tr).trLickTrj)
                    aboveThresI = tbytDat(tr).trLickTrj(1, :)>lickTrjThres;
                    thresCrossTs = tbytDat(tr).trLickTrj(2,aboveThresI);
                    if ~isempty(thresCrossTs)
                        thresCrossTsVal = thresCrossTs([true diff(thresCrossTs)>0.1]);
                        tbytDat(tr).periToneLicksVid = thresCrossTsVal(thresCrossTsVal<2.2); % 2.2s relative to tone-onset
                        tbytDat(tr).postToneLicksVid = thresCrossTsVal(thresCrossTsVal>2.2);
                    end
                end
            end
        end

        %% ---------- Save ----------
        saveName = [header, '_tbytDat_alignedPupilOrofacial.mat'];
        if iscell(filePath_mat)
            saveDir = filePath_mat{1};
        else
            saveDir = filePath_mat;
        end
        save(fullfile(saveDir, saveName), 'tbytDat');
    end
end

%% ===== Local helpers =====
    function sigOut = localMatchLength(sigIn, tgtTime, maxDiff)
        % Harmonize the length of sigIn to match tgtTime by interpolation or smart extension.
        % Returns [] (and warns) if the mismatch is too large so caller can skip gracefully.

        Nsig = numel(sigIn);
        Ntgt = numel(tgtTime);
        sigIn = sigIn(:);
        tgtTime = tgtTime(:);

        if Nsig == Ntgt
            sigOut = sigIn;
            return
        end

        if Nsig <= 1
            warning('localMatchLength:TooShort','Signal too short to resample.');
            sigOut = [];
            return
        end

        % Acceptable mismatch if within maxDiff or < 20% of sigIn
        diffLen = abs(Nsig - Ntgt);
        allowExtension = Ntgt < Nsig && (diffLen <= max(maxDiff, 0.2 * Nsig));

        if diffLen > maxDiff && ~allowExtension
            warning('localMatchLength:TooDifferent', ...
                'Length mismatch too large (sig=%d, tgt=%d). Skipping.', Nsig, Ntgt);
            sigOut = [];  % <-- graceful failure
            return
        end

        % Extend tgtTime if too short and allowed
        if allowExtension
            dt = mode(round(diff(tgtTime)*1e3)/1e3);  % round mode to ms precision
            extraTs = tgtTime(end) + dt*(1:(Nsig - Ntgt))';
            tgtTime = [tgtTime; extraTs];
        end

        % Interpolate signal to tgtTime
        tSrc = linspace(tgtTime(1), tgtTime(end), Nsig);
        sigOut = interp1(tSrc, sigIn, tgtTime, 'linear', 'extrap');
    end


    function [v, tMid] = localVelocity(sig, t)
        % Centered-difference velocity and mid-point timestamps
        sig = sig(:); t = t(:);
        v    = diff(sig);
        tMid = t(1:end-1) + 0.5*median(diff(t), 'omitnan');
        % Make sure it's finite for later envelope
        v(~isfinite(v)) = NaN;
        v = fillmissing(v,'pchip');
        v = fillmissing(v,'nearest','EndValues','nearest');
    end

    function env = localEnvelope(sig)
        % Envelope via Hilbert; requires finite inputs
        sig(~isfinite(sig)) = NaN;
        sig = fillmissing(sig,'pchip');
        sig = fillmissing(sig,'nearest','EndValues','nearest');
        env = abs(hilbert(sig));
    end

end
