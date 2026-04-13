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
