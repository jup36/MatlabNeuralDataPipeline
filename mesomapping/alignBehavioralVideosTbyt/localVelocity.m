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