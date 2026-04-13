function env = localEnvelope(sig)
    % Envelope via Hilbert; requires finite inputs
    sig(~isfinite(sig)) = NaN;
    sig = fillmissing(sig,'pchip');
    sig = fillmissing(sig,'nearest','EndValues','nearest');
    env = abs(hilbert(sig));
end