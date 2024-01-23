function baseSubSig = baseSubNorm(sig, sigTs, baseWin)
baseMean = nanmean(sig(sigTs >= min(baseWin) & sigTs <= max(baseWin)));
baseSubSig = sig-baseMean;
end