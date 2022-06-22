function [spkCntPerTrial] = TNC_QuantSpksPerTrial(sumFlag,alignedRaster,window)

if sumFlag==1
    spkCntPerTrial = sum(alignedRaster(:,window(1,1):window(1,2)),2);
else
    spkCntPerTrial = trapz(alignedRaster(:,window(1,1):window(1,2)),2);
end
 