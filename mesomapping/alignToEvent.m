function [evtAlignedDffTs, evtAlignedTs] = alignToEvent(dffTs, eventTimeToAlign, frameT, timeWin)
        
        win = timeWin + eventTimeToAlign; 

        if min(win) >= min(frameT) && max(win) <= max(frameT)    
          frameI = frameT >= min(win) & frameT <= max(win); 
          evtAlignedDffTs = dffTs(frameI); 
          evtAlignedDffTs = evtAlignedDffTs(:).'; % make it a row vector
          evtAlignedTs = frameT(frameI);  
          evtAlignedTs = evtAlignedTs(:).'; % make it a row vector  
        else
          warning("The peri-event window goes out of bound!")
          evtAlignedDffTs = []; 
          evtAlignedTs = []; 
        end
       
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [imagStackMaskedMean, imgStackMasked] = apply2DMaskTo3DStack(imgStack, mask2D)
%     % Replicate the 2D mask to match the 3D stack dimensions
%     replicatedMask = repmat(mask2D, [1, 1, size(imgStack, 3)]);
% 
%     % Apply the mask: replace values in the stack where the mask is zero (or false) with NaN
%     imgStack(~replicatedMask) = NaN;
% 
%     imgStackMasked = imgStack;
%     
%     imagStackMaskedMean = squeeze(nanmean(nanmean(imgStackMasked, 1), 2)); 
% 
% end

end