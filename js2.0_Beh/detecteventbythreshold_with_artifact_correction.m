function [evt_id, sig_after] = detecteventbythreshold_with_artifact_correction(sig_before, nSamp, decimate_logic)
%This function detects events separately for signal segments that are below
% or above a threshold. The detected events are combined and sorted. This function 
% can be useful for event detections in the presence of regular pulse-like
% artifacts. 

% set the threshold
thres = mean(sig_before);

% get thresholded indices
abv_ind = find(sig_before > thres);
blw_ind = find(sig_before <= thres);

if decimate_logic
    % detect events separately for segments below or abv the threshold with decimation
    sig_abv = decimate(sig_before(abv_ind), nSamp/1000);
    sig_blw = decimate(sig_before(blw_ind), nSamp/1000);
    abv_ind_deci = decimate(abv_ind, nSamp/1000);
    blw_ind_deci = decimate(blw_ind, nSamp/1000); 

    sig_abv_ind = detecteventbythreshold(sig_abv, 1000, 20, 'stdFactor',3, 'plotRez',false, 'chunkPulses', false);
    sig_blw_ind = detecteventbythreshold(sig_blw, 1000, 20, 'stdFactor',3, 'plotRez',false, 'chunkPulses', false);

    evt_id = find(ismember(sort([abv_ind_deci, blw_ind_deci]), [abv_ind_deci(sig_abv_ind), blw_ind_deci(sig_blw_ind)]));
    
    sig_abv_ms = sig_abv - median(sig_abv); 
    
    sig_after_srt = sortrows([[sig_abv_ms; abv_ind_deci], [sig_blw; blw_ind_deci]]', 2);
    sig_after = sig_after_srt(:, 1)'; 


else
    % detect events separately for segments below or abv the threshold without decimation
    sig_abv   = sig_before(abv_ind);
    sig_blw = sig_before(blw_ind);

    sig_abv_ind = detecteventbythreshold(sig_abv, nSamp, 20, 'stdFactor',3, 'plotRez',false, 'chunkPulses', false);
    sig_blw_ind = detecteventbythreshold(sig_blw, nSamp, 20, 'stdFactor',3, 'plotRez',false, 'chunkPulses', false);

    evt_id = find(ismember(sort([abv_ind, blw_ind]), [abv_ind(sig_abv_ind), blw_ind(sig_blw_ind)]));

    sig_after_srt = sortrows([[sig_abv; abv_ind], [sig_blw; blw_ind]]', 2);
    sig_after = sig_after_srt(:, 1)'; 

end

