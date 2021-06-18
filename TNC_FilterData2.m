 function [hiBandData] = TNC_FilterData2(data,par)
% FUNCTION DETAILS: general utility to separate data into a pair of bandwidths for spike sorting and continuous analysis. LowBand is 2-0.1k; HiBand: 0.7k-7k
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

tstart = tic;

hiBandData = spike_detection_filter(data, par);
    
telapsed = toc(tstart);

% report the status and elapsed time
disp(sprintf('     EVENT FILTER | Data size (seconds): %g | Elapsed time (seconds): %g',length(data)./par.sr,telapsed));


%Taken from: https://github.com/csn-le/wave_clus/blob/4c74f120c304948ad23973af816978a974af22aa/Batch_files/spike_detection_filter.m

function xf_detect = spike_detection_filter(x, par)
%this function filter the signal, using the detection filter. Is used in the
%readInData class. The filtered data will be downsampled and returned.

sr = par.sr;
fmin_detect = par.detect_fmin;
fmax_detect = par.detect_fmax;


% HIGH-PASS FILTER OF THE DATA
if exist('ellip','file')                         %Checks for the signal processing toolbox
    [b,a] = ellip(2,0.1,40,[fmin_detect fmax_detect]*2/sr);
    if exist('FiltFiltM','file')
        xf_detect = FiltFiltM(b, a, x);
    else
        xf_detect = filtfilt(b, a, x);
    end
else
    xf_detect = fix_filter(x);                   %Does a bandpass filtering between [300 3000] without the toolbox.
end