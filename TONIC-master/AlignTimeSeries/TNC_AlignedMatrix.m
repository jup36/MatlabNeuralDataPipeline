function [alignedMatrix] = TNC_AlignedMatrix(eventTimes,channels,spikeTimes,window,sampling)
% FUNCTION DETAILS: For each channel create a field in the structure 'alignedMatrix' that has an matrix of SPIKETIMES aligned to the list of EVENTIMES passed in and that extends from -WINDOW(1,1) to WINDOW(1,2) with SAMPLING precision
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: josh@dudmanlab.org
% CONTRIBUTIONS: dudmanlab.org/html/projects.html
% _________________________________________________________________________

% cycle through the list eventTimes