% FUNCTION DETAILS: a script to run through the set of analysis routines to examine the evoked field potential response to stimulation.
% This script will operate on whatever data is stored in the structure "data". To use with multiple data files it is recommended that one loads as many structures into memory as desired and stores inidivudal data structures as "data" for the purpose of executing the script. It will overwrite variables however so one should be careful.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% 

% first, normalize the data to a given reference channel and convert to z-scores
normedData = TNC_NormDataMat(data.Data,1,1);

% using the stimulus data obtained from an analog channel, collect aligned
% responses **Might need to check this and ensure that there are similar patterns of stim events**
[stim] = TNC_ConvertContToEvent(data.Data(68,:),68);

% build up a cell array of repeated trial responses
[PeWin] = TNC_BuildSegsFromMemory(normedData,stim.indices,[0,1500],1:1:64);

% calculate the mean response for each matrix of trials in PeWin
[meanDataMat] = TNC_GetMeanData(PeWin);

% examine the mean response as a function of time across the population
% NOTE: here some channel remapping may be critical, but not yet positive it is correctly implemented
[waveRep] = TNC_EvolveContOverTime([8 8],meanDataMat,[1 1000]);

% Some plotting functions to examine the mean response and how it evolves in time across the array
TNC_PlotEachCh([8 8],meanDataMat);
TNC_WatchAllChanEvolve(waveRep);
