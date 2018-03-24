function [ evtST, evtSC ] = perievtSpikeRaster( st, evtts, win, edges )
%perieventFR utilizes the chronux function createdatamatpt to get binned 
% perievent spike counts and firing rates. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
%   - st: spike timestamps of the current unit during the entire recording
%   session
%   - eventmat: event timestamps  
%   - win: uv for createdatamatpt specifying the size of perievent windows
%   - edges: edges for histc
%   - binsize: size of the time bins
% Output: 
%   - evtST: corrected spike timestamps within the perievent timewindows 
%   - evtSC: binned spike counts within the perievent timestamps
%            transformed (needs to be decoded using full)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evtST = cell(1,1);    % event spike timestamps  
evtSC = NaN(length(edges),1);    % event spike counts

evtST = createdatamatpt(st,evtts,win);    % st: spike timestamps, eventts: event timestamps, win: time window
evtST = evtST.times - win(1);    % correct timestamps to align the event time at 0

tempSCbin = histc(evtST,edges);     % get binned spike counts using histc
if isempty(tempSCbin) == 0     % this is to avoid filling in emptry cells
    evtSC = tempSCbin;    % fill evtSC mat with transposed evtSC
    
elseif isnan(evtts) == 0    % in case, the current trial is defined, meaning that event occurred but the unit didn't fire
    evtSC(:,1) = 0;
    
elseif isnan(evtts) == 1    % in case, the current trial is undefined (NaN), most likely to be the missing FTPK due to the end of the block
    % do not feed 0 in this case
end
clearvars temp*


