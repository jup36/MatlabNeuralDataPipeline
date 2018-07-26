function [data] = TNC_LoadData(channel, analog, fileNameStr)
% FUNCTION DETAILS: Function that loads multiple data formats into TONIC; Can load: Blackrock, MAT structures; FUTURE SUPPORT: Neuralynx, Plexon?, APIG, NeuroExplorer, OFS
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
%
% TO DO: add in some checking to make sure this is data that is prepped for
% TONIC analysis

loadTrue = 0;

if channel==0
    % read all channels into memory
    allChanFlag = 1;
else
    if analog==0
        electrodeChannels = sprintf('e:%g',channel);
    else
        electrodeChannels = sprintf('ainp:%g',channel);
    end

end

lastloaded = fileNameStr;

strindex=findstr(fileNameStr,'.mat');
if ~isempty(strindex);
    data    = load(fileNameStr);
end


strindex=findstr(fileNameStr,'.nex');
if ~isempty(strindex);
    %     disp('By default nex loading is for the entire data structure; channel by channel selection should occur using extraction protocols');
    data = readNexFile(fileNameStr);
    loadTrue = 1;
end

strindex=findstr(fileNameStr,'.nev');
if ~isempty(strindex);
    data = openNEV(fileNameStr,'read','noparse','nomat','nosave');
    loadTrue = 1;
end

strindex=findstr(fileNameStr,'.ns3');
if ~isempty(strindex);
    data = openNSx(fileNameStr,'read','report');
    loadTrue = 1;
end


strindex=findstr(fileNameStr,'.ns4');
if ~isempty(strindex);
    data = openNSx(fileNameStr,'read','report');
    loadTrue = 1;
end

strindex=findstr(fileNameStr,'.ns5');
if ~isempty(strindex);
    if allChanFlag==0
        data = openNSx(fileNameStr,'read','report',electrodeChannels);
        data.channel = channel;
    else
        data = openNSx(fileNameStr,'read','report');
    end
    loadTrue = 1;
end

if (~loadTrue)
    error('The file you have attempted to load is not currently supported on this system');
end

data.fileName = fileNameStr;