% _________________________________________________________________________
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / Janelia
% 
% BUG REPORTING: josh@dudmanlab.org
% FURTHER INFORMATION: www.dudmanlab.org
% Copyright (C) 2015 by Howard Hughes Medical Institute.
% _________________________________________________________________________
% _________________________________________________________________________
% 
% DETAILS:
% Use online sorted data to generate a template and then re-search the data to try and extract the template matches in the data to get a full sorted sessions
% 
% INPUTS:
% segArray: array of segments that sorting will be applied over
% seed_segs: array of segments that should be used to generate the template
% shank_num: what shank to perform sorting on
% file_Path: basename of the files used ('_ft_tns.mat' '_ss.mat' '_ft.mat' are required)
% updateTnsFlag: 1=overwrite tns file with new sorted indices
% 
% EXAMPLE:
% [shank] = TNC_ExtendManualSort([1:2],[1],5,'DA_F06_20140815-i16',1);
% 
% OUTPUTS:
% Can overwrite the spike ids file (*_tns) used for sorting in GUI

function [] = TNC_ExtendOnlineSort(list_of_units,baseFileNameStr,updateTnsFlag)

baseFileNameStr = 'Gad2t_150202002';

% call the previewer function
[PrevData] = TNC_SSPL_NevPreviewer([baseFileNameStr '.nev'],0);

% load the original data
[data] = TNC_LoadData(0, 0, [baseFileNameStr '.nev']);
    
% ask the user for the list of valid looking units to go through and try to improve
list_of_units = [1 2 3 5 6 7 8 9 11 13 14 17 19 20]

% use the mean waveform as a template and run through the data trying to
% extract all the spikes matching that template

for i=1:numel(list_of_units)
    
    % what unit number?
    m = list_of_units(i)
    
    
    % what electrode?
    elecNumb = PrevData.unit(m).el;
    
%     % load the continuous data (.ns5)
%     [data_30k] = TNC_LoadData(['e:' num2str(elecNumb,0,[baseName '.ns5']);
    
    % get all events
    toMatchInds = find(data.Data.Spikes.Electrode==elecNumb);
    
    % get a matrix of all waveforms
    toMatch = double(data.Data.Spikes.Waveform(toMatchInds,:));
    template = mean(PrevData.unit(m).wfs);
    templateE= std(PrevData.unit(m).wfs,[],1);
%     template2 = sinc([-15:32]./2);
    
    % convolve with the template
    heuristics = toMatch*template';
%     heuristics2 = -750 .* (toMatch*template2');
    amplitude = log10( max(toMatch,[],2) - min(toMatch,[],2) );
    
    % look for a cluster in that distribution
    figure(21); hold off; shadedErrorBar([-15:32] , template , templateE, {'Color' , [0.5 0.5 0.5]}); 
    figure(18); hist(amplitude,2:0.001:4);
    figure(17); hist(heuristics,10000);
    
%     figure(22); plot(amplitude , heuristics , 'k.');
    
    % ask user for the threshold to apply...
    [userThresh,Y] = ginput(2); % user should specify the bounds
    
    % take all suprathreshold events and create a new 
    matchedInds = toMatchInds(find(heuristics>userThresh(1) & heuristics<userThresh(2)));
    
    % update the data matrix
    PrevData2.unit(m).wfs = double(data.Data.Spikes.Waveform(matchedInds,:));
    PrevData2.unit(m).ts  = ceil(double(data.Data.Spikes.Timestamps(matchedInds) ./ 30));
    PrevData2.unit(m).el  = PrevData.unit(m).el;
    PrevData2.unit(m).un  = PrevData.unit(m).un;   
    
    figure(21); hold off; shadedErrorBar([-15:32] , template , templateE, {'Color' , [0.5 0.5 0.5]}); 
    hold on; 
    shadedErrorBar([-15:32] , mean(PrevData2.unit(m).wfs,1) , std(PrevData2.unit(m).wfs,[],1) , {'Color' , [1 0 0]});

    
end

%% Also could just use this as a way to scan all channels?