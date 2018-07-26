function rasterJP(neural, behav, row, col) 
% neural - spike data unit spike times, behav - behavioral event times
% row and col are the num_rows and num_cols of the subplot
% behav = the behav variable you want to align to, can be a matrix
% rasterLogic is boolean; if true draw a raster plot
% Things to do: raster plots with trials sorted (reference 'newperievtraster.m')

popData = neural;           % neural population data
units = 1:length(popData);  % # of unit
alignVar1 = behav;          % the behavioral events to which neural data are aligned

%% Produce raster plot
% Set a gaussian kernel
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15; % gaussian std
[currParams.filter.kernel]     = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

for i = 1:size(alignVar1,2) % # of events alignVar needs to be trials-by-events (of different kinds)
    
    for j = 1:length(units) % # of units
        numStamps = length(popData{1,j});               % # of spikes of the current unit
        delta = zeros(1,ceil(popData{1,j}(numStamps))); % ceil: rounds each element to the nearest integer greater than or equal to that element
        delta(1,ceil(popData{1,j})) = 1;                % delta function of the spike train
        tmpSmooth = conv(delta,currParams.filter.kernel,'same'); % . convultion of the delta function with the Gaussian kernel
        
        psthWin = [1e3,1e3]; % 2 sec window
        [respCSS] = TNC_AlignRasters(tmpSmooth, popData{1,j}, -1, alignVar1(:,i), psthWin, 1, 1); % actual zscore normalized psth % alignVar1(:,i) needs to be all trials of an event of interest
        subplot(row,col,j);
        hold on;
        unitname = num2str(j);
        plotName = strcat('u#', unitname); % subplot name
        
        if j == 1 % for the first unit label the axes
            ylabel('Trials');
            xlabel('Time (ms)');
        end
        
        title(plotName);
        
        % raster plot
        for tr = 1:length(respCSS.raster.trial) % # of trials
            hold on;
            for o = 1:length(respCSS.raster.trial(tr).ts) % # of spikes per trial
                line([respCSS.raster.trial(tr).ts(o) respCSS.raster.trial(tr).ts(o)],[tr-1 tr],'LineWidth',2.0)      % put line-segment at each timestamp
            end
            hold off;
        end
        xlim([-1e3,1e3]) % 2 sec window
        ylim([-1 length(respCSS.raster.trial)+1])
        
    end
    hold off
end

end


