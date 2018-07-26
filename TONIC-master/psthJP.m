function [respCSS, psthZscore] = psthJP(neural, behav, row, col, psthPlotFlag) % neural spike data unit spike times, behav behavioral event times
%This function takes spike times data (a cell array) and the behavioral
% timestamps (a vector or a matrix) and generates psths aligned to the
% behavioral timestamps, and also generates subplots each depicting each
% unit's psth
% row and col are the num_rows and num_cols of the subplots each depicting
% each unit's psth
% behav = the behav variable you want to align to, can be a matrix
% rasterLogic is boolean; if true draw a raster plot
% modified from Josh's PSTH code
% Things to do: 
%

units = 1:length(neural);  % # of unit
alignVar1 = behav;          % the behavioral events to which neural data are aligned

%% Produce PSTH
% Set a gaussian kernel
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15; % gaussian std
[currParams.filter.kernel]     = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

for j = 1:length(units) % # of units
    numStamps = length(neural{1,j});               % # of spikes of the current unit
    delta = zeros(1,ceil(neural{1,j}(numStamps))); % delta function of the spike train
    delta(1,ceil(neural{1,j})) = 1;                % delta function of the spike train
    tmpSmooth = conv(delta,currParams.filter.kernel,'same'); % convultion of the delta function with the Gaussian kernel
    
    psthWin = [1e3,1e3]; % 2 sec window
    [respCSS, psthZscore(j,:)] = TNC_AlignRasters(tmpSmooth, neural{1,j}, -1, alignVar1, psthWin, 1, 1); % actual zscore normalized psth % alignVar1(:,i) needs to be all trials of an event of interest
    % TNC_BinSpikeCount
    
    if psthPlotFlag % if the boolean for psth plot is true, draw psths
        subplot(row,col,j);
        hold on;
        unitname = num2str(j);
        plotName = strcat('u#', unitname); % subplot name
        plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2])
        ylim([-.3 1.5]); % ylim zscores
        
        if j == 1 % for the first unit label the axes
            ylabel('Norm. Firing Rate');
            xlabel('Time (ms)');
        end
        
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');  % get shadedErrorBar function
        
        line([0 0],[-.5 2.5]);
        title(plotName);
    else    % do not draw psths
    end
end
hold off

end


