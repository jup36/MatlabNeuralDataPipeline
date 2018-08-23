%This script performs DA cell classification

clear all; close all; clc;
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab\Functions'));

cd('C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT')
load('VTA_sua_dat_022115.mat')

cd('C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT\VTA_cell_classification')
load('SUA_AGB_VTA_cellclass_50ms_step_092315','EM')

%% these are variables that the function will take passed in
Fs = 40000;     % 40 kHz sampling rate
wf = cell2mat(dat(:,7)); % waveforms
avgFR = cell2mat(dat(:,5)); % waveforms
frThresh = 10;
scrsz = get(0,'ScreenSize');

%% find firing rates less than classification threshold
wfClass = avgFR <= frThresh; % cells with their FRs < 10 Hz

%% manually select the start, peak, valley, and end of waveforms
counter = 0;
for u = 1:size(wf,1)    % # of units
    %% select the waveform and the zero crossings of the derivative
    counter = counter +1;
    spike = wf(u,:);     % the mean spike waveform
    slope = diff(spike); % differences between elements (t1 - t0)
    inflec = crossing(slope); % Does the slope change its sign? (From positive to negative or vice versa)
    idx = spike >= 0;    % index for the positive voltage values in the spike waveform
    xValVec = [1:length(spike)]; % x value vector
    
    %% instruct user at the command line
    disp(['Determine 1) Triphasic 2) Non-triphasic'])
    disp(['For Triphasic WF Select: 1) Start 2) Peak 3) Valley 4) 2nd Peak 5) Stop'])
    disp(['For Non-Triphasic WF Select: 1) Start 2) Peak 3) Valley 4) Stop'])
    disp(['Once you select critical points, they turn blue. Press space to accept.Plus to redo.'])
    disp(['Unit #',num2str(u)])
    disp(['Inflections = ',num2str(inflec)])
    
    %% display figures
    close all;
    figure('Position',[scrsz])
    subplot(2,1,1)
    plot(xValVec(idx),spike(idx),'ko',xValVec(~idx),spike(~idx),'co','LineWidth',3); % waveform in black, negative points in cyan
    hold on
    plot(spike,'k','LineWidth',2) % waveform in black
    hold on
    plot(inflec+1,spike(inflec+1),'ro','LineWidth',4)
    subplot(2,1,2)
    plot(spike,'k-o','LineWidth',2);
    hold on
    plot(inflec+1,spike(inflec+1),'ro','LineWidth',4)
    line([1 length(spike)],[0 0]);
    ylim([-.2 .2])
    
    %% get user inputs and verify them
    tnt = input('Enter 1 for Triphasic WF, Enter 2 for Non-triphasic WF!'); % Enter 1 for triphasic, or 2 for non-triphasic
    
    if tnt == 1 % in case of a triphasic waveform
        [x,y] = ginput(5); % graphical input from mouse: 1) start, 2) 1st peak, 3) valley, 4) 2nd peak, 5) stop
        fail = true;
        while fail == true;
            if y(2)< 0 || y(3) > 0 || ~isempty(find(diff(x)<0)) % in case of an erroneous input
                disp('ERROR! Select: 1) Start 2) Peak 3) Valley 4) 2nd Peak 5) Stop')
                clear x y
                fail = true;
                [x,y] = ginput(5);
            else % in case of a valid input, show the users their choices
                estTime = round(x);
                subplot(2,1,1)
                hold on
                plot(estTime,spike(estTime),'bo','LineWidth',4)
                subplot(2,1,2)
                hold on
                plot(estTime,spike(estTime),'bo','LineWidth',4)
                
                % ask user to verify that the points are correct
                disp(['You selected time points:   ',num2str(estTime')]);
                result = input('Do you accept these measurements (press enter to accept and X to invalidate start/stop times)','s');
                if strcmp(result,'') % Accept the input and stop the loop
                    fail = false;
                elseif strcmp(result,'x') % NaN out the input and stop the loop
                    disp(['You elected to not use the waveform measurements for unit:',num2str(u)]);
                    pause
                    fail = false;
                    x(:)=NaN;    % put NaNs for time points
                else    % Request the user to re-select the input
                    fail = true;
                    close all;clc;
                    disp(['Unit #',num2str(u)])
                    disp(['Inflections = ',num2str(inflec)])
                    % display figures
                    figure('Position',[scrsz])
                    subplot(2,1,1)
                    plot(xValVec(idx),spike(idx),'ko',xValVec(~idx),spike(~idx),'co','LineWidth',3); % waveform in black, negative points in cyan
                    hold on
                    plot(spike,'k','LineWidth',2) % waveform in black
                    hold on
                    plot(inflec+1,spike(inflec+1),'ro','LineWidth',4)
                    subplot(2,1,2)
                    plot(spike,'k-o','LineWidth',2);
                    hold on
                    plot(inflec+1,spike(inflec+1),'ro','LineWidth',4)
                    line([1 length(spike)],[0 0]);
                    ylim([-.2 .2])
                    [x,y] = ginput(5);
                end
            end
        end
        
        % once user accepts the points, store them in
        startTime(counter,1) = round(x(1));
        peakTime(counter,1) = round(x(2));
        valleyTime(counter,1) = round(x(3));
        peakTime2(counter,1) = round(x(4)); % There's 2nd peak in a non-triphasic waveform
        endTime(counter,1) = round(x(5));
        allPoints(counter,:) = x;
        pkVal(counter,1) = spike(peakTime(counter));
        valVal(counter,1) = spike(valleyTime(counter));
        
        if isempty(crossing(spike(peakTime2(counter,1)),spike(endTime(counter,1))))     % if there's no zero crossing between the 2nd peak and the endpoint, this indicates the waveform curtailing due to the small window size
            estslope = (spike(endTime(counter,1))-spike(peakTime2(counter,1)))/(endTime(counter,1)-peakTime2(counter,1));
            estx = -(spike(peakTime2(counter,1))/estslope);
            peak2enddist = round(estx); % estimated distance from the peak 2 to the endpoint
            if peak2enddist > (endTime(counter,1)-peakTime2(counter,1)) % in case the estimated peak2 to endpoint duration is greater than that of the curtailed waveform
               estendTime(counter,1) = peakTime2(counter,1) + peak2enddist; % estimated endTime;
            else
               estendTime(counter,1) = round(x(5)); % estimated endTime 
            end
        else
            estendTime(counter,1) = round(x(5)); % estimated endTime
        end
            
    elseif tnt == 2 % in case of a non-triphasic waveform
        [x,y] = ginput(4); % graphical input from mouse 1) start, 2) peak, 3) valley, 4) stop
        fail = true;
        while fail == true;
            if y(2)< 0 || y(3) > 0 || ~isempty(find(diff(x)<0)) % in case of an erroneous input
                disp('ERROR! Select: 1) Start 2) Peak 3) Valley 4) Stop')
                clear x y
                fail = true;
                [x,y] = ginput(4);
            else % in case of a valid input, show the users their choices
                estTime = round(x);
                subplot(2,1,1)
                hold on
                plot(estTime,spike(estTime),'bo','LineWidth',4)
                subplot(2,1,2)
                hold on
                plot(estTime,spike(estTime),'bo','LineWidth',4)
                
                % ask user to verify that the points are correct
                disp(['You selected time points:   ',num2str(estTime')]);
                result = input('Do you accept these measurements (press enter to accept and X to invalidate start/stop times)','s');
                if strcmp(result,'') % Accept the input and stop the loop
                    fail = false;
                elseif strcmp(result,'x') % NaN out the input and stop the loop
                    disp(['You elected to not use the waveform measurements for unit:',num2str(u)]);
                    pause
                    fail = false;
                    x(:)=NaN;    % put NaNs for time points
                else    % Request the user to re-select the input
                    fail = true;
                    close all;clc;
                    disp(['Unit #',num2str(u)])
                    disp(['Inflections = ',num2str(inflec)])
                    % display figures
                    figure('Position',[scrsz])
                    subplot(2,1,1)
                    plot(xValVec(idx),spike(idx),'ko',xValVec(~idx),spike(~idx),'co','LineWidth',3); % waveform in black, negative points in cyan
                    hold on
                    plot(spike,'k','LineWidth',2) % waveform in black
                    hold on
                    plot(inflec+1,spike(inflec+1),'ro','LineWidth',4)
                    subplot(2,1,2)
                    plot(spike,'k-o','LineWidth',2);
                    hold on
                    plot(inflec+1,spike(inflec+1),'ro','LineWidth',4)
                    line([1 length(spike)],[0 0]);
                    ylim([-.2 .2])
                    [x,y] = ginput(4);
                end
            end
        end
        
        % once user accepts the points, store them in
        startTime(counter,1) = round(x(1));
        peakTime(counter,1) = round(x(2));
        valleyTime(counter,1) = round(x(3));
        peakTime2(counter,1) = NaN; % There's no 2nd peak in a non-triphasic waveform
        endTime(counter,1) = round(x(4));
        allPoints(counter,:) = [x(1), x(2), x(2), NaN, x(4)];
        pkVal(counter,1) = spike(peakTime(counter));
        valVal(counter,1) = spike(valleyTime(counter));
        estendTime(counter,1) = round(x(4)); % estimated endTime
    end
end
close all; clc;
clear counter s u spike slope inflec x y fail estTime idx result scrsz xValVec

%% update endTime
upendTime = nan(size(endTime,1),size(endTime,2));
for u = 1:length(upendTime)
    if estendTime(u) > endTime(u) && estendTime(u) < 100 % in case the estimated endtime is greater than the endtime, and the estimated endtime is smaller than 100   
        upendTime(u,1) = estendTime(u); % use the estimated endtime!
    else
        upendTime(u,1) = endTime(u); % use the original endtime!
    end
end
clearvars u

%% calculate waveform statistics
totalDur = (upendTime - startTime)*(1/Fs)*1000; % in msec
% startValleyDur = (valleyTime - startTime)*(1/Fs);
% peakValleyDur = (valleyTime - peakTime)*(1/Fs);
% valleyEndDur = (endTime - valleyTime) * (1/Fs);

%% Split Variables in Regions
durThresh = 0.0014;
DAIdx1 = totalDur >= 1.2 & avgFR <= 15; % DA cell criteria: SpikeWidth greater than 1.2 ms, avgFR less than 10 Hz
DAIdx2 = EM.c'==1; % DA cell criteria based on the reward response 
DAIdx = DAIdx1 & DAIdx2; % DAIdx, intersection of two criteria

%% Save
save('VTA_SUA_unitClassification_120115')

%% Plot
hold on;
plot(totalDur(DAIdx),avgFR(DAIdx),'o','MarkerSize',5,'MarkerFaceColor',[76 135 198]./255,'MarkerEdgeColor','none')
plot(totalDur(~DAIdx),avgFR(~DAIdx),'o','MarkerSize',5,'MarkerFaceColor',[220 30 62]./255,'MarkerEdgeColor','none')
xlim([0.6 2.3])
ylim([0 50]) 
set(gca,'TickDir','out')
grid on
%set(gca,'YTick',[0:5:50],'FontSize',14);




