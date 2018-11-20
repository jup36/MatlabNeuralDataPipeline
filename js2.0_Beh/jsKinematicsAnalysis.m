
%filePath='/Volumes/RAID2/parkj/NeuralData/js2.0/WR25/110718_LowHighShift';  
filePath = 'Z:\parkj\NeuralData\js2.0\WR25\110718_LowHighShift'; 

cd(filePath)
load('BehVariablesJs.mat', 'jsTime1k', 'p')

% get meanMass
if isfield(p.Results,'meanMass')
    meanMass = p.Results.meanMass;
else
    meanMass(1,:) = [10 20 30 40 50 60 70 80 90 100]; 
    meanMass(2,:) = [2.95 3.52 4.10 4.67 5.25 5.82 6.40 6.97 7.55 8.12]; 
end
    
S = jsTime1k; 

[cmap1] = cbrewer('div','Spectral',10); 
[cmapg] = cbrewer('seq','Greys',10); 

% demarcate trials
b = 1;
for t = 1:length(S)
    if t==1
        S(t).tBound=b;
    else % t>1
        if ~(S(t).reachP1==S(t-1).reachP1 && S(t).pull_torque==S(t-1).pull_torque) 
            b = b+1; 
        else
        end
        S(t).tBound=b; % trial boundaries
    end
    
end
clearvars t b

tType = {jsTime1k(:).trialType};

%% reaction time
%movKinsPlot(jsTime1k(2).movKins)
[S(:).rt]=deal(NaN); 
actTrIdx = zeros(1,length(S));  
actTrIdx(find([S(:).rewarded]==1,1,'first'):find([S(:).rewarded]==1,1,'last'))=1;

pullTqs = sort(unique([S(:).pull_torque]),'descend'); % pull torques
pullTqsCmap = cmap1(round(linspace(1, length(cmap1), length(pullTqs))),:); 

reachPs = sort(unique([S(:).reachP1])); % reach position 1
reachPsCmap = cmapg(2:length(reachPs)+1,:);  

% detect shift in pullTorque, reachP1
reachP1 = [S(:).reachP1]; 
reachP1s = unique(reachP1); 
reachP1shiftPts = [1 find([reachP1(1) reachP1(1:end-1)]-[S(:).reachP1]~=0)]; 

% draw the reachPosition1 shifts
figure; hold on; 
for s = 1:length(reachP1shiftPts)
    if s<length(reachP1shiftPts)
        tempColor = reachPsCmap(reachP1s==S(reachP1shiftPts(s)).reachP1,:); 
        tempX1 = S(reachP1shiftPts(s)).trStart; 
        tempX2 = max(S(reachP1shiftPts(s+1)-1).trEnd, S(reachP1shiftPts(s+1)-1).rewardT); 
        fill([tempX1 tempX2 tempX2 tempX1]./1000, [0 0 p.Results.trialTimeout p.Results.trialTimeout], tempColor, 'EdgeColor','none'); 
    elseif s==length(reachP1shiftPts) 
        tempColor = reachPsCmap(reachP1s==S(reachP1shiftPts(s)).reachP1,:); 
        tempX1 = S(reachP1shiftPts(s)).trStart; 
        tempX2 = max(S(length(S)).trEnd, S(length(S)).rewardT); 
        fill([tempX1 tempX2 tempX2 tempX1]./1000, [0 0 p.Results.trialTimeout p.Results.trialTimeout], tempColor, 'EdgeColor','none');          
    end
end
clearvars s 
set(gca,'TickDir','out')
%set(gca,'YScale','log')

rtCollect = []; % just collect RT across all active trials 
% get trial-by-trial rt and plot them 
for t = 1:length(S)
    switch S(t).trialType
        case 'sp'
            S(t).rt = S(t).movKins.pullStart;             
        case 'ps'
            S(t).rt = S(t).movKins.pushStart; 
        case 'to'
            S(t).rt = p.Results.trialTimeout; 
        case 'pmpp'
            S(t).rt = S(t).movKins.pull.startI; 
    end
    if actTrIdx(t)
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        if S(t).rewarded
           %figure(1) % 
           %plot(S(t).trJsReady,S(t).pull_torque,'+','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',7); hold on; % Js ready
           %plot(S(t).trJsReady+S(t).rt,S(t).pull_torque,'o','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',7) % response
           %figure(2) %
           plot((S(t).trJsReady+S(t).rt)/1000, S(t).rt,'o','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',6); 
           rtCollect = [rtCollect; S(t).trJsReady+S(t).rt, S(t).rt]; 
        else % not rewarded
           %figure(1) %
           %plot(S(t).trJsReady,S(t).pull_torque,'+','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',7); hold on; % Js ready
           %plot(S(t).trJsReady+S(t).rt,S(t).pull_torque,'o','MarkerEdgeColor',tempColor,'MarkerSize',7) % response 
           %figure(2) %
           plot((S(t).trJsReady+S(t).rt)/1000, S(t).rt,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);
           rtCollect = [rtCollect; S(t).trJsReady+S(t).rt, S(t).rt]; 
        end
    end
end
clearvars t
plot(rtCollect(:,1)/1000,rtCollect(:,2),':k')

rtCell = cell(length(unique([S(:).tBound])),2); 
for tB = unique([S(:).tBound]) 
    tempColor = pullTqsCmap(S(find([S(:).tBound]==tB,1)).pull_torque == pullTqs,:);
    temptB = [S(:).tBound]==tB & actTrIdx; 
    rtCell{tB,1} = [S(temptB).trJsReady]; 
    rtCell{tB,2} = [S(temptB).rt];
    bRt = bar(nanmean(rtCell{tB,1})/1000,nanmean(rtCell{tB,2}),'barWidth',10);
    set(bRt,'FaceColor',tempColor)
    ebRt = errorbar(nanmean(rtCell{tB,1})/1000,nanmean(rtCell{tB,2}), nanstd(rtCell{tB,2})/sqrt(length(rtCell{tB,2})),'CapSize',10,'color',tempColor); 
    set(ebRt,'Color',tempColor)
end
clearvars tB

ylabel('RT(ms)')
xlabel('Time(s)')
print('rt','-dpdf','-bestfit') % save the RT plot 

ylabel('RT(ms)')
xlabel('Time(s)')
set(gca,'yscale','log')
print('rtlogy','-dpdf','-bestfit') % save the RT (log scale) plot
hold off; 

%% reach duration 
[S(:).rd]=deal(NaN);
rdCollect = []; % just collect RD across all active trials 
spIdx = cellfun(@(x)strcmp(x,'sp'),tType); % successful pull index
% get trial-by-trial rd and plot them 
% draw the reachPosition1 shifts
figure; hold on; 
for s = 1:length(reachP1shiftPts)
    if s<length(reachP1shiftPts)
        tempColor = reachPsCmap(reachP1s==S(reachP1shiftPts(s)).reachP1,:); 
        tempX1 = S(reachP1shiftPts(s)).trStart; 
        tempX2 = max(S(reachP1shiftPts(s+1)-1).trEnd, S(reachP1shiftPts(s+1)-1).rewardT); 
        fill([tempX1 tempX2 tempX2 tempX1]./1000, [0 0 1000 1000], tempColor, 'EdgeColor','none'); 
    elseif s==length(reachP1shiftPts) 
        tempColor = reachPsCmap(reachP1s==S(reachP1shiftPts(s)).reachP1,:); 
        tempX1 = S(reachP1shiftPts(s)).trStart; 
        tempX2 = max(S(length(S)).trEnd, S(length(S)).rewardT); 
        fill([tempX1 tempX2 tempX2 tempX1]./1000, [0 0 1000 1000], tempColor, 'EdgeColor','none');          
    end
end
clearvars s 
set(gca,'TickDir','out')

% plot trial-by-trial rd
for t = 1:length(S)
    if strcmp(S(t).trialType,'sp')
       tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:); 
       S(t).rd = S(t).movKins.pullStop-S(t).movKins.pullStart;
       rdCollect = [rdCollect; t, S(t).trJsReady+S(t).rt, S(t).rd]; 
       plot((S(t).trJsReady+S(t).rt)/1000, S(t).rd,'o','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',6);
    end
end
clearvars t
plot(rdCollect(:,2)/1000,rdCollect(:,3),':k')

rdCell = cell(length(unique([S(:).tBound])),2); 
for tB = unique([S(:).tBound]) 
    tempColor = pullTqsCmap(S(find([S(:).tBound]==tB,1)).pull_torque == pullTqs,:);
    temptB = [S(:).tBound]==tB & actTrIdx; 
    rdCell{tB,1} = [S(temptB).trJsReady]; 
    rdCell{tB,2} = [S(temptB).rd];
    bRd = bar(nanmean(rdCell{tB,1})/1000,nanmean(rdCell{tB,2}),'barWidth',10);
    set(bRd,'FaceColor',tempColor)
    ebRd = errorbar(nanmean(rdCell{tB,1})/1000,nanmean(rdCell{tB,2}), nanstd(rdCell{tB,2})/sqrt(length(rdCell{tB,2})),'CapSize',10,'color',tempColor); 
    set(ebRd,'Color',tempColor)
end
clearvars tB
hold off; 
ylabel('RD(ms)')
xlabel('Time(s)')
ylim([0 max(rdCollect(:,3))+10])
print('sp_rd','-dpdf','-bestfit') % save the RT (log scale) plot

%% Max force generated  
[S(:).maxForcePull]=deal(NaN);
mfCollect = []; % just collect RD across all active trials 
% get trial-by-trial rd and plot them 
% draw the reachPosition1 shifts
figure; hold on; 
for s = 1:length(reachP1shiftPts)
    if s<length(reachP1shiftPts)
        tempColor = reachPsCmap(reachP1s==S(reachP1shiftPts(s)).reachP1,:); 
        tempX1 = S(reachP1shiftPts(s)).trStart; 
        tempX2 = max(S(reachP1shiftPts(s+1)-1).trEnd, S(reachP1shiftPts(s+1)-1).rewardT); 
        fill([tempX1 tempX2 tempX2 tempX1]./1000, [-10 -10 1000 1000], tempColor, 'EdgeColor','none'); 
    elseif s==length(reachP1shiftPts) 
        tempColor = reachPsCmap(reachP1s==S(reachP1shiftPts(s)).reachP1,:); 
        tempX1 = S(reachP1shiftPts(s)).trStart; 
        tempX2 = max(S(length(S)).trEnd, S(length(S)).rewardT); 
        fill([tempX1 tempX2 tempX2 tempX1]./1000, [-10 -10 1000 1000], tempColor, 'EdgeColor','none');          
    end
end
clearvars s 
set(gca,'TickDir','out')

% plot trial-by-trial max Force in the PULL DIRECTION
for t = 1:length(S)
    if strcmp(S(t).trialType,'sp')
       tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:); 
       S(t).maxForcePull = S(t).movKins.maxForce;
       mfCollect = [mfCollect; t, S(t).trJsReady+S(t).rt, S(t).maxForcePull]; 
       plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxForcePull,'o','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',6);
    elseif strcmp(S(t).trialType,'to')
       tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:); 
       if ~isempty(find(S(t).movKins.periodicAbsVelSum(1:1000)>50,1,'last'))
           absltStillPt = find(S(t).movKins.periodicAbsVelSum(1:1000)>50,1,'last')+50; 
       else
           absltStillPt = 1;
       end
       S(t).maxForcePull = min(S(t).movKins.forceMN(1,absltStillPt:end));  
       mfCollect = [mfCollect; t, S(t).trJsReady+S(t).rt, S(t).maxForcePull]; 
       plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxForcePull,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);    
       
    elseif strcmp(S(t).trialType,'pmpp')
       tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:); 
       S(t).maxForcePull = S(t).movKins.pull.maxForce;
       mfCollect = [mfCollect; t, S(t).trJsReady+S(t).rt, S(t).maxForcePull]; 
       plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxForcePull,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);
    end
end
clearvars t
plot(mfCollect(:,2)/1000,-mfCollect(:,3),':k')

mfCell = cell(length(unique([S(:).tBound])),2); 
for tB = unique([S(:).tBound]) 
    tempColor = pullTqsCmap(S(find([S(:).tBound]==tB,1)).pull_torque == pullTqs,:);
    temptB = [S(:).tBound]==tB & actTrIdx; 
    mfCell{tB,1} = [S(temptB).trJsReady]; 
    mfCell{tB,2} = -[S(temptB).maxForcePull];
    bMf = bar(nanmean(mfCell{tB,1})/1000,nanmean(mfCell{tB,2}),'barWidth',10);
    set(bMf,'FaceColor',tempColor)
    ebMf = errorbar(nanmean(mfCell{tB,1})/1000,nanmean(mfCell{tB,2}), nanstd(mfCell{tB,2})/sqrt(length(mfCell{tB,2})),'CapSize',10,'color',tempColor); 
    set(ebMf,'Color',tempColor)
end
clearvars tB
hold off; 
ylabel('max Force (mN)')
xlabel('Time(s)')
ylim([0 -min(mfCollect(:,3))+10])
print('mf','-dpdf','-bestfit') % save the RT (log scale) plot










