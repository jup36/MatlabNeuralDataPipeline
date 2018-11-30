function jsKinematicsAnalysis(filePath)
%jsKinematicsAnalysis loads the 'jsTime1k' and 'p', and computes
% trial-by-trial kinetic and kinematic variables such as rt, rd, velocity, 
% acceleration, force etc., and generates figures illustrating these variables.  

%filePath='/Volumes/RAID2/parkj/NeuralData/js2.0/WR28/111118';
%filePath = 'Z:\parkj\NeuralData\js2.0\WR25\110718_LowHighShift';

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
[S(:).smJsTrajFull] = deal(NaN); % smoothed Js trajectory 1 to pull/push stop
[S(:).smJsTrajStartMaxV] = deal(NaN);
[S(:).smJsVelStartStop] = deal(NaN);
[S(:).smJsVelStartMaxV] = deal(NaN);
[S(:).smJsAclFull] = deal(NaN);
[S(:).smJsAclStartStop] = deal(NaN);
[S(:).smJsAclStartMaxV] = deal(NaN);
[S(:).smJsTrajStartStop] = deal(NaN);
[S(:).smJsVelFull] = deal(NaN);

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

% trial type indices
tType = {S(:).trialType};
spIdx = cellfun(@(x)strcmp(x,'sp'),tType); % successful pull index
toIdx = cellfun(@(x)strcmp(x,'to'),tType); % timeout index
psIdx = cellfun(@(x)strcmp(x,'ps'),tType); % push index
pmIdx = cellfun(@(x)strcmp(x,'pm'),tType); % premature pull
pmppIdx = cellfun(@(x)strcmp(x,'pmpp'),tType); % premature pull and push index
actIdx = zeros(1,length(S)); % active trial index
actIdx(find([S(:).rewarded]==1,1,'first'):find([S(:).rewarded]==1,1,'last'))=1; % active trial index

%% reaction time
%movKinsPlot(jsTime1k(2).movKins)
[S(:).rt]=deal(NaN);

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
    if actIdx(t)
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        if S(t).rewarded
            plot((S(t).trJsReady+S(t).rt)/1000, S(t).rt,'o','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',6);
            rtCollect = [rtCollect; S(t).trJsReady+S(t).rt, S(t).rt];
        else % not rewarded
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
    temptB = [S(:).tBound]==tB & actIdx;
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
axis tight
ylim([0 p.Results.trialTimeout])
print('rt','-dpdf','-bestfit') % save the RT plot

ylabel('RT(ms)')
xlabel('Time(s)')
ylim([0 p.Results.trialTimeout])
set(gca,'yscale','log')
print('rtlogy','-dpdf','-bestfit') % save the RT (log scale) plot
hold off;

%% reach duration
[S(:).rd]=deal(NaN);
rdCollect = []; % just collect RD across all active trials
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
    temptB = [S(:).tBound]==tB & actIdx;
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
axis tight
ylim([0 max(rdCollect(:,3))+10])
print('sp_ReachDur','-dpdf','-bestfit') % save the RT (log scale) plot

%% Max force generated in the PULL DIRECTION
[S(:).maxForcePull]=deal(NaN);
mfCollect = []; % just collect max forces across all active trials
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
    if spIdx(t)&&actIdx(t)
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).maxForcePull = S(t).movKins.maxForce;
        mfCollect = [mfCollect; t, S(t).trJsReady+S(t).rt, S(t).maxForcePull];
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxForcePull,'o','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',6);
    elseif toIdx(t)&&actIdx(t)
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        if ~isempty(find(S(t).movKins.periodicAbsVelSum(1:1000)>50,1,'last'))
            absltStillPt = find(S(t).movKins.periodicAbsVelSum(1:1000)>50,1,'last')+50;
        else
            absltStillPt = 1;
        end
        S(t).maxForcePull = nanmin(S(t).movKins.forceMN(1,absltStillPt:end));
        mfCollect = [mfCollect; t, S(t).trJsReady+S(t).rt, S(t).maxForcePull];
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxForcePull,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);
        
    elseif pmppIdx(t)&&actIdx(t)
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).maxForcePull = S(t).movKins.pull.maxForce;
        mfCollect = [mfCollect; t, S(t).trJsReady+S(t).rt, S(t).maxForcePull];
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxForcePull,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);
    end
end
clearvars t
plot(mfCollect(:,2)/1000,-mfCollect(:,3),':k')

% plot errorbars per block
maCell = cell(length(unique([S(:).tBound])),2);
for tB = unique([S(:).tBound])
    tempColor = pullTqsCmap(S(find([S(:).tBound]==tB,1)).pull_torque == pullTqs,:);
    temptB = [S(:).tBound]==tB & actIdx;
    maCell{tB,1} = [S(temptB).trJsReady];
    maCell{tB,2} = -[S(temptB).maxForcePull];
    bMf = bar(nanmean(maCell{tB,1})/1000,nanmean(maCell{tB,2}),'barWidth',10);
    set(bMf,'FaceColor',tempColor)
    ebMf = errorbar(nanmean(maCell{tB,1})/1000,nanmean(maCell{tB,2}), nanstd(maCell{tB,2})/sqrt(length(maCell{tB,2})),'CapSize',10,'color',tempColor);
    set(ebMf,'Color',tempColor)
end
clearvars tB
hold off;
ylabel('max Force (mN)')
xlabel('Time(s)')
axis tight
ylim([0 -min(mfCollect(:,3))+10])
print('maxF','-dpdf','-bestfit') % save the RT (log scale) plot

%% Max acceleration generated in the PULL direction
[S(:).maxAccelPull]=deal(NaN);
maCollect = []; % just collect max accelerations across all active trials
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

% plot trial-by-trial max acceleration in the PULL DIRECTION
for t = 1:length(S)
    if spIdx(t)&&actIdx(t)
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).maxAccelPull = S(t).movKins.maxForce/S(t).movKins.mass;
        maCollect = [maCollect; t, S(t).trJsReady+S(t).rt, S(t).maxAccelPull];
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxAccelPull,'o','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',6);
    elseif toIdx(t)&&actIdx(t)
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        if ~isempty(find(S(t).movKins.periodicAbsVelSum(1:1000)>50,1,'last'))
            absltStillPt = find(S(t).movKins.periodicAbsVelSum(1:1000)>50,1,'last')+50;
        else
            absltStillPt = 1;
        end
        S(t).maxAccelPull = nanmin(S(t).movKins.forceMN(1,absltStillPt:end))/S(t).movKins.mass;
        maCollect = [maCollect; t, S(t).trJsReady+S(t).rt, S(t).maxAccelPull];
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxAccelPull,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);
        
    elseif pmppIdx(t)&&actIdx(t)
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).maxAccelPull = S(t).movKins.pull.maxForce/S(t).movKins.mass;
        maCollect = [maCollect; t, S(t).trJsReady+S(t).rt, S(t).maxAccelPull];
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxAccelPull,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);
    end
end
clearvars t
plot(maCollect(:,2)/1000,-maCollect(:,3),':k')

% plot errorbars per block
maCell = cell(length(unique([S(:).tBound])),2);
for tB = unique([S(:).tBound])
    tempColor = pullTqsCmap(S(find([S(:).tBound]==tB,1)).pull_torque == pullTqs,:);
    temptB = [S(:).tBound]==tB & actIdx;
    maCell{tB,1} = [S(temptB).trJsReady];
    maCell{tB,2} = -[S(temptB).maxAccelPull];
    bMf = bar(nanmean(maCell{tB,1})/1000,nanmean(maCell{tB,2}),'barWidth',10);
    set(bMf,'FaceColor',tempColor)
    ebMf = errorbar(nanmean(maCell{tB,1})/1000,nanmean(maCell{tB,2}), nanstd(maCell{tB,2})/sqrt(length(maCell{tB,2})),'CapSize',10,'color',tempColor);
    set(ebMf,'Color',tempColor)
end
clearvars tB
hold off;
ylabel('max Accel (m/s^2)')
xlabel('Time(s)')
axis tight
ylim([0 -min(maCollect(:,3))+5])
print('maxA','-dpdf','-bestfit') % save the RT (log scale) plot

%% Max velocity generated in the PULL direction
[S(:).maxVelPull]=deal(NaN);
mvCollect = []; % just collect max velocitys across all active trials
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

% plot trial-by-trial max velocity in the PULL DIRECTION
for t = 1:length(S)
    if spIdx(t)&&actIdx(t) % successful pull trials
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).maxVelPull = S(t).movKins.pullMaxVel;
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxVelPull,'o','MarkerFaceColor',tempColor,'MarkerEdgeColor',tempColor,'MarkerSize',6);
    elseif toIdx(t)&&actIdx(t) % timeout trials
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        if ~isempty(find(S(t).movKins.periodicAbsVelSum(1:1000)>50,1,'last'))
            absltStillPt = find(S(t).movKins.periodicAbsVelSum(1:1000)>50,1,'last')+50;
        else
            absltStillPt = 1;
        end
        S(t).maxVelPull = nanmin(S(t).movKins.smJsVel(1,absltStillPt:end));
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxVelPull,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);
    elseif pmppIdx(t)&&actIdx(t) % premature pull and push
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).maxVelPull = S(t).movKins.pull.maxVel;
        plot((S(t).trJsReady+S(t).rt)/1000, -S(t).maxVelPull,'o','MarkerEdgeColor',tempColor,'MarkerSize',6);
    end
    mvCollect = [mvCollect; t, S(t).trJsReady+S(t).rt, S(t).maxVelPull];
end
clearvars t
plot(mvCollect(:,2)/1000,-mvCollect(:,3),':k')

% plot errorbars per block
mvCell = cell(length(unique([S(:).tBound])),2);
for tB = unique([S(:).tBound])
    tempColor = pullTqsCmap(S(find([S(:).tBound]==tB,1)).pull_torque == pullTqs,:);
    temptB = [S(:).tBound]==tB & actIdx;
    mvCell{tB,1} = [S(temptB).trJsReady];
    mvCell{tB,2} = -[S(temptB).maxVelPull];
    bMf = bar(nanmean(mvCell{tB,1})/1000,nanmean(mvCell{tB,2}),'barWidth',10);
    set(bMf,'FaceColor',tempColor)
    ebMf = errorbar(nanmean(mvCell{tB,1})/1000,nanmean(mvCell{tB,2}), nanstd(mvCell{tB,2})/sqrt(length(mvCell{tB,2})),'CapSize',10,'color',tempColor);
    set(ebMf,'Color',tempColor)
end
clearvars tB
hold off;
ylabel('max Vel (mm/s)')
xlabel('Time(s)')
axis tight
ylim([0 -min(mvCollect(:,3))+5])
print('maxVel','-dpdf','-bestfit') % save the RT (log scale) plot

%% Collect/plot pos, velocity, acceleration curves trial-by-trial
% Position 3-d plot
spTrCnt = 0;
figure; hold on;
for t = 1:length(S)
    if spIdx(t)&&actIdx(t) % successful pull trials
        spTrCnt = spTrCnt + 1;
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).smJsTrajFull = S(t).movKins.smJsTraj(1:S(t).movKins.pullStop);
        S(t).smJsTrajStartStop = S(t).movKins.smJsTraj(S(t).movKins.pullStart:S(t).movKins.pullStop);
        S(t).smJsTrajStartMaxV = S(t).movKins.smJsTraj(S(t).movKins.pullStart:S(t).movKins.pullMaxVelI);
        
        tempMaxVelAlignedPos = S(t).smJsTrajFull(S(t).movKins.pullMaxVelI-min(200,S(t).movKins.pullMaxVelI)+1:S(t).movKins.pullMaxVelI);
        smJsPosMatMaxVelAlign(spTrCnt,200-length(tempMaxVelAlignedPos)+1:200)=tempMaxVelAlignedPos;
        
        if isequal(S(t).pull_torque == pullTqs,[true,false])
            plot3(1:length(S(t).smJsTrajStartStop),ones(1,length(S(t).smJsTrajStartStop))*t,S(t).smJsTrajStartStop,'Color',tempColor)
        else
            plot3(1:length(S(t).smJsTrajStartStop),-ones(1,length(S(t).smJsTrajStartStop))*t,S(t).smJsTrajStartStop,'Color',tempColor)
        end
        
    end
end
hold off;
grid on

% Velocity 3-d plot
smJsVelMatMaxVelAlign = [];
spTrCnt = 0;
figure; hold on;
for t = 1:length(S)
    if spIdx(t)&&actIdx(t) % successful pull trials
        spTrCnt = spTrCnt + 1;
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).smJsVelFull = S(t).movKins.smJsVel(1:S(t).movKins.pullStop);
        S(t).smJsVelStartStop = S(t).movKins.smJsVel(S(t).movKins.pullStart:S(t).movKins.pullStop);
        S(t).smJsVelStartMaxV = S(t).movKins.smJsVel(S(t).movKins.pullStart:S(t).movKins.pullMaxVelI);
        
        tempMaxVelAlignedVel = S(t).smJsVelFull(S(t).movKins.pullMaxVelI-min(100,S(t).movKins.pullMaxVelI)+1:S(t).movKins.pullMaxVelI);
        smJsVelMatMaxVelAlign(spTrCnt,100-length(tempMaxVelAlignedVel)+1:100)=-tempMaxVelAlignedVel;
        
        if isequal(S(t).pull_torque == pullTqs,[true,false])
            plot3(1:length(S(t).smJsVelStartStop),ones(1,length(S(t).smJsVelStartStop))*t,-S(t).smJsVelStartStop,'Color',tempColor)
        else
            plot3(1:length(S(t).smJsVelStartStop),-ones(1,length(S(t).smJsVelStartStop))*t,-S(t).smJsVelStartStop,'Color',tempColor)
        end
        
    end
end
hold off;

figure; hold on
surf(smJsVelMatMaxVelAlign)
caxis([0 400])
plot3(ones(1,size(smJsVelMatMaxVelAlign,1))*100,1:size(smJsVelMatMaxVelAlign,1),600+[S(spIdx&actIdx).pull_torque]) % indicate pull torque
title('Velocity aligned to max vel')
xlabel('Time (ms)')
ylabel('Trial#')
zlabel('Vel (mm/s)')
shading interp
grid on
hold off

% Acceleration 3-d plot
figure; hold on;
for t = 1:length(S)
    if spIdx(t)&&actIdx(t) % successful pull trials
        tempColor = pullTqsCmap(S(t).pull_torque == pullTqs,:);
        S(t).smJsAclFull = S(t).movKins.smJsAcl(1:S(t).movKins.pullStop);
        S(t).smJsAclStartStop = S(t).movKins.smJsAcl(S(t).movKins.pullStart:S(t).movKins.pullStop);
        S(t).smJsAclStartMaxV = S(t).movKins.smJsAcl(S(t).movKins.pullStart:S(t).movKins.pullMaxVelI);
        
        tempMaxVelAlignedAccel = S(t).smJsAclFull(S(t).movKins.pullMaxVelI-min(100,S(t).movKins.pullMaxVelI)+1:S(t).movKins.pullMaxVelI);
        smJsAccelMatMaxVelAlign(spTrCnt,100-length(tempMaxVelAlignedVel)+1:100)=-tempMaxVelAlignedVel;
        
        if isequal(S(t).pull_torque == pullTqs,[true,false])
            plot3(1:length(S(t).smJsAclStartStop),ones(1,length(S(t).smJsAclStartStop))*t,-S(t).smJsAclStartStop,'Color',tempColor)
        else
            plot3(1:length(S(t).smJsAclStartStop),-ones(1,length(S(t).smJsAclStartStop))*t,-S(t).smJsAclStartStop,'Color',tempColor)
        end
        
    end
end
hold off;
grid on

%%
jsTime1k_K = S; 
save('jsTime1k_kinematics','jsTime1k_K')

end

