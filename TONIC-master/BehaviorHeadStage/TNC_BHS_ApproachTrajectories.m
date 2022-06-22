%% Script for analysis of behavior headstage data: TNC_BHS_ApproachBehavior


ns3File = 'F06-120222trace-t500-005-01test-02.ns3'
nevFile = 'F06-120222trace-t500-005-01test-02.nev'
csvFile = 'F062012222192GA.csv'
nexFile = 'F06-120222trace-t500-005-01test-02.nex'

trajWin = [10 400];

%% LOAD THE CSV DATA
    
    % Number of samples in running median
    runMedWin   = 5;

    fileNameBase = csvFile(1:numel(csvFile)-4);
    year = fileNameBase(strfind(fileNameBase,'201'):strfind(fileNameBase,'201')+3);
    month = fileNameBase(strfind(fileNameBase,'201')+4:strfind(fileNameBase,'201')+4);
    if month==1
        month = fileNameBase(strfind(fileNameBase,'201')+4:strfind(fileNameBase,'201')+5);
        day = fileNameBase(strfind(fileNameBase,'201')+6:strfind(fileNameBase,'201')+7);
    else
        month = ['0' fileNameBase(strfind(fileNameBase,'201')+4:strfind(fileNameBase,'201')+4)];
        day = fileNameBase(strfind(fileNameBase,'201')+5:strfind(fileNameBase,'201')+6);
    end

    disp(' ');
    disp(['Filename: ' fileNameBase]);
    disp(['Date created: ' year '/' month '/' day]);    
    disp(' ');

    data.creation.year  = year;
    data.creation.month = month;
    data.creation.day   = day;
    
    
    % Script to load csv behavior data
    clear disp* trial* 
 
    tmpData = dlmread(csvFile,',',2,0);
    numSamps = size(tmpData,1);
    
    %___________________________________________________________
    % TIMESTAMPS (MS)
    data.tStamps = tmpData(:,1);
    
    %___________________________________________________________
    % SMOOTH AND NICE UP THE ACCEL DATA
    data.accel.x = tmpData(:,2);
    data.accel.y = tmpData(:,3);
    data.accel.z = tmpData(:,4);
    data.accel.mag=sqrt(data.accel.x.^2+data.accel.y.^2+data.accel.z.^2);
    
    %___________________________________________________________
    % SMOOTH AND NICE UP THE GYRO DATA
    data.gyro.x = tmpData(:,5);
    data.gyro.y = tmpData(:,6);
    data.gyro.z = tmpData(:,7);
    data.gyro.mag = sqrt(data.gyro.x.^2+data.gyro.y.^2+data.gyro.z.^2);
 
    %___________________________________________________________
    % DIGITAL SIGNALS FROM BEHAVIOR SYSTEM
    data.dig    = tmpData(:,8);
 
    %___________________________________________________________
    % SMOOTH AND NICE UP THE POSITION DATA
    data.posRaw.xA = tmpData(:,9);
    data.posRaw.yA = tmpData(:,10);
    data.posRaw.xB = tmpData(:,11);
    data.posRaw.yB = tmpData(:,12);
 
    %___________________________________________________________
    % POSITION DATA APPEARS TO BE OVERSAMPLED AND SO I NEED TO REMOVE ZEROS
    % AND THEN CREATE THE CONTINUOUS FUNCTION AGAIN WITH SMOOTHING
    data.posSmth.xA = sgolayfilt(data.posRaw.xA,11,31);
    data.posSmth.yA = sgolayfilt(data.posRaw.yA,11,31);
    data.posSmth.xB = sgolayfilt(data.posRaw.xB,11,31);
    data.posSmth.yB = sgolayfilt(data.posRaw.yB,11,31);
 
    data.posSmth.x  = (data.posSmth.xA+data.posSmth.xB)./2;
    data.posSmth.y  = (data.posSmth.yA+data.posSmth.yB)./2;    
    
    %___________________________________________________________
    % FURTHER SMOOTH POSITION DATA BY CREATING RUNNING MEDIANS    
    data.posSmth.Xrmed=medfilt1(data.posSmth.x,runMedWin);
    data.posSmth.Yrmed=medfilt1(data.posSmth.y,runMedWin);
            
    %___________________________________________________________
    % EUCLIDEAN DISTANCE WITH SMOOTHED POSITION DATA
    i=2:numel(data.posSmth.Xrmed);
    data.velocity.vel(i) = sqrt( (data.posSmth.Xrmed(i)-data.posSmth.Xrmed(i-1)).^2 + (data.posSmth.Yrmed(i)-data.posSmth.Yrmed(i-1)).^2 );
 
    %___________________________________________________________
    % FURTHER SMOOTH VELOCITY DATA 
    data.velocity.velRM     = abs(sgolayfilt(data.velocity.vel,5,11));
    data.velocity.velRML(i) = real(log(data.velocity.velRM(i)));
    
   %___________________________________________________________
   % TAKE CUMULATIVE SUMS OF DISTANCE TO GET DISTANCE FROM ORIGIN
    data.velocity.cumdistance = (cumsum(data.velocity.vel))';
    
    clear tmpData;

   %___________________________________________________________
   % FIND THE TRIAL STAMPS
    data.trialInds  = find(diff(data.dig)==1);
    data.numTrials  = numel(data.trialInds);
    
   %___________________________________________________________
   % COMPUTE THE PERI-TRIAL TRAJECTORIES 
    [trajectories] = TNC_AlignTraj2d(data.posSmth.x,data.posSmth.y,data.velocity.velRM,data.trialInds,trajWin);
%     [trajectoriesA] = TNC_AlignTraj2d(data.posSmth.xA,data.posSmth.yA,data.velocity.velRM,data.trialInds,trajWin);
%     [trajectoriesB] = TNC_AlignTraj2d(data.posSmth.xB,data.posSmth.yB,data.velocity.velRM,data.trialInds,trajWin);
    
sessData.data = data;

%% LOAD THE NEV/NEX DATA

    dataEvents = openNEV(nevFile,'read','nosave','nomat');

    csInds = find(dataEvents.Data.Spikes.Electrode==137 & dataEvents.Data.Spikes.Unit==1); % cs onset
        ContData.behavior.csInds = double(round(dataEvents.Data.Spikes.Timestamps(csInds)./30));

    lickInds = find(dataEvents.Data.Spikes.Electrode==138); % licking
        ContData.behavior.lickInds = double(round(dataEvents.Data.Spikes.Timestamps(lickInds)./30));

        
%% Check that position data makes sense with lickinds

    for j=1:numel(ContData.behavior.lickInds)
        % for every lick time find the correspoding index in the 
        tmp = find((data.tStamps-data.tStamps(1))>ContData.behavior.lickInds(j),1);
        ContData.behavior.lickPosInd(j) = tmp;
    end
    
    figure(4);
    plot(data.posSmth.xB,data.posSmth.yB,'k',data.posSmth.xB(ContData.behavior.lickPosInd),data.posSmth.yB(ContData.behavior.lickPosInd),'bo');
%     plot(data.tStamps,data.posSmth.y,'k',data.tStamps(ContData.behavior.lickPosInd),data.posSmth.y(ContData.behavior.lickPosInd),'bo');
        
%% GET FIRST LICK LATENCY FOR EACH TRIAL

    for i=1:numel(ContData.behavior.csInds)
        
        thisTrialLicks = double( ContData.behavior.lickInds( find(ContData.behavior.lickInds > ContData.behavior.csInds(i) , 1) ));
        if numel(thisTrialLicks)>0
            ContData.behavior.firstLick(i) = double(thisTrialLicks(1)-ContData.behavior.csInds(i));
        else
            disp(['no FL on trial ' num2str(i)]);
            ContData.behavior.firstLick(i) = 10000;
        end        
    end
    
    [vals,inds] = sort(ContData.behavior.firstLick);

%% DISPLAY APPROACH SORTED BY LATENCY
count = 0;
colInd = TNC_CreateRBColormap(numel(inds),'cpb');
figure(1); clf;
for j=numel(inds)-1:-1:1

    i = inds(j);

    portSamp = find( data.tStamps > (data.tStamps(data.trialInds(i)) + vals(j)) , 1);

    % check for a bogus trial where position tracking was shit
    tmp = find((data.tStamps-data.tStamps(1))>ContData.behavior.csInds(i)+ContData.behavior.firstLick(i),1);
    if data.posSmth.y(tmp) < 150 % legit trial; if greater it means that is some b.s. jumping around of position data
        count = count+1;
        numSamps = portSamp - data.trialInds(i);
        if numSamps>400
            numSamps=400;
        end
        figure(1);
%             figure(1); subplot(4,5,j);
        plot(trajectories.x(i,1),trajectories.y(i,1),'o','MarkerSize',10,'LineWidth',2,'Color',colInd(i,:)); hold on;
        plot(trajectories.x(i,1:numSamps),trajectories.y(i,1:numSamps),'-','LineWidth',1,'Color',colInd(i,:));
%             plot(trajectoriesB.x(i,1:numSamps),trajectoriesB.y(i,1:numSamps),'-','LineWidth',1,'Color',[0.5 0 0]);
        plot([min(min(trajectories.x)) max(max(trajectories.x)) max(max(trajectories.x)) min(min(trajectories.x)) min(min(trajectories.x))],[min(min(trajectories.y)) min(min(trajectories.y)) max(max(trajectories.y)) max(max(trajectories.y)) min(min(trajectories.y))],'-','LineWidth',1,'Color',[0.5 0.5 0.5]);
        axis off; axis([0 400 0 400]); 
%         title(num2str(ContData.behavior.firstLick(i)));
    end        
end

disp([num2str(count) ' of ' num2str(numel(inds)) ' trials had valid position data']);
    
sessData.trajectories = trajectories;

%% LOOK AT GYRO AND ACCEL SIGNALS SORTED BY THE LATENCY TO THE FIRST LICK
    
[gyroAlign] = TNC_AlignTraj2d(data.gyro.x,data.gyro.y,data.gyro.mag,data.trialInds,trajWin);
[accelAlign] = TNC_AlignTraj2d(data.accel.x,data.accel.y,data.accel.z,data.trialInds,trajWin);

aScale = 150;
gScale = 2000;
numPnts = 50;

figure(2); clf;
for j=1:numel(inds)-1

    i = inds(j);
% 
% figure(2); subplot(231);
% plot(1:numPnts,accelAlign.x(i,1:numPnts)-accelAlign.x(i,1)+(j.*aScale),'k'); hold on;
% figure(2); subplot(232);
% plot(1:numPnts,accelAlign.y(i,1:numPnts)-accelAlign.y(i,1)+(j.*aScale),'k'); hold on;
% figure(2); subplot(233);
figure(2); subplot(121);
plot(1:numPnts,abs(accelAlign.v(i,1:numPnts)-accelAlign.v(i,1))+(j.*aScale),'k',10,accelAlign.v(i,10)-accelAlign.v(i,1)+(j.*aScale),'ro'); hold on;

% figure(2); subplot(234);
% plot(1:numPnts,gyroAlign.x(i,1:numPnts)-gyroAlign.x(i,1)+(j.*gScale),'k'); hold on;
% figure(2); subplot(235);
% plot(1:numPnts,gyroAlign.y(i,1:numPnts)-gyroAlign.y(i,1)+(j.*gScale),'k'); hold on;
figure(2); subplot(122);
plot(1:numPnts,gyroAlign.v(i,1:numPnts)-gyroAlign.v(i,1)+(j.*gScale),'k',10,gyroAlign.v(i,10)-gyroAlign.v(i,1)+(j.*gScale),'ro'); hold on;

end

figure(3); clf
shadedErrorBar([1:numPnts].*33,mean(gyroAlign.v(:,1:numPnts),1)./10,std(gyroAlign.v(:,1:numPnts)./sqrt(size(gyroAlign.v,1))./10,[],1),'k')
hold on;
shadedErrorBar([1:numPnts].*33,mean(accelAlign.v(:,1:numPnts),1),std(accelAlign.v(:,1:numPnts)./sqrt(size(accelAlign.v,1)),[],1),'r')
shadedErrorBar([1:numPnts].*33,mean(accelAlign.x(:,1:numPnts),1),std(accelAlign.x(:,1:numPnts)./sqrt(size(accelAlign.x,1)),[],1),'c')
shadedErrorBar([1:numPnts].*33,mean(accelAlign.y(:,1:numPnts),1),std(accelAlign.y(:,1:numPnts)./sqrt(size(accelAlign.y,1)),[],1),'g')
plot([10 10].*33,[-50 250],'k--');
plot([14 14].*33,[-50 250],'r--');

sessData.accelAlign = accelAlign;
sessData.gyroAlign = gyroAlign;

%% COMPARE BEHAVIOR AS A FUNCTION OF FL (not sure this makes sense)

% Analysis idea:
% if I want to show that FL does not predict approach i need the pairwise correlation in approach trajectory as a function of the difference in FL delay.
% for gyro/accel that is easy just take the first x milliseconds and compare. 
% for trajectories it is trickier, but I think I can just resample the trajectory with interp to be the same number of points and then do pairwise correlation. 
count = 0;
for i=2:numel(inds)
    for j=1:i-1
        
        count = count+1;
        
        % difference in approach time
        sessData.appDiff(count) = abs(ContData.behavior.firstLick(i) - ContData.behavior.firstLick(j));
        
        % pairwise difference in gyro
        sessData.gyroCorr(count) = corr(gyroAlign.v(i,trajWin(1):numPnts)',gyroAlign.v(j,trajWin(1):numPnts)');
        
        % pairwise difference in accel
        sessData.accCorr(count) = corr2( [accelAlign.x(i,trajWin(1):numPnts)' accelAlign.y(i,trajWin(1):numPnts)' accelAlign.v(i,trajWin(1):numPnts)'] , ...
                                [accelAlign.x(j,trajWin(1):numPnts)' accelAlign.y(j,trajWin(1):numPnts)' accelAlign.v(j,trajWin(1):numPnts)'] ...
                                );
                            
        % pairwise correlation of trajectory
        sessData.trajCorr(count) = corr2([trajectories.x(i,trajWin(1):numPnts)',trajectories.y(i,trajWin(1):numPnts)'],[trajectories.x(j,trajWin(1):numPnts)',trajectories.y(j,trajWin(1):numPnts)']);

    end
end
    
figure(201); 

subplot(131);
semilogx(sessData.appDiff,sessData.gyroCorr,'k.'); axis([1 1e5 -1 1]); 
title(num2str(corr(double(sessData.appDiff'),double(sessData.gyroCorr'))));
xlabel('Difference in approach time (ms)'); ylabel('Correlation in gyro signal');

subplot(132);
semilogx(sessData.appDiff,sessData.accCorr,'k.'); axis([1 1e5 -1 1]);
title(num2str(corr(double(sessData.appDiff'),double(sessData.accCorr'))));
xlabel('Difference in approach time (ms)'); ylabel('Correlation in accelerometer signal');

subplot(133);
semilogx(sessData.appDiff,sessData.trajCorr,'k.'); axis([1 1e5 -1 1]);
title(num2str(corr(double(sessData.appDiff'),double(sessData.trajCorr'))));
xlabel('Difference in approach time (ms)'); ylabel('Correlation in trajectory');

%% save the output
newFile = [csvFile(1:numel(csvFile)-4) '_sd']
save(newFile,'sessData');