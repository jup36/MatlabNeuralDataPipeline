function [ContData] = TNC_MoverBehaviorExtract(filenamestr,targetName,dataRate,chan)

%% FILE NAMES
    
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  

    dispOn = 0;
    
%% LOAD EVENT DATA
    
    filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'nev'];
    dataEvents = openNEV(filenamestrE,'read','nosave','nomat');

    rewIndsTmp = find(dataEvents.Data.Spikes.Electrode==chan.rew); % reward
        % eliminate any pulses that were detected twice
        rewTsTmp = round(dataEvents.Data.Spikes.Timestamps(rewIndsTmp)./30);
        validRewInds = find(diff(rewTsTmp)>200);
        % write valid data to structure
        ContData.behavior.rewardInds = rewTsTmp(validRewInds);

    thrIndsTmp = find(dataEvents.Data.Spikes.Electrode==chan.thr); % threshold crossing
        % eliminate any pulses that were detected twice
        thrTsTmp = round(dataEvents.Data.Spikes.Timestamps(thrIndsTmp)./30);
        validThrInds = find(diff(thrTsTmp)>10);
        % write valid data to structure
        ContData.behavior.threshInds = thrTsTmp(validThrInds);

disp(['Number of rewards delivered: ' num2str(numel(validRewInds)) '(' num2str(numel(rewTsTmp)) ' inds)']);
disp(' ');disp(' '); 
        
    lickInds = find(dataEvents.Data.Spikes.Electrode==chan.lick); % licking
        ContData.behavior.lickInds = round(dataEvents.Data.Spikes.Timestamps(lickInds)./30);
        
%% EXTRACT BLOCK STRUCTURE
 
    % work backwards from last trial
    numTrials = numel(ContData.behavior.threshInds);

    for p=numTrials:-1:1
        ContData.behavior.block(p) = 7-floor(p./15);
    end

%% LOAD CONTINUOUS LEVER DATA
    
switch dataRate
    
    case 'ns4'    
        filenamestrE = [filenamestr(1,1:length(filenamestr)-3) dataRate];
        Ns4DATA = openNSx('report','read',filenamestrE);
        clear leverData sLeverData tmpLeverData
        
        xChan = find(Ns4DATA.MetaTags.ChannelID==chan.x);
        yChan = find(Ns4DATA.MetaTags.ChannelID==chan.y);
        lChan = find(Ns4DATA.MetaTags.ChannelID==chan.lick);
        
        leverData(1,:) = decimate(Ns4DATA.Data(xChan,:),10);
        leverData(2,:) = decimate(Ns4DATA.Data(yChan,:),10);
        rawLick = decimate(Ns4DATA.Data(lChan,:),10);

    case 'ns3'
        filenamestrE = [filenamestr(1,1:length(filenamestr)-3) dataRate];
        Ns4DATA = openNSx('report','read',filenamestrE);
        clear leverData sLeverData tmpLeverData
        leverData(1,:) = decimate(Ns4DATA.Data(chan.x,:),2);
        leverData(2,:) = decimate(Ns4DATA.Data(chan.y,:),2);
        rawLick = decimate(Ns4DATA.Data(3,:),10);
        
        
    case 'ns2'
        filenamestrE = [filenamestr(1,1:length(filenamestr)-3) dataRate];
        Ns2DATA = openNSx('report','read',filenamestrE);
        clear leverData sLeverData tmpLeverData
        xChan = find(Ns2DATA.MetaTags.ChannelID==chan.x);
        yChan = find(Ns2DATA.MetaTags.ChannelID==chan.y);
        lChan = find(Ns2DATA.MetaTags.ChannelID==chan.lick);
        leverData(1,:)  = Ns2DATA.Data(xChan,:);
        leverData(2,:)  = Ns2DATA.Data(yChan,:);
        rawLick         = Ns2DATA.Data(lChan,:);
        
end
    
    sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);
    sLeverData(2,:) = sgolayfilt(leverData(2,:),9,101);

    ContData.behavior.sLeverData = sLeverData;

    difLick = diff(rawLick);
    evLick=zeros(1,numel(rawLick));

    ContData.behavior.evLick    = evLick;
    ContData.behavior.rawLick   = rawLick;

%% EXTRACT VELOCITY DATA FROM CONT LEVER DATA

    numSamples = size(ContData.behavior.sLeverData,2);
    tmpLeverData(1,:) = sgolayfilt(ContData.behavior.sLeverData(1,:),3,151);
    tmpLeverData(2,:) = sgolayfilt(ContData.behavior.sLeverData(2,:),3,151);

    sLeverV = zeros(1,numSamples);

    disp(' ');disp(' ');disp('Extracting velocity...');

    dX = diff(tmpLeverData(1,:));
    dY = diff(tmpLeverData(2,:));

    sLeverV = sqrt( dX.^2 + dY.^2 );

    disp(' ');disp(' Complete. ');disp(' ');

    ContData.behavior.sLeverV = sLeverV;
    ContData.behavior.sLeverVm = sgolayfilt(ContData.behavior.sLeverV,3,501);
    
%     figure(1); plot(sLeverV);
    
    clear sLeverV;

%% FIND MOVEMENTS
    method = 'vel';
    numC = 5;
    clear progSt* reach

    currX       = ContData.behavior.sLeverData(1,:);
    currY       = ContData.behavior.sLeverData(2,:);
    currV       = ContData.behavior.sLeverVm; 

    pre     = 10;
    post    = 10;       
    minSpace = 250;
    count = 1;

    % threshold the velocities
    switch dataRate
        case 'ns4'
            allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>0.35);
        case 'ns2'
            allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>1.5);
    end
    
    if numel(allValidSamps)>0

        switch method

            case 'vel'

    %             progStartTMP(count,1)  = currStamps(allValidSamps(1));
                progStartTMP(count,2)  = currX(allValidSamps(1));
                progStartTMP(count,3)  = currY(allValidSamps(1));
                progStartTMP(count,4)  = allValidSamps(1);

                for j=2:numel(allValidSamps)

                    if allValidSamps(j)>allValidSamps(j-1)+minSpace

                        postN = find(currV(allValidSamps(j-1):allValidSamps(j))<1,1,'first');
                        if numel(postN)<1
                            figure(3); clf;
                            plot(progStartTMP(count,4):allValidSamps(j),currV(progStartTMP(count,4):allValidSamps(j)),'k');
                            hold on; plot(progStartTMP(count,4),currV(progStartTMP(count,4)),'bo',allValidSamps(j-1),currV(allValidSamps(j-1)),'ro');
                            pause(0.1);

                            disp('Cannot find stop');                        
                            postN=post;
                        end
    %                     progStopTMP(count,1)   = currStamps(allValidSamps(j-1)+post);
                        progStopTMP(count,2)   = currX(allValidSamps(j-1)+postN);
                        progStopTMP(count,3)   = currY(allValidSamps(j-1)+postN);
                        progStopTMP(count,4)   = allValidSamps(j-1)+postN;
                        count                  = count+1;

                        preN = find(currV(allValidSamps(j-1) : allValidSamps(j))<1,1,'last');
                        if numel(preN)<1
                            disp('Cannot find start');
    %                     progStartTMP(count,1)  = currStamps(allValidSamps(j)-pre);
                            progStartTMP(count,2)  = currX(allValidSamps(j)-pre);
                            progStartTMP(count,3)  = currY(allValidSamps(j)-pre);
                            progStartTMP(count,4)  = allValidSamps(j)-pre;                    
                        else
    %                     progStartTMP(count,1)  = currStamps(allValidSamps(j-1)+preN);
                            progStartTMP(count,2)  = currX(allValidSamps(j-1)+preN);
                            progStartTMP(count,3)  = currY(allValidSamps(j-1)+preN);
                            progStartTMP(count,4)  = allValidSamps(j-1)+preN;
                        end
                    end

                    if j==numel(allValidSamps)
    %                     post = find(currV(allValidSamps(j):allValidSamps(j)+minSpace)<0.5,1,'first');
    %                     progStopTMP(count,1)   = currStamps(allValidSamps(j)+post);
                        progStopTMP(count,2)   = currX(allValidSamps(j)+post);
                        progStopTMP(count,3)   = currY(allValidSamps(j)+post);
                        progStopTMP(count,4)   = allValidSamps(j)+post;
                    end

                end

                count = 1;
                for k = 1:size(progStartTMP,1)

                    if k==1
                        reach.init = 1;
                    end

                    % reaches must be at least 50 ms long
                    if progStopTMP(k,4)-progStartTMP(k,4)>=90 & progStartTMP(k,4)>minSpace

                        trajAngle   = atan2(progStopTMP(k,3)-progStartTMP(k,3),progStopTMP(k,2)-progStartTMP(k,2));

                        if (pdist2([progStopTMP(k,2),progStopTMP(k,3)],[mean(currX),mean(currY)]) > pdist2([progStartTMP(k,2),progStartTMP(k,3)],[mean(currX),mean(currY)]))
                            reach.out(count) = 1;
                        else
                            reach.out(count) = 0;
                        end
                        velTraj = ContData.behavior.sLeverV(progStartTMP(k,4) : progStopTMP(k,4));
                        xVals = ContData.behavior.sLeverData(1,progStartTMP(k,4) : progStopTMP(k,4));
                        yVals = ContData.behavior.sLeverData(2,progStartTMP(k,4) : progStopTMP(k,4));

                        reach.start(count,:)  = progStartTMP(k,:);
                        reach.stop(count,:)   = progStopTMP(k,:);
                        reach.angle(count,1)  = trajAngle;
                        reach.dist(count,1)   = trapz(velTraj);
                        reach.dist(count,2)   = pdist2(progStartTMP(k,2:3) , progStopTMP(k,2:3));

                        tmp = findpeaks(velTraj);
                        reach.numpks(count,1) = numel(tmp.loc);
                        reach.dur(count,1)    = progStopTMP(k,4) - progStartTMP(k,4);
                        reach.vel(count,1)   = max(velTraj);
                        reach.vel(count,2)   = trapz(velTraj) ./ reach.dur(count,1);
                        reach.vel(count,3)   = var(velTraj);
                        reach.vel(count,4)   = find(velTraj==max(velTraj),1);

                        reach.acc(count,1)   = max(diff(velTraj));
                        reach.acc(count,2)   = mean(diff(velTraj));
                        reach.acc(count,3)   = max(diff(velTraj(1:90))); % max in first 90 ms of movement

                        reach.tort(count,1)  = reach.dist(count,1) ./ pdist2([progStopTMP(k,2),progStopTMP(k,3)],[progStartTMP(k,2),progStartTMP(k,3)]);

                        % find max displacement of the reach
                        xVals = xVals - xVals(1);
                        yVals = yVals - yVals(1);
                        
                        
                        
                        valThrInd = find(ContData.behavior.threshInds>reach.start(count,4) & ContData.behavior.threshInds<reach.stop(count,4));
                        if numel(valThrInd)>0
                            reach.rewarded(count)= 1;
                        else
                            reach.rewarded(count)= 0;                       
                        end

                        count                 = count+1;

                    end            
                end


        end

        disp(['Num valid reaches: ' num2str(count-1)]);

        level = ones(1,size(progStartTMP,1)).*0.35;
        figure(3); clf;
        plot(ContData.behavior.sLeverV,'k'); hold on;
        plot(progStartTMP(:,4),level,'r^');
        plot(progStopTMP(:,4),level,'bo');

        [vals,inds] = sort(reach.dur,'ascend');    

        reach.numReaches = size(reach.start,1);

        dims = floor(sqrt(reach.numReaches));
        scaler = 3000;

        if dispOn
            
            figure(5); clf;
        
            for p = 1:1:dims.^2

                l = inds(p);

                offset = l.*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
                [x,y] = ind2sub(dims,p);

                figure(5);
                plot((x.*scaler),(y.*scaler),'k.','MarkerSize',8);
                hold on;
                plot(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4))+(x.*scaler),currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4))+(y.*scaler) , 'k','LineWidth',1,'Color',[p/reach.numReaches 0.67-(p/reach.numReaches)*0.67 1-(p/reach.numReaches)]);
                if i==dims.^2
                    axis tight; axis off;
                end

                startInd = round(reach.start(l,4)./2);
                stopInd = round(reach.stop(l,4)./2);

            end
        end
        
        if reach.numReaches>5
                figure(7); clf;
                subplot(311);
                p = polyfit(reach.dist(:,1),reach.vel(:,1),1);
                plot(reach.dist,reach.vel(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Peak Velocity');

                subplot(312);
                p = polyfit(reach.dist(:,1),reach.dur,2);
                plot(reach.dist,reach.dur,'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Movement Duration (ms)');

                subplot(313);
                p = polyfit(reach.dist(:,1),reach.acc(:,3),1);
                plot(reach.dist,reach.acc(:,3),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Initial acceleration');

                figure(8); clf;
                rose(reach.angle);

                drawnow;
        end
        
    else
        
        reach.numReaches=0;
    
    end
    

    
%% DIMENSION REDUCTION OF REACH PARAMETERS

if reach.numReaches > 5
    clear matFor*

    matForDR = [reach.vel reach.acc reach.tort reach.dur reach.dist];

    for j=1:size(matForDR,2)
        matForCov(:,j) = (matForDR(:,j) - mean(matForDR(:,j))) ./ std(matForDR(:,j)) ;
    end

    figure(11);
    imagesc(corr(matForCov),[0 1]);
    map = colormap(TNC_CreateRBColormap(1024,'cpb'));
%     colormap(1-map);

    [eV,eD] = eig(cov(matForCov));
    figure(12); plot(cumsum(diag(eD))./sum(diag(eD)),'ko-');

    for m=1:size(matForDR,1)
        reach.pca(1,m) = dot(matForCov(m,:),eV(:,size(matForCov,2))');
        reach.pca(2,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-1)');
        reach.pca(3,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-2)');
    end
end

%% CALCULATE META-REACH PARAMETERS
    if reach.numReaches > 5
        [winReach] = TNC_ReachVigorWindow(reach.pca(1,:),reach.numReaches,9); % capture fluctuations along the maximally variant dimension
        ContData.behavior.winReach = winReach;
    end
    
%% WAS THE REACH REWARDED?
if reach.numReaches > 5
    for p=1:reach.numReaches

        tmp = find(ContData.behavior.threshInds>reach.start(p,4) & ContData.behavior.threshInds<reach.stop(p,4));
        if numel(tmp)>0
            reach.rewarded(p)=1;
        else
            reach.rewarded(p)=0;
        end
        
    end
end    
%% WRITE THE REACH STRUCTURE OF THE CONTDATA STRUCTURE
    ContData.behavior.reach = reach;
    
%% COMPLETED ALL ANALYSIS. SAVE AND START OVER

    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Completed file: ' filenamestr(1,1:length(filenamestr)-3)]);
    
    % save the data from this file
    disp(['saved as ' targetName '_bh.mat']);
    save([targetName '_bh.mat'],'ContData');
    
    disp('%-------------------------------------------------------------------');
    disp(' ');
    
