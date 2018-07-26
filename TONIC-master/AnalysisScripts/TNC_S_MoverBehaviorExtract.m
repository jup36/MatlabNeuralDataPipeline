% TNC_S_MoverBehaviorExtract    
fileList = dir('*.ns5');
numFiles = size(fileList,1)

for q = 1:numFiles
    
%% FILE NAMES
    filenamestr = fileList(q).name;
    
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  

%% LOAD EVENT DATA
    
    filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'nev']
    dataEvents = openNEV(filenamestrE,'read','nosave','nomat');

    rewIndsTmp = find(dataEvents.Data.Spikes.Electrode==142); % reward
        % eliminate any pulses that were detected twice
        rewTsTmp = round(dataEvents.Data.Spikes.Timestamps(rewIndsTmp)./30);
        validRewInds = find(diff(rewTsTmp)>10);
        % write valid data to structure
        ContData.behavior.rewardInds = rewTsTmp(validRewInds);

    thrIndsTmp = find(dataEvents.Data.Spikes.Electrode==143); % threshold crossing
        % eliminate any pulses that were detected twice
        thrTsTmp = round(dataEvents.Data.Spikes.Timestamps(thrIndsTmp)./30);
        validThrInds = find(diff(thrTsTmp)>10);
        % write valid data to structure
        ContData.behavior.threshInds = thrTsTmp(validThrInds);

    lickInds = find(dataEvents.Data.Spikes.Electrode==139); % licking
        ContData.behavior.lickInds = round(dataEvents.Data.Spikes.Timestamps(lickInds)./30);
        
%% EXTRACT BLOCK STRUCTURE
 
    % work backwards from last trial
    numTrials = numel(ContData.behavior.threshInds)

    for p=numTrials:-1:1
        ContData.behavior.block(p) = 7-floor(p./15);
    end

%% LOAD CONTINUOUS LEVER DATA

    filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'ns4']
    Ns4DATA = openNSx('report','read',filenamestrE);
    clear leverData sLeverData tmpLeverData 

    leverData(1,:) = decimate(Ns4DATA.Data(1,:),10);
    leverData(2,:) = decimate(Ns4DATA.Data(2,:),10);

    sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);
    sLeverData(2,:) = sgolayfilt(leverData(2,:),9,101);

    ContData.behavior.sLeverData = sLeverData;

    if strcmp(filenamestr(1),'M')
        figure(55); clf;
    else
        figure(5); clf;
    end
    plotArray = 1:250000;
    plot3(sLeverData(1,plotArray),sLeverData(2,plotArray),plotArray,'k');

    rawLick = decimate(Ns4DATA.Data(3,:),10);
    difLick = diff(rawLick);
    evLick=zeros(1,numel(rawLick));

    ContData.behavior.evLick    = evLick;
    ContData.behavior.rawLick   = rawLick;

%% EXTRACT VELOCITY DATA FROM CONT LEVER DATA

    numSamples = size(ContData.behavior.sLeverData,2);
    tmpLeverData(1,:) = sgolayfilt(ContData.behavior.sLeverData(1,:),3,201);
    tmpLeverData(2,:) = sgolayfilt(ContData.behavior.sLeverData(2,:),3,201);

    sLeverV = zeros(1,numSamples);

    disp(' ');disp(' ');disp('Extracting velocity...');

    dX = diff(tmpLeverData(1,:));
    dY = diff(tmpLeverData(2,:));

    sLeverV = sqrt( dX.^2 + dY.^2 );

    disp(' ');disp(' Complete. ');disp(' ');

    ContData.behavior.sLeverV = sLeverV;
    ContData.behavior.sLeverVm = sgolayfilt(ContData.behavior.sLeverV,3,501);
    clear sLeverV;

%% FIND MOVEMENTS
    method = 'vel';
    numC = 5;
    clear progSt* reach

    currX       = ContData.behavior.sLeverData(1,:);
    currY       = ContData.behavior.sLeverData(2,:);
    currV       = ContData.behavior.sLeverVm; 
    

    switch method

        case 'vel'

            pre     = 10;
            post    = 10;       
            minSpace = 250;
            count = 1;
            
            % threshold the velocities
            allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>2);

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
                    reach.dist(count,:)   = trapz(velTraj);
                    
                    tmp = findpeaks(velTraj);
                    reach.numpks(count,1) = numel(tmp.loc);
                    reach.dur(count,1)    = progStopTMP(k,4) - progStartTMP(k,4);
                    reach.vel(count,1)   = max(velTraj);
                    reach.vel(count,2)   = trapz(velTraj) ./ reach.dur(count,1);
                    reach.vel(count,3)   = var(velTraj);

                    reach.acc(count,1)   = max(diff(velTraj));
                    reach.acc(count,2)   = mean(diff(velTraj));
                    reach.acc(count,3)   = max(diff(velTraj(1:90))); % max in first 90 ms of movement
                    
                    reach.tort(count,1)  = reach.dist(count,:) ./ pdist2([progStopTMP(k,2),progStopTMP(k,3)],[progStartTMP(k,2),progStartTMP(k,3)]);
                    
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

    level = ones(1,size(progStartTMP,1));
    figure(3); clf;
    plot(ContData.behavior.sLeverV,'k'); hold on;
    plot(progStartTMP(:,4),level,'r^');
    plot(progStopTMP(:,4),level,'bo');

%     [vals,inds] = sort(reach.dist+(reach.maxV(:,3).*1e3),'ascend');        
%     [vals,inds] = sort(reach.maxV(:,3),'ascend');    
    [vals,inds] = sort(reach.dur,'ascend');    

    numReaches = size(reach.start,1)

    figure(5); clf;
    figure(7); clf;

    dims = floor(sqrt(numReaches));
    scaler = 3000;
    
    for p = 1:1:dims.^2
        
        l = inds(p);
        
    %   offset = reach.start(l,4).*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
    %   plot3(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4)), currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4)) , offset , 'r','LineWidth',1);
        
        offset = l.*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
        [x,y] = ind2sub(dims,p);
        
        figure(5);
        plot((x.*scaler),(y.*scaler),'k.','MarkerSize',8);
        hold on;
        plot(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4))+(x.*scaler),currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4))+(y.*scaler) , 'k','LineWidth',1,'Color',[p/numReaches 0.67-(p/numReaches)*0.67 1-(p/numReaches)]);
        if i==dims.^2
            axis tight; axis off;
        end
                
%         figure(7);
        startInd = round(reach.start(l,4)./2);
        stopInd = round(reach.stop(l,4)./2);
        
%         meanField = sqrt(ContData.shank(2).contData.betaP(startInd-(minSpace/2):startInd+600));
%         for d=3:5
%             meanField = meanField + sqrt(ContData.shank(d).contData.betaP(startInd-(minSpace/2):startInd+600));
%         end
%         meanField = meanField ./ 4;
            
%         plot((x.*scaler)+([startInd-(minSpace/2):startInd+600]-startInd).*2 , meanField + (y.*scaler./20),'k','LineWidth',1,'Color',[p/numReaches 0.67-(p/numReaches)*0.67 1-(p/numReaches)]);
%         hold on;
% %         plot((x.*scaler)+([startInd-(minSpace/2):startInd+600]-startInd).*2 , ContData.shank(4).contData.beta(startInd-(minSpace/2):startInd+600) + (y.*scaler./20),'k--','LineWidth',1,'Color',[p/numReaches 0.67-(p/numReaches)*0.67 1-(p/numReaches)]);
%         plot((x.*scaler)+[reach.start(l,4)-minSpace:reach.start(l,4)+1200]-reach.start(l,4) , ContData.behavior.sLeverV(reach.start(l,4)-minSpace:reach.start(l,4)+1200).*5 + (y.*scaler./20),'k','LineWidth',2);
%         plot((x.*scaler),(y.*scaler./20),'k.','MarkerSize',8);
%         if i==dims.^2
%             axis tight; axis off;
%         end
    end
    
    figure(7); clf;
    subplot(311);
    p = polyfit(reach.dist,reach.vel(:,1),1);
    plot(reach.dist,reach.vel(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Peak Velocity');
    
    subplot(312);
    p = polyfit(reach.dist,reach.dur,2);
    plot(reach.dist,reach.dur,'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Movement Duration (ms)');
    
    subplot(313);
    p = polyfit(reach.dist,reach.acc(:,3),1);
    plot(reach.dist,reach.acc(:,3),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Initial acceleration');
    
    figure(8); clf;
    rose(reach.angle);
    
    reach.numReaches = numReaches;
            
    drawnow;
    
%% DIMENSION REDUCTION OF REACH PARAMETERS

    clear matFor*

    matForDR = [reach.vel reach.acc reach.tort reach.dur reach.dist];

    for j=1:size(matForDR,2)
        matForCov(:,j) = (matForDR(:,j) - mean(matForDR(:,j))) ./ std(matForDR(:,j)) ;
    end

    figure(11);
    imagesc(corr(matForCov),[0 1]);
    map = colormap(pink);
    colormap(1-map);

    [eV,eD] = eig(cov(matForCov));
    figure(12); plot(cumsum(diag(eD))./sum(diag(eD)),'ko-');

    for m=1:size(matForDR,1)
        reach.pca(1,m) = dot(matForCov(m,:),eV(:,size(matForCov,2))');
        reach.pca(2,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-1)');
        reach.pca(3,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-2)');
    end

%% CALCULATE META-REACH PARAMETERS
    if numReaches>0
        [winReach] = TNC_ReachVigorWindow(reach.pca(1,:),reach.numReaches,9); % capture fluctuations along the maximally variant dimension
        ContData.behaviora.winReach = winReach;
    end
    
%% WAS THE REACH REWARDED?
    for p=1:reach.numReaches

        tmp = find(ContData.behavior.threshInds>reach.start(p,4) & ContData.behavior.threshInds<reach.stop(p,4));
        if numel(tmp)>0
            reach.rewarded(p)=1;
        else
            reach.rewarded(p)=0;
        end
        
    end
    
%% WRITE THE REACH STRUCTURE OF THE CONTDATA STRUCTURE
    ContData.behavior.reach = reach;
    
%% COMPLETED ALL ANALYSIS. SAVE AND START OVER

    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Completed file: ' filenamestr(1,1:length(filenamestr)-3)]);

    ContData
    
    % save the data from this file
    disp(['save ' filenamestr(1,1:length(filenamestr)-3) 'mat ContData']);
    eval(['save ' filenamestr(1,1:length(filenamestr)-3) 'mat ContData']);
    
    disp('%-------------------------------------------------------------------');
    disp(' ');
    
    % clear memory and start again
    clear ContData prog* sLev* lowB* Ns* dat* reach
    
end