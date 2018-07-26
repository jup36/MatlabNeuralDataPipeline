
fileList = dir('*.ns5');
numFiles = size(fileList,1);

for q = 1:numFiles

    filenamestr = fileList(q).name;
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  
    
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % LOAD EVENT DATA
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    
    filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'nev']
    dataEvents = openNEV(filenamestrE,'read','nosave','nomat');

    nevRes = dataEvents.MetaTags.SampleRes;

    rewardInds = find(dataEvents.Data.Spikes.Electrode==142); % reward
    ContData.behavior.rewardInds = round(dataEvents.Data.Spikes.Timestamps(rewardInds)./30);

    threshInds = find(dataEvents.Data.Spikes.Electrode==143); % threshold crossing
    ContData.behavior.threshInds = round(dataEvents.Data.Spikes.Timestamps(threshInds)./30);

    lickInds = find(dataEvents.Data.Spikes.Electrode==139); % licking
    ContData.behavior.lickInds = round(dataEvents.Data.Spikes.Timestamps(lickInds)./30);
    
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % LOAD CONTINUOUS LEVER DATA
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------

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

%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % EXTRACT VELOCITY DATA FROM CONT LEVER DATA
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------

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

%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % FIND MOVEMENTS
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
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
                    
                    count                 = count+1;

                end            
            end


    end

    level = ones(1,size(progStartTMP,1));
    figure(3); clf;
    plot(ContData.behavior.sLeverV,'k'); hold on;
    plot(progStartTMP(:,4),level,'r^');
    plot(progStopTMP(:,4),level,'bo');

    numReaches = size(reach.start,1)
    
    figure(q+10); clf;
    subplot(411);
    p = polyfit(reach.dist,reach.vel(:,1),2);
    plot(reach.dist,reach.vel(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Peak Velocity');
    axis([0 2e4 0 50]);
    
    subplot(412);
    p = polyfit(reach.dist,reach.dur,2);
    plot(reach.dist,reach.dur,'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Movement Duration (ms)');
    axis([0 2e4 0 4e3]);
    
    subplot(413);
    p = polyfit(reach.dist,reach.acc(:,3),1);
    plot(reach.dist,reach.acc(:,3),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Initial acceleration');
    axis([0 2e4 0 0.5]);

    subplot(414);
    rose(reach.angle);
    
    reach.numReaches = numReaches;
    ContData.behavior.reach = reach;
    drawnow;
    
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % COMPLETED ALL ANALYSIS. SAVE AND START OVER
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
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