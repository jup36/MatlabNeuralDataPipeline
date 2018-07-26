function [reachAnalysis] = TNC_NewLeverTrajectoryAnalysis(ContData,dataType,figNum)
% FUNCTION DETAILS: Analysis of continuous position and velocity data acquired through the MOVER behavior program.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN (dudmanlab.org)
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% _________________________________________________________________________
%
% ContData is a structure produced by loading data through:
% 1) TNC_MoverBehaviorExtract for data acquired with Blackrock
% 2) for data acquired directly with the Processing code
%
% dataType is a string that can take one of two values
% 1) 'blackrock'
% 2) 'processing'
%
% figNum is the beginning figure number. To suppress plotting use figNum=0
%

%% Pre analysis set up...
disp(' ');
disp(' ');

if figNum==0
    displayOnline = 0;
else
    displayOnline=1;
end

switch dataType

    case 'blackrock'
        maxReach = size(ContData.behavior.reach.vel,1)
    case 'processing'
        maxReach = size(ContData.data.reach.vel,1)
        minThresh = min(ContData.data.trialParams(:,6));
        maxThresh = max(ContData.data.trialParams(:,6));
end

%% Per reach analysis
if displayOnline==1
    figure(figNum+10); 
    clf;
end

clear forFit;

for reachNum = 1:maxReach;
    
    switch dataType
        
        case 'blackrock'
            xVals = ContData.behavior.sLeverData(1, ContData.behavior.reach.start(reachNum,4):ContData.behavior.reach.stop(reachNum,4));
            yVals = ContData.behavior.sLeverData(2, ContData.behavior.reach.start(reachNum,4):ContData.behavior.reach.stop(reachNum,4));
            vVals = ContData.behavior.sLeverV(1, ContData.behavior.reach.start(reachNum,4):ContData.behavior.reach.stop(reachNum,4));
            startR = ContData.behavior.reach.start(reachNum,4);
            stopR = ContData.behavior.reach.stop(reachNum,4);
            xVals = xVals - xVals(1);
            yVals = yVals - yVals(1);

            trialThresh = reachNum;
            
            reachRange = [-6000 6000 -6000 6000];
            multiplier = 1;
            posMult = 1;
            currColor = [(reachNum./maxReach) 0.67.*(1-(reachNum./maxReach)) 1-(reachNum./maxReach)];

        case 'processing'    
            totalTrials = size(ContData.data.trialTimes,1);
            startR  = find(ContData.data.trajectories.contXY(:,1) >= ContData.data.reach.start(reachNum) , 1);
            stopR   = find(ContData.data.trajectories.contXY(:,1) >= ContData.data.reach.stop(reachNum) , 1);
            trialNum= find(ContData.data.trialTimes(:,2) >= ContData.data.reach.start(reachNum) , 1);
                if numel(trialNum)==1
                    trialThresh = ContData.data.trialParams(trialNum,6);
                else
                    trialThresh = ContData.data.trialParams(size(ContData.data.trialParams,1),6);
                end
%             yVals = ContData.data.trajectories.contXYsmooth(startR:stopR, 3);
%             xVals = ContData.data.trajectories.contXYsmooth(startR:stopR, 2);

            xVals = sgolayfilt(ContData.data.trajectories.contXY(startR:stopR, 2) , 5 , 11);
            yVals = sgolayfilt(ContData.data.trajectories.contXY(startR:stopR, 3) , 5 , 11);
            
            vVals = ContData.data.trajectories.sLeverVm(startR:stopR) .* 100;
            xVals = xVals - xVals(1);
            yVals = yVals - yVals(1);   
            
            reachRange = [-100 100 -100 100];            
            multiplier = 33;
            posMult = 100;
            
            currColor = [(trialThresh./maxThresh) 0.67.*(1-(trialThresh./maxThresh)) 1-(trialThresh./maxThresh)];
    end
    
    try
        k = convhull(xVals,yVals);
    catch ME
        k = 1:numel(xVals);
    end

    hullPntsX = xVals(k);
    hullPntsY = yVals(k);
    hullDists = sqrt( hullPntsX.^2 + hullPntsY.^2 );
    maxDisplace = find(hullDists == max(hullDists),1);
    maxDispTotalTime = find(xVals==hullPntsX(maxDisplace) & yVals==hullPntsY(maxDisplace));

        changeInX = diff(xVals(1:maxDispTotalTime)).^2;
        changeInY = diff(yVals(1:maxDispTotalTime)).^2;
        changeInPosition = sqrt(changeInX + changeInY);
            
        forFit(reachNum,1) = sum(changeInPosition).*posMult;
        forFit(reachNum,2) = max(hullDists).*posMult;
        forFit(reachNum,3) = maxDispTotalTime.*multiplier;
        forFit(reachNum,4) = trialThresh;
        forFit(reachNum,5) = polyarea( hullPntsX , hullPntsY );
        forFit(reachNum,6) = stopR - startR;
    
    if displayOnline==1
        
        figure(figNum+10); 

        
        subplot(3,3,[1 2 4 5 7 8]);
                scatter(xVals,yVals,50,vVals,'filled'); hold on; 
                colormap(TNC_CreateRBColormap(1024,'wred'))
                plot(xVals(1),yVals(1),'ro', 'MarkerSize',20, 'LineWidth', 3); 
                plot(xVals(k),yVals(k),'k-');
                plot(hullPntsX(maxDisplace),hullPntsY(maxDisplace),'bo','MarkerSize',20, 'LineWidth', 3); 
                axis(reachRange); 
                hold off;

        subplot(3,3,[3]); 
                plot([0 60],[0 60],'k--'); hold on;
                plot( max(hullDists).*posMult, sum(changeInPosition).*posMult, 'o' , 'Color' , currColor); hold on; 
                axis tight;


        subplot(3,3,[6]); % fitts law type plot (log( 1+ (amplitude/displacement) )  vs. duration )
                semilogy( forFit(reachNum,6) , forFit(reachNum,5) , 'o' , 'Color' , currColor ); hold on;
                if reachNum==maxReach
                    p0 = polyfit(forFit(:,6),forFit(:,5),1);
                    plot(min(forFit(:,6)):max(forFit(:,6)),polyval(p0,min(forFit(:,6)):max(forFit(:,6))),'r-');
                end
                axis tight;

         subplot(3,3,[9]); 
                plot([sum(changeInPosition).*posMult max(hullDists).*posMult],[maxDispTotalTime.*multiplier maxDispTotalTime.*multiplier],'ko-', 'Color' , currColor);
                hold on; 
                axis tight;

                if reachNum==maxReach
                    p1 = polyfit(forFit(:,1),forFit(:,2),2);
                    p2 = polyfit(forFit(:,1),forFit(:,3),2);
                    plot( min(forFit(:,1)):max(forFit(:,1)), polyval(p1,min(forFit(:,1)):max(forFit(:,1))) ,'r-', min(forFit(:,1)):max(forFit(:,1)), polyval(p2,min(forFit(:,1)):max(forFit(:,1))) ,'r--');
                end

    end

end

reachAnalysis.fitPerReach   = forFit;
reachAnalysis.maxReach      = maxReach;

disp('Completed the per reach analysis...');

%% Per reward analysis

if displayOnline==1
    figure(figNum); 
    clf;
end

clear dispData;

switch dataType

    case 'blackrock'
        spacing         = 0.05;
        numTrials       = numel(ContData.behavior.rewardInds);

        rewardTimes     = ContData.behavior.rewardInds;
        eventTimes      = ContData.behavior.threshInds;
        offset = 1;
        
        for i=1+offset:numel(rewardTimes)
    
            startR      = rewardTimes(i-1);
            stopR       = rewardTimes(i);
            (stopR-startR)
            eventTime   = find(eventTimes >= startR, 1) - startR;

            xVals       = ContData.behavior.sLeverData(1,startR:stopR);
            yVals       = ContData.behavior.sLeverData(2,startR:stopR);
            vVals       = ContData.behavior.sLeverV(1,startR:stopR);

            totalDisp                   = trapz(vVals);
            currColor                   = [(i./numel(rewardTimes)) 0.67.*(1-(i./numel(rewardTimes))) 1-(i./numel(rewardTimes))];

            dispData.disp(i-offset)     = totalDisp;
            dispData.rrate(i-offset)    = 1000 ./ double(stopR-startR);
            dispData.thresh(i-offset)   = 10;
            dispData.vel(i-offset)      = max(vVals);

            if displayOnline==1
                figure(figNum);
                subplot(1,5,1);
                plot(yVals,xVals+(i*500),'-','Color',currColor,'LineWidth',1); hold on;
                subplot(1,5,2:3);
                plot([1:numel(vVals)],vVals+(i*8),'-','Color',currColor,'LineWidth',1); hold on;
                subplot(1,5,4);
                semilogx(totalDisp,i*spacing,'-','Color',currColor,'LineWidth',2); hold on;
                subplot(1,5,5);
                plot(stopR-startR,i*spacing,'-','Color',currColor,'LineWidth',2); hold on;
            end    
        end
        
    case 'processing'
        spacing         = 0.05;
        minThresh       = min(ContData.data.trialParams(:,6));
        maxThresh       = max(ContData.data.trialParams(:,6));
        % shift           = ContData.data.trajectories.contXY(1) - ContData.data.trialTimes(1,2) + 1000;
        numTrials       = size(ContData.data.trialTimes,1);
        shift           = mean(ContData.data.trialTimes(2:numTrials,2) - ContData.data.trialTimes(1:numTrials-1,4));

        rewardTimes     = ContData.data.trialTimes(:,4) - shift;
        eventTimes      = ContData.data.trialTimes(:,3) - shift;
        startTimes      = ContData.data.trialTimes(:,2) - shift;
        offset = 9;
        
        for i=1+offset:numel(rewardTimes)
    
            startR      = find(ContData.data.trajectories.contXY(:,1) >= rewardTimes(i-1) , 1)-10;
            stopR       = find(ContData.data.trajectories.contXY(:,1) >= rewardTimes(i) , 1);
            trialThresh = ContData.data.trialParams(i,6);
            eventTime   = find(ContData.data.trajectories.contXY(:,1) >= eventTimes(i) , 1) - startR;
            startTime   = find(ContData.data.trajectories.contXY(:,1) >= startTimes(i) , 1) - startR;

            yVals       = sgolayfilt(ContData.data.trajectories.contXY(startR:stopR, 3) , 5 , 11);
            xVals       = sgolayfilt(ContData.data.trajectories.contXY(startR:stopR, 2) , 5 , 11); 

            dX = diff(xVals);
            dY = diff(yVals);
            vVals = sqrt( dX.^2 + dY.^2 );

            totalDisp                   = trapz(vVals);
            currColor                   = [(trialThresh./maxThresh) 0.67.*(1-(trialThresh./maxThresh)) 1-(trialThresh./maxThresh)];

            dispData.disp(i-offset)     = totalDisp;
            dispData.thresh(i-offset)   = trialThresh;
            dispData.rrate(i-offset)    = 1000 ./ ((stopR-startR).*33);
            dispData.vel(i-offset)      = max(vVals);

            if displayOnline==1
                figure(figNum);
                subplot(1,5,1);
                plot(yVals,xVals+(i*5),'-','Color',currColor,'LineWidth',1); hold on; axis([-150 150 0 (numel(rewardTimes)+3)*5]);
                subplot(1,5,2:3);
                plot([1:numel(vVals)].*33,vVals+(i*8),'-','Color',currColor,'LineWidth',1); hold on;  axis([0 300.*33 0 (numel(rewardTimes)+3)*8]);
                plot(eventTime.*33,i*8,'o',startTime.*33,i*8,'o','Color',currColor,'LineWidth',2);
                subplot(1,5,4);
                semilogx(totalDisp,i*spacing,'-','Color',currColor,'LineWidth',2); hold on;  axis([0 1500 0 (numel(rewardTimes)+3)*spacing]);
                subplot(1,5,5);
                plot(stopR-startR,i*spacing,'-','Color',currColor,'LineWidth',2); hold on;  axis([0 400 0 (numel(rewardTimes)+3)*spacing]);
            end    
        end
        
end


% Plot the outputted results

[rho,h] = corr(dispData.rrate',dispData.disp')
reachAnalysis.dispRewCorr.rho = rho;
reachAnalysis.dispRewCorr.prob = h;

[rho,h] = corr(dispData.thresh',dispData.vel')
reachAnalysis.thrVelCorr.rho = rho;
reachAnalysis.thrVelCorr.prob = h;
p = polyfit(dispData.thresh,dispData.vel,1);

if displayOnline==1
    figure(figNum+5); clf; subplot(211); plot(1:numel(dispData.disp),dispData.disp,'o-'); hold on; plot(1-3:numel(dispData.disp)-3,medfilt1(dispData.rrate,5).*1000,'g'); plot(1:numel(dispData.disp),medfilt1(dispData.disp,4),'r-'); plot(1:numel(dispData.thresh),dispData.thresh.*10,'k-');
    axis([-20 160 0 2000]);

    figure(figNum+5); subplot(212); plot(dispData.thresh, dispData.vel, 'ko', 10:55, polyval(p,10:55) , 'r--');
    axis([0 60 0 80]);
end

disp('Completed the per reward analysis...');

reachAnalysis.dispData = dispData;


disp(' ');
disp(' ');
