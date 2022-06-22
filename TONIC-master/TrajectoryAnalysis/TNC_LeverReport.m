function [sessAnalysis] = TNC_LeverReport(dataStructure,figNum)

%% Create a local version of the position, lick, and trial times
% restrict to "range" of recording
range(1)    = find(dataStructure.data.cont.ts>dataStructure.data.trialData.eTs(10),1);
maxTrials   = numel(dataStructure.data.trialData.eTs);
range(2)    = find(dataStructure.data.cont.ts<dataStructure.data.trialData.eTs(maxTrials),1,'last');

currX       = dataStructure.data.cont.xSR(range(1):range(2));
currY       = dataStructure.data.cont.ySR(range(1):range(2));
currV       = dataStructure.data.cont.velRM(range(1):range(2));
currA       = dataStructure.data.cont.angVS(range(1):range(2));
currStamps  = dataStructure.data.cont.ts(range(1):range(2));
currLick    = dataStructure.data.cont.lick(range(1):range(2));

currTrialTs = dataStructure.data.trialData.eTs;
currRewTs   = dataStructure.data.trialData.rTs;
currTrialTh = dataStructure.data.trialData.thresh;

%% Get block structure timestamps
blockData.inds  = [15 ; find(abs(diff(currTrialTh))>0)];
blockData.ts    = currTrialTs(blockData.inds);

for p = 1:numel(blockData.inds)
    blockData.cInds(p) = currStamps(find(currStamps>blockData.ts(p),1,'first'));
end
 blockData.cInds
%% Main figure for output display
currFileName = [dataStructure.data.mouse ' ... ' dataStructure.data.year '.' dataStructure.data.month '.' dataStructure.data.day];
h = figure(figNum);
clf;
set(h,'Name',currFileName,'Color',[1 1 1]);
debugD = 0;

%% Find movements

% method 1 is based upon 'center changes'
% method 2 is based upon thresholding the velocities

method = 'vel';
numC = 5;

switch method
    
    case 'cc'
        cradius=5;
        clength=0;
        cindex=1;

        centers=[];

        for i=2:numel(currX)
            cjump   = pdist2([currX(i),currY(i)], [currX(cindex),currY(cindex)]);
            caway   = pdist2([currX(i),currY(i)], [0,0]) - pdist2([currX(cindex),currY(cindex)],[0,0]);
            clength = clength + pdist2([currX(i-1),currY(i-1)],[currX(i),currY(i)]);

            if cjump>cradius
                ctheta      = atan2(currY(i)-currY(cindex),currX(i)-currX(cindex));
                ctime       = currStamps(i)-currStamps(cindex);
                centers     =[centers; i currX(i) currY(i) cjump clength ctheta ctime caway];
                clength     =0;
                cindex      =i;
            end
        end

        if debugD==1
            figure(21); clf; hold on;
            plot(1:numel(currX),currX,'b');
            plot(1:numel(currY),currY,'k');
            plot(centers(:,1),centers(:,2),'ro');
            plot(centers(:,1),centers(:,3),'ro')
        end

        disp('Completed calculating center changes.');

        % Plot the number of center changes as a function of sample
        numSamples  = numel(currX);
        binWidth    = (30*60*2);
        binCnt = 0;
        clear numChanges %= zeros(1,numSamples./binWidth);

        for j=1:binWidth:numSamples
            binCnt = binCnt+1;
            numChanges(binCnt,1) = numel(find(centers(:,1)>=j & centers(:,1)<j+binWidth-1));
        end

        figure(figNum);
        subplot(6,4,1);
        plot(numChanges,'k.-','MarkerSize',8); axis([-10 numel(numChanges)+10 0 max(numChanges)+20]);
        grid on; xlabel(['Bins (' num2str(ceil(binWidth.*0.033)) ' sec/bin)']); ylabel('Center changes');

        % Look for clusters
        count = 0;
        progressions = [];
        tmpPlot = zeros(size(centers,1),3);
        % figure(2); subplot(211); hold on;
        % figure(2); subplot(212); hold on;

        for i = numC:1:size(centers,1)

            tmpPlot(i,1) = mean(centers(i-(numC-1):i,5));
            tmpPlot(i,2) = var(centers(i-(numC-1):i,6));
            tmpPlot(i,3) = log(mean(centers(i-(numC-1):i,7)));
            tmpPlot(i,4) = mean(centers(i-(numC-1):i,8));

            if mean(centers(i-(numC-1):i,8))>0 %&& log(mean(centers(i-(numC-1):i,7)))<=6.5
                progressions = [progressions; centers(i-(numC-1),1), centers(i-(numC-1),2), centers(i-(numC-1),3), i-(numC-1), currStamps(centers(i-(numC-1),1))];
            end

        end

        figure(figNum);
        subplot(6,4,2);
        plot(tmpPlot(:,4),tmpPlot(:,3),'k.','MarkerSize',2);
        ylabel('Movement duration');
        xlabel('Centripetal / Centrifugal');
        % zlabel('Time between center changes');
        grid on; 
        % axis([-1 max(tmpPlot(:,1)) 0.001 max(tmpPlot(:,2)).*10]);

        % figure(1); clf; hold on;
        % plot(currX,currY,'-','Color',[0, 0, 0],'LineWidth',2);
        % plot(centers(:,2),centers(:,3),'go','LineWidth',2);
        % plot(progressions(:,2),progressions(:,3),'bo','LineWidth',2);

        % Find the beginning and end of progressions
        numProg = 0;
        count = 0;
        clear progS*
        contig = diff(progressions(:,4));

        for i=2:numel(contig)-1
            if contig(i)==1 && contig(i-1)>1 && contig(i+1)==1
                count = count+1;
                progStarts(count,1) = currStamps(progressions(i,1));
                progStarts(count,2) = currX(progressions(i,1));
                progStarts(count,3) = currY(progressions(i,1));
                progStarts(count,4) = progressions(i,1);
            else
                if count>0
                    if contig(i-1)==1 && contig(i)==1 && contig(i+1)>1
                        progStops(count,1) = currStamps(progressions(i+1,1));
                        progStops(count,2) = currX(progressions(i+1,1));
                        progStops(count,3) = currY(progressions(i+1,1));
                        progStops(count,4) = progressions(i+1,1);
                    end
                end
            end
        end

        % must add in the final point
        i = numel(contig);
        progStops(count,1) = currStamps(progressions(i+1,1));
        progStops(count,2) = currX(progressions(i+1,1));
        progStops(count,3) = currY(progressions(i+1,1));
        progStops(count,4) = progressions(i+1,1);

    case 'vel'
        
        % threshold the velocities
        allValidSamps   = find(currV>0.5);
        
        pre = 5;
        post = 3;
        
        count = 1;
        progStartTMP(count,1)  = currStamps(allValidSamps(1));
        progStartTMP(count,2)  = currX(allValidSamps(1));
        progStartTMP(count,3)  = currY(allValidSamps(1));
        progStartTMP(count,4)  = allValidSamps(1);
        
        for j=2:numel(allValidSamps)
            if allValidSamps(j)>allValidSamps(j-1)+5
                progStopTMP(count,1)   = currStamps(allValidSamps(j-1)+post);
                progStopTMP(count,2)   = currX(allValidSamps(j-1)+post);
                progStopTMP(count,3)   = currY(allValidSamps(j-1)+post);
                progStopTMP(count,4)   = allValidSamps(j-1)+post;

                count                  = count+1;

                progStartTMP(count,1)  = currStamps(allValidSamps(j)-pre);
                progStartTMP(count,2)  = currX(allValidSamps(j)-pre);
                progStartTMP(count,3)  = currY(allValidSamps(j)-pre);
                progStartTMP(count,4)  = allValidSamps(j)-pre;                
            end
            
            if j==numel(allValidSamps)
                progStopTMP(count,1)   = currStamps(allValidSamps(j)+post);
                progStopTMP(count,2)   = currX(allValidSamps(j)+post);
                progStopTMP(count,3)   = currY(allValidSamps(j)+post);
                progStopTMP(count,4)   = allValidSamps(j)+post;
            end
            
        end
        
        count = 1;
        for k = 1:size(progStartTMP,1)
            
            % 'progressions' must be centripetal and at least multiple samples 
            if progStopTMP(k,4)-progStartTMP(k,4)>=4             

                trajAngle   = atan2(progStopTMP(k,3)-progStartTMP(k,3),progStopTMP(k,2)-progStartTMP(k,2)) + pi;
                returnAngle = atan2(-progStartTMP(k,3),-progStartTMP(k,2)) + pi;
                
                if abs(trajAngle-returnAngle)>0.22 | (pdist2([progStopTMP(k,2),progStopTMP(k,3)],[0,0]) > pdist2([progStartTMP(k,2),progStartTMP(k,3)],[0,0]))% +/-22.5 degrees
                    progStartTMP2(count,:) = progStartTMP(k,:);
                    progStopTMP2(count,:)  = progStopTMP(k,:);
                    count               = count+1;
                end
            end            
        end
        
        count = 1;
        numCentrip = size(progStartTMP2,1)
        k = 1;
        
        while k < numCentrip
            
            if pdist2([progStartTMP2(k,2),progStartTMP2(k,3)],[0,0])<10

                % valid start found
                progStarts(count,:) = progStartTMP2(k,:);

                % now search for the end
                m = k+1;
                
                while m < numCentrip

                    if pdist2([progStartTMP2(m,2),progStartTMP2(m,3)],[0,0])<10
                        % m-1 is the end of a complete trajectory
                        origStart   = k;
                        mergeEnd    = m-1;
                        disp(['Original start index: ' num2str(origStart)]);
                        disp(['Merged final index: ' num2str(mergeEnd)]);
                        break % the while loop on m
                    else
                        m = m + 1;
                    end
                    
                end

                progStarts(count,:) = progStartTMP2(origStart,:);
                progStops(count,:)  = progStopTMP2(mergeEnd,:);
                count               = count+1;
                k                   = m;
            end
            
        end
        
        progressions(:,1) = currStamps(allValidSamps);
        progressions(:,2) = currX(allValidSamps);
        progressions(:,3) = currY(allValidSamps);
        
        % Plot the number of movements as a function of sample
        numSamples  = numel(currX);
        binWidth    = (30*60*3);
        binCnt = 0;
        clear numChanges %= zeros(1,numSamples./binWidth);

        for j=1:numel(blockData.cInds)
            if j==numel(blockData.cInds)
                numChanges(j,1) = numel(find(progStartTMP2(:,1)>=blockData.cInds(j) & progStartTMP2(:,1)<currStamps(numel(currStamps))));                
                numChanges(j,1) = 1000 .* numChanges(j,1) ./ (currStamps(numel(currStamps))-blockData.cInds(j))
            else
                numChanges(j,1) = numel(find(progStartTMP2(:,1)>=blockData.cInds(j) & progStartTMP2(:,1)<blockData.cInds(j+1)));
                numChanges(j,1) = 1000 .* numChanges(j,1) ./ (blockData.cInds(j+1)-blockData.cInds(j))
            end
        end

        figure(figNum);
        subplot(6,4,1);
        numChanges
        plot(numChanges,'k.-','MarkerSize',8); axis([-10 numel(numChanges)+10 0 max(numChanges)+(max(numChanges).*0.1)]);
        grid on; xlabel(['Blocks']); ylabel('Mvmt Rate (Hz)');

        
end

figure(figNum);
subplot(6,4,[5,9,13,17,21]); 
hold off;
plot3(currX,currY,currStamps,'-','Color',[0.75 0.75 0.75],'LineWidth',1);
hold on;
plot3(progressions(:,2),progressions(:,3),progressions(:,1),'Color',[0 0.67 1],'LineWidth',1);
plot3(progStarts(:,2),progStarts(:,3),progStarts(:,1),'o','Color',[1 0 0],'LineWidth',1);
plot3(progStops(:,2),progStops(:,3),progStops(:,1),'s','Color',[0 0.33 0.5],'LineWidth',1);

for p = 1:numel(blockData.inds)
%     pf = patch([-80 -80 80 80],[-80 80 80 -80],blockData.cInds(p).*ones(4,1),'k');
%     set(pf,'FaceColor',[0.1 0.1 0.1],'FaceAlpha',0.1);
    plot3([80 80],[80 -80],blockData.cInds(p).*ones(2,1),'k','LineWidth',2);
end

grid on; view([-35 5]); axis([-80 80 -80 80 dataStructure.data.cont.ts(range(1)) dataStructure.data.cont.ts(range(2))]); 
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Time (samples)');

%% Plot all progressions
numProgs = size(progStarts,1);
twoD = 0;

prev = 0;

figure(figNum);
subplot(6,4,[5,9,13,17,21]+1); 
hold off;

maxLength=log(max(progStops(:,4)-progStarts(:,4)));

for i=1:numProgs
    if twoD
        subplot(ceil(numProgs./8),8,i);
    end
    currStart = progStarts(i,4);
    currStop = progStops(i,4);
    if twoD
        plot(currX(currStart-prev:currStop)-currX(currStart),currY(currStart-prev:currStop)-currY(currStart),'Color',[i./numProgs,0.6.*(1-(i./numProgs)),1-(i./numProgs)],'LineWidth',2);
        hold on;
        %axis([-80 80 -80 80]); grid on;
    else
        plot3(currX(currStart-prev:currStop),currY(currStart-prev:currStop),currStamps(currStart-prev:currStop),'Color',[1-(log(currStop-currStart)./maxLength), 0.6.*(log(currStop-currStart)./maxLength), log(currStop-currStart)./maxLength],'LineWidth',2);
        hold on;
        grid on; view([-35 5]);
    end
    
    
end
figure(figNum);
plot3([0,0],[0,0],[currStamps(1),currStamps(numel(currStamps))],'k','LineWidth',1);
for p = 1:numel(blockData.inds)
%     pf = patch([-80 -80 80 80],[-80 80 80 -80],blockData.cInds(p).*ones(4,1),'k');
%     set(pf,'FaceColor',[0.1 0.1 0.1],'FaceAlpha',0.1);
    plot3([0 0],[80 -80],blockData.cInds(p).*ones(2,1),'k--','LineWidth',1);
end
axis([-80 80 -80 80 dataStructure.data.cont.ts(range(1)) dataStructure.data.cont.ts(range(2))]);
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Time (samples)');

% figure(4); clf;
% plot(currStamps,currV,'k.-'); hold on;
% plot(currStamps(allValidSamps),currV(allValidSamps),'c.');
% for i=1:numProgs
%     currStart = progStarts(i,4);
%     currStop = progStops(i,4);
%     plot(currStamps(currStart-prev:currStop),currV(currStart-prev:currStop),'ro');hold on;
% end

%% Histogram of progLengths

clear progLengths
numProgs = size(progStarts,1)
for i=1:numProgs
    currStart   = progStarts(i,4);
    currStop    = progStops(i,4);

    progLengths.time(i) = (progStops(i,1)-progStarts(i,1))./1000;
    
    incrDist = 0;
    clear vel;
    for j=1:(currStop-currStart)
        vel(j)   = pdist2([currX(currStart+(j-1)),currY(currStart+(j-1))],[currX(currStart+(j)),currY(currStart+(j))]);
        incrDist = incrDist + vel(j);
    end
    
    progLengths.disp(i) = incrDist;
    progLengths.tort(i) = incrDist ./ pdist2([currX(currStart),currY(currStart)] , [currX(currStop),currY(currStop)]);
    progLengths.vmax(i) = max(currV(currStart:currStop));
    progLengths.jerk(i) = sum(diff(diff(currV(currStart:currStop))).^2) ./ progLengths.disp(i);
    progLengths.amax(i) = max(currA(currStart:currStop));
    progLengths.asum(i) = trapz(currA(currStart:currStop));
    progLengths.pfra(i) = numel(find(currV(currStart:currStop)<0.2))./(currStop-currStart);
    
    
%     figure(24);
%     subplot(311);
%     plot(vel);
%     subplot(312);
%     plot(diff(vel));
%     subplot(313);
%     plot(diff(diff(vel)));
%     drawnow; pause(0.1);

end

progLengths.timeHx = 10.^(-1:0.1:1);
progLengths.dispHx = 10.^(0:0.1:3);
progLengths.tortHx = 10.^(-0.1:0.1:2);
progLengths.vmaxHx = 10.^(-1:0.1:2);
progLengths.jerkHx = 10.^(-3:0.1:2);
progLengths.amaxHx = 10.^(-1:0.1:2);
progLengths.asumHx = 10.^(-1:0.1:3);
progLengths.pfraHx = 0:0.05:1;

progLengths.timeH = hist(progLengths.time,progLengths.timeHx);
progLengths.dispH = hist(progLengths.disp,progLengths.dispHx);
progLengths.tortH = hist(progLengths.tort,progLengths.tortHx);
progLengths.vmaxH = hist(progLengths.vmax,progLengths.vmaxHx);
progLengths.jerkH = hist(progLengths.jerk,progLengths.jerkHx);
progLengths.amaxH = hist(progLengths.amax,progLengths.amaxHx);
progLengths.asumH = hist(progLengths.asum,progLengths.asumHx);
progLengths.pfraH = hist(progLengths.pfra,progLengths.pfraHx);

figure(figNum);
subplot(6,4,3);
semilogx(progLengths.disp,progLengths.pfra,'k.');
axis([0 max(progLengths.disp) 0 1]);
xlabel('Displacement'); ylabel('Fraction paused');

figure(figNum);
subplot(6,4,4);
semilogx(progLengths.timeHx,progLengths.timeH,'k');
axis([min(progLengths.timeHx) max(progLengths.timeHx) 0 max(progLengths.timeH)]);
xlabel('Time'); ylabel('Count');

figure(figNum);
subplot(6,4,7);
semilogx(progLengths.tortHx,progLengths.tortH,'k');
axis([min(progLengths.tortHx) max(progLengths.tortHx) 0 max(progLengths.tortH)])
xlabel('Tortuosity'); ylabel('Count');

figure(figNum);
subplot(6,4,8);
semilogx(progLengths.dispHx,progLengths.dispH./max(progLengths.dispH),'k');
axis([min(progLengths.dispHx) max(progLengths.dispHx) 0 1])
xlabel('Displacement'); ylabel('Count');

figure(figNum);
subplot(6,4,11);
semilogx(progLengths.vmaxHx,progLengths.vmaxH,'k');
axis([min(progLengths.vmaxHx) max(progLengths.vmaxHx) 0 max(progLengths.vmaxH)])
xlabel('Max Velocity'); ylabel('Count');

figure(figNum);
subplot(6,4,12);
semilogx(progLengths.jerkHx,progLengths.jerkH,'k');
axis([min(progLengths.jerkHx) max(progLengths.jerkHx) 0 max(progLengths.jerkH)])
xlabel('Total Jerk / displacement'); ylabel('Count');

figure(figNum);
subplot(6,4,15);
semilogx(progLengths.amaxHx,progLengths.amaxH,'k');
axis([min(progLengths.amaxHx) max(progLengths.amaxHx) 0 max(progLengths.amaxH)])
xlabel('Max Angular Velocity'); ylabel('Count');

figure(figNum);
subplot(6,4,16);
semilogx(progLengths.asumHx,progLengths.asumH,'k');
axis([min(progLengths.asumHx) max(progLengths.asumHx) 0 max(progLengths.asumH)])
xlabel('Total Angular Velocity'); ylabel('Count');

%% Find all progressions that move away from the center


%% Look at lever movements around trials
clear trialData
numTrials   = numel(currTrialTs);
window      = -7:66;

% figure(22); clf;
% figure(23); clf;

for j=16:numTrials-1
    
    eventIndex                  = find(currStamps > currTrialTs(j), 1);
    progIndex                   = find(progStarts(:,1) < currTrialTs(j), 1, 'last');

    trialData.trajX(j-15,:)     = currX(window+eventIndex);
    trialData.trajY(j-15,:)     = currY(window+eventIndex);
    trialData.vel(j-15,:)       = currV(window+eventIndex);
    trialData.lick(j-15,:)      = currLick(window+eventIndex);
    
%     figure(22); subplot(6,15,j-15); hold on;
%     plot([0 0],[-50 50],'k--',[-50 50],[0 0],'k--');
%     plot(trialData.trajX(j-15,1:7),trialData.trajY(j-15,1:7),'r.-');
%     plot(trialData.trajX(j-15,7:40),trialData.trajY(j-15,7:40),'k.-');
%     plot(trialData.trajX(j-15,40:73),trialData.trajY(j-15,40:73),'c.-');
%     axis([-50 50 -50 50]);
% 
%     figure(23); subplot(6,15,j-15); hold on;
%     plot([0 0],[-50 50],'k--',[1000 1000],[-50 50],'k--');
%     plot(window.*33,trialData.vel(j-15,:).*2,'r',window.*33,trialData.lick(j-15,:)./10,'b');
%     axis([-7*33 66*33 -1 50]);

end

%% Pack output
% sessAnalysis.centers = centers;
sessAnalysis.progressions = progressions;
% sessAnalysis.position = position;
sessAnalysis.progLengths = progLengths;
% sessAnalysis.accelProgs = accelProgs;
% sessAnalysis.gyroProgs = gyroProgs;
sessAnalysis.progStarts = progStarts;
sessAnalysis.progStops = progStops;
sessAnalysis.trialData = trialData;
sessAnalysis.blockData = blockData;

% sessAnalysis.range = range;
