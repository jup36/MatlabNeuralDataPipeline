function [sessAnalysis] = TNC_OpenFieldReport(dataStructure,figNum,range)

%% Extract the local version of the position, accel, and gyro data from the passed structure
% restrict to 45 minutes of recording

currX       = dataStructure.data.posSmth.Xrmed(range(1):range(2));
currY       = dataStructure.data.posSmth.Yrmed(range(1):range(2));
currVel     = dataStructure.data.velocity.velRM(range(1):range(2));
currStamps  = dataStructure.data.tStamps(range(1):range(2));
currAx      = dataStructure.data.accel.x(range(1):range(2));
currAy      =  dataStructure.data.accel.y(range(1):range(2));
currAz      =  dataStructure.data.accel.z(range(1):range(2));
currGx      =  dataStructure.data.gyro.x(range(1):range(2));
currGy      =  dataStructure.data.gyro.y(range(1):range(2));
currGz      =  dataStructure.data.gyro.z(range(1):range(2));

%% Main figure for output display
currFileName = [dataStructure.data.year ' | ' dataStructure.data.month ' | ' dataStructure.data.day];
h = figure(figNum);
clf;
set(h,'Name',currFileName,'Color',[1 1 1]);
debugD = 0;

%% Find center changes    
cradius=15;
clength=0;
cindex=1;
   
centers=[];

for i=2:numel(currX)
    cjump   = pdist2([currX(i),currY(i)], [currX(cindex),currY(cindex)]);
    clength = clength + pdist2([currX(i-1),currY(i-1)],[currX(i),currY(i)]);

    if cjump>cradius
        ctheta      = atan2(currY(i)-currY(cindex),currX(i)-currX(cindex));
        ctime       = currStamps(i)-currStamps(cindex);
        centers     =[centers; i currX(i) currY(i) cjump clength ctheta ctime];
        clength     =0;
        cindex      =i;
    end
end

if debugD==1
    figure(1); clf; hold on;
    plot(currX,currY,'k');
    plot(centers(:,2),centers(:,3),'ro')
end

disp('Completed calculating center changes.');
    
%% Plot the number of center changes as a function of sample

numSamples  = numel(currX);
binWidth    = (30*60);
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

%% Look for clusters
numC=5;
count = 0;
progressions = [];
tmpPlot = zeros(size(centers,1),3);
% figure(2); subplot(211); hold on;
% figure(2); subplot(212); hold on;

for i = numC:1:size(centers,1)
    
%     subplot(211); plot(var(centers(i-(numC-1):i,6)),log(mean(centers(i-(numC-1):i,7))),'k.');
%     xlabel('Variance in center change heading');
%     ylabel('Time between center changes');
%     subplot(212); plot(mean(centers(i-(numC-1):i,5)),var(centers(i-(numC-1):i,6)),'k.');
%     ylabel('Variance in center change heading');
%     xlabel('Distance travelled');     

    tmpPlot(i,1) = mean(centers(i-(numC-1):i,5));
    tmpPlot(i,2) = var(centers(i-(numC-1):i,6));
    tmpPlot(i,3) = log(mean(centers(i-(numC-1):i,7)));
        
    if mean(centers(i-(numC-1):i,5))<=28 && var(centers(i-(numC-1):i,6))<=1 && mean(centers(i-(numC-1):i,7))<=400
        progressions = [progressions; centers(i-(numC-1),1), centers(i-(numC-1),2), centers(i-(numC-1),3), i-(numC-1), currStamps(centers(i-(numC-1),1))];
    end
    
end

figure(figNum);
subplot(6,4,2);
loglog(tmpPlot(:,1),tmpPlot(:,2),'k.','MarkerSize',2);
ylabel('Variance in center change heading');
xlabel('Distance travelled');
% zlabel('Time between center changes');
grid on; axis([10 max(tmpPlot(:,1)) 0.001 max(tmpPlot(:,2)).*10]);

% figure(1); clf; hold on;
% plot(currX,currY,'-','Color',[0, 0, 0],'LineWidth',2);
% plot(centers(:,2),centers(:,3),'go','LineWidth',2);
% plot(progressions(:,2),progressions(:,3),'bo','LineWidth',2);

%% Find the beginning and end of progressions
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

% figure(3); clf;
% plot(1:numel(contig),contig,'k',progStarts,ones(1,numel(progStarts)),'ro',progStops,ones(1,numel(progStops)),'bo');

figure(figNum);
subplot(6,4,[5,9,13,17,21]); 
hold off;
plot3(currX,currY,currStamps,'-','Color',[0.75 0.75 0.75],'LineWidth',1);
hold on;
plot3(progressions(:,2),progressions(:,3),progressions(:,5),'o','Color',[0 0.67 1],'LineWidth',1);
plot3(progStarts(:,2),progStarts(:,3),progStarts(:,1),'o','Color',[1 0 0],'LineWidth',2);
plot3(progStops(:,2),progStops(:,3),progStops(:,1),'s','Color',[0 0.33 0.5],'LineWidth',2);
grid on; view([-35 5]); axis([100 450 50 350 0 max(progressions(:,5))]); 
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Time (samples)');

%% Plot all progressions
numProgs = size(progStarts,1);
twoD = 0;

figure(figNum);
subplot(6,4,[5,9,13,17,21]+1); 
hold off;

maxLength=max(progStops(:,4)-progStarts(:,4));

for i=1:numProgs
    if twoD
        subplot(ceil(numProgs./8),8,i);
    end
    currStart = progStarts(i,4);
    currStop = progStops(i,4);
    if twoD
        plot(currX(currStart-20:currStop)-currX(currStart),currY(currStart-20:currStop)-currY(currStart),'Color',[i./numProgs,0.6.*(1-(i./numProgs)),1-(i./numProgs)],'LineWidth',2);
        hold on;
        axis([-250 250 -250 250]); grid on;
    else
%         plot3(currX(currStart-20:currStop),currY(currStart-20:currStop),currStamps(currStart-20:currStop),'Color',[i./numProgs,0.6.*(1-(i./numProgs)),1-(i./numProgs)],'LineWidth',2);
        plot3(currX(currStart-20:currStop),currY(currStart-20:currStop),currStamps(currStart-20:currStop),'Color',[(currStop-currStart)./maxLength, 0.6.*(1-((currStop-currStart)./maxLength)), 1-((currStop-currStart)./maxLength)],'LineWidth',2);
        hold on;
        grid on; view([-35 5]);
    end
end
axis([100 450 50 350 0 max(progressions(:,5))]);
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Time (samples)');

%% Histogram of progLengths
clear progLengths 
numProgs = size(progStarts,1);
for i=1:numProgs
    currStart   = progStarts(i,4);
    currStop    = progStops(i,4);

    progLengths.time(i) = currStop-currStart;
    
    incrDist = 0;
    for j=1:(currStop-currStart)
        incrDist = incrDist + pdist2([currX(currStart+(j-1)),currY(currStart+(j-1))],[currX(currStart+(j)),currY(currStart+(j))]);
    end
    progLengths.disp(i) = incrDist;

end

progLengths.timeH = hist(progLengths.time,0:3:150);
progLengths.dispH = hist(progLengths.disp,0:20:500);
progLengths.timeHx = 0:3:150;
progLengths.dispHx = 0:20:500;

figure(figNum);
subplot(6,4,4);
bar(progLengths.timeHx,progLengths.timeH,'k');
axis([0 200 0 max(progLengths.timeH)+10]);
xlabel('Time'); ylabel('Count');
subplot(6,4,8);
bar(progLengths.dispHx,progLengths.dispH,'k');
axis([0 750 0 max(progLengths.dispH)+10])
xlabel('Displacement'); ylabel('Count');

figure(figNum);
subplot(6,4,[3,7]);
plot(progLengths.time,progLengths.disp,'k.','MarkerSize',8);
axis([0 200 0 750])
xlabel('Time'); ylabel('Displacement');

%% Look at accelerations prior to progressions

numProgs = size(progStarts,1)
sampThresh = 20;
dispThresh = 80;
count = 0;

clear accelProgs gyroProgs position

for i=1:numProgs
    currStart = progStarts(i,4);
    currStop  = progStops(i,4);

%     if progLengths.time(i) > sampThresh
    if progLengths.disp(i) > dispThresh & currStart+80<numel(currAx)
        count = count+1;
        accelProgs.x(count,:) = currAx(currStart-20:currStart+80);
        accelProgs.y(count,:) = currAy(currStart-20:currStart+80);
        accelProgs.z(count,:) = currAz(currStart-20:currStart+80);

        gyroProgs.x(count,:) = abs(currGx(currStart-20:currStart+80)-mean(currGx(currStart-20:currStart+80)));
        gyroProgs.y(count,:) = abs(currGy(currStart-20:currStart+80)-mean(currGy(currStart-20:currStart+80)));
        gyroProgs.z(count,:) = abs(currGz(currStart-20:currStart+80)-mean(currGz(currStart-20:currStart+80)));
        
        position.x(count,:)  = currX(currStart-20:currStart+80) - currX(currStart);
        position.y(count,:)  = currY(currStart-20:currStart+80) - currY(currStart);
        
        for j=-20:80
            if j<0
                incrDist = pdist2([currX(currStart),currY(currStart)],[currX(currStart+(j)),currY(currStart+(j))]);
                position.d(count,j+21) = -incrDist;
            else
                incrDist = pdist2([currX(currStart),currY(currStart)],[currX(currStart+(j)),currY(currStart+(j))]);
                position.d(count,j+21) = incrDist;
            end
        end
        
        position.t = -20:80;
        currProgInstVel             = currVel(currStart:currStop);
        progLengths.curVelProf(i).v = currVel(currStart:currStop);        
        progLengths.maxVel(i)       = max(currProgInstVel);
        progLengths.maxVelInd(i)    = find(currProgInstVel==max(currProgInstVel),1);
        progLengths.maxVelIndF(i)   = find(currProgInstVel==max(currProgInstVel),1) ./ (currStop-currStart);
%         figure(5)
%         plot(progLengths.curVelProf(i).v,'k');
%         drawnow;
    else
        progLengths.curVelProf(i).v = 0;        
        progLengths.maxVel(i)       = 0;
        progLengths.maxVelInd(i)    = 0;
        progLengths.maxVelIndF(i)   = 0;
    end
end

figure(figNum);
subplot(6,4,[11,15]); hold off;
plot(progLengths.maxVel,progLengths.maxVelIndF,'k.','MarkerSize',8);
ylabel('Index of Max Velocity');
xlabel('Max Velocity');
axis([0 30 0 1]);
subplot(6,4,[19,23]); hold off;
plot(progLengths.maxVel, progLengths.disp, 'k.','MarkerSize',8);
ylabel('Total Displacement');
xlabel('Max Velocity');
axis([0 30 1 750]);
% patch(-80:60,mean(accelProgs.x,1) + (std(accelProgs.x,[],1)./sqrt(numProgs)),'k')

figure(figNum); 
subplot(6,4,[11,15]+1);
% hold off;
% plot(mean(position.d,1),mean(accelProgs.x,1)- mean(mean(accelProgs.x(:,1:20),1)),'LineWidth',2,'Color',[0.5 0.5 0.5])
% hold on;
% plot(mean(position.d,1),mean(accelProgs.y,1)- mean(mean(accelProgs.y(:,1:20),1)),'r','LineWidth',2)
% plot(mean(position.d,1),mean(accelProgs.z,1)- mean(mean(accelProgs.z(:,1:20),1)),'LineWidth',2,'Color',[0 0.67 1])
% plot(mean(position.d,1),mean(accelProgs.x,1)- mean(mean(accelProgs.x(:,1:20),1)) + (std(accelProgs.x,[],1)./sqrt(numProgs)),'k',mean(position.d,1),mean(accelProgs.x,1)- mean(mean(accelProgs.x(:,1:20),1)) - (std(accelProgs.x,[],1)./sqrt(numProgs)),'k','LineWidth',1,'Color',[0.5 0.5 0.5])
% plot(mean(position.d,1),mean(accelProgs.y,1)- mean(mean(accelProgs.y(:,1:20),1)) + (std(accelProgs.y,[],1)./sqrt(numProgs)),'r',mean(position.d,1),mean(accelProgs.y,1)- mean(mean(accelProgs.y(:,1:20),1)) - (std(accelProgs.y,[],1)./sqrt(numProgs)),'r','LineWidth',1)
% plot(mean(position.d,1),mean(accelProgs.z,1)- mean(mean(accelProgs.z(:,1:20),1)) + (std(accelProgs.z,[],1)./sqrt(numProgs)),'b',mean(position.d,1),mean(accelProgs.z,1)- mean(mean(accelProgs.z(:,1:20),1)) - (std(accelProgs.z,[],1)./sqrt(numProgs)),'b','LineWidth',1,'Color',[0 0.67 1])
% plot(mean(position.d,1),zeros(1,size(position.d,2)),'k--');
plot([0 0],[-150 150],'k--');
xlabel('Displacement from start of progression (pixels)');

hold off;
plot(mean(position.t,1),mean(accelProgs.x,1)- mean(mean(accelProgs.x(:,1:20),1)),'LineWidth',2,'Color',[0.5 0.5 0.5])
hold on;
plot(mean(position.t,1),mean(accelProgs.y,1)- mean(mean(accelProgs.y(:,1:20),1)),'r','LineWidth',2)
plot(mean(position.t,1),mean(accelProgs.z,1)- mean(mean(accelProgs.z(:,1:20),1)),'LineWidth',2,'Color',[0 0.67 1])
plot(mean(position.t,1),mean(accelProgs.x,1)- mean(mean(accelProgs.x(:,1:20),1)) + (std(accelProgs.x,[],1)./sqrt(numProgs)),'k',mean(position.t,1),mean(accelProgs.x,1)- mean(mean(accelProgs.x(:,1:20),1)) - (std(accelProgs.x,[],1)./sqrt(numProgs)),'k','LineWidth',1,'Color',[0.5 0.5 0.5])
plot(mean(position.t,1),mean(accelProgs.y,1)- mean(mean(accelProgs.y(:,1:20),1)) + (std(accelProgs.y,[],1)./sqrt(numProgs)),'r',mean(position.t,1),mean(accelProgs.y,1)- mean(mean(accelProgs.y(:,1:20),1)) - (std(accelProgs.y,[],1)./sqrt(numProgs)),'r','LineWidth',1)
plot(mean(position.t,1),mean(accelProgs.z,1)- mean(mean(accelProgs.z(:,1:20),1)) + (std(accelProgs.z,[],1)./sqrt(numProgs)),'b',mean(position.t,1),mean(accelProgs.z,1)- mean(mean(accelProgs.z(:,1:20),1)) - (std(accelProgs.z,[],1)./sqrt(numProgs)),'b','LineWidth',1,'Color',[0 0.67 1])
plot(mean(position.t,1),zeros(1,size(position.t,2)),'k--');

plot([0 0],[-150 150],'k--');
xlabel('Time from start of progression (pixels)');

ylabel('Mean Accelerometer Value (a.u.)');

figure(figNum);
subplot(6,4,[19,23]+1);
% hold off;
% plot(mean(position.d,1),mean(gyroProgs.x,1),'k','LineWidth',2,'Color',[0.5 0.5 0.5])
% hold on;
% plot(mean(position.d,1),1000 + mean(gyroProgs.y,1),'r','LineWidth',2)
% plot(mean(position.d,1),2000 + mean(gyroProgs.z,1),'b','LineWidth',2,'Color',[0 0.67 1])
% plot(mean(position.d,1),mean(gyroProgs.x,1) + (std(gyroProgs.x,[],1)./sqrt(numProgs)),'k',mean(position.d,1),mean(gyroProgs.x,1) - (std(gyroProgs.x,[],1)./sqrt(numProgs)),'k','LineWidth',1,'Color',[0.5 0.5 0.5])
% plot(mean(position.d,1),1000 + mean(gyroProgs.y,1) + (std(gyroProgs.y,[],1)./sqrt(numProgs)),'r',mean(position.d,1),1000 + mean(gyroProgs.y,1) - (std(gyroProgs.y,[],1)./sqrt(numProgs)),'r','LineWidth',1)
% plot(mean(position.d,1),2000 + mean(gyroProgs.z,1) + (std(gyroProgs.z,[],1)./sqrt(numProgs)),'b',mean(position.d,1),2000 + mean(gyroProgs.z,1) - (std(gyroProgs.z,[],1)./sqrt(numProgs)),'b','LineWidth',1,'Color',[0 0.67 1])
% xlabel('Displacement from start of progression (pixels)');

hold off;
plot(mean(position.t,1),mean(gyroProgs.x,1),'k','LineWidth',2,'Color',[0.5 0.5 0.5])
hold on;
plot(mean(position.t,1),1000 + mean(gyroProgs.y,1),'r','LineWidth',2)
plot(mean(position.t,1),2000 + mean(gyroProgs.z,1),'b','LineWidth',2,'Color',[0 0.67 1])
plot(mean(position.t,1),mean(gyroProgs.x,1) + (std(gyroProgs.x,[],1)./sqrt(numProgs)),'k',mean(position.t,1),mean(gyroProgs.x,1) - (std(gyroProgs.x,[],1)./sqrt(numProgs)),'k','LineWidth',1,'Color',[0.5 0.5 0.5])
plot(mean(position.t,1),1000 + mean(gyroProgs.y,1) + (std(gyroProgs.y,[],1)./sqrt(numProgs)),'r',mean(position.t,1),1000 + mean(gyroProgs.y,1) - (std(gyroProgs.y,[],1)./sqrt(numProgs)),'r','LineWidth',1)
plot(mean(position.t,1),2000 + mean(gyroProgs.z,1) + (std(gyroProgs.z,[],1)./sqrt(numProgs)),'b',mean(position.t,1),2000 + mean(gyroProgs.z,1) - (std(gyroProgs.z,[],1)./sqrt(numProgs)),'b','LineWidth',1,'Color',[0 0.67 1])
xlabel('Time from start of progression (pixels)');

plot([0 0],[0 3500],'k--');
ylabel('Mean Gyro Value (a.u.)');

%% Pack output
sessAnalysis.centers = centers;
sessAnalysis.progressions = progressions;
sessAnalysis.position = position;
sessAnalysis.progLengths = progLengths;
sessAnalysis.accelProgs = accelProgs;
sessAnalysis.gyroProgs = gyroProgs;
sessAnalysis.progStarts = progStarts;
sessAnalysis.progStops = progStops;
sessAnalysis.range = range;
