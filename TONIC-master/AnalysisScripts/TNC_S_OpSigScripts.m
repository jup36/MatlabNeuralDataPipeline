%% TNC_S_OpSigScripts

%% GENERAL >>>  INITIALIZE PARAMETERS
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 100;
currParams.filter.causal       = 1;

[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

if currParams.filter.causal
    currParams.filter.kernel(1:currParams.smthParams.decay.*15) = 0;
end
currParams.offset = 0;

currParams.winParams.prior     = 2e3;
currParams.winParams.after     = 5e3; 


PopData.currParams = currParams;

disp(' ');
disp(' ');
disp('________________________________');
disp(' ');
disp('Initialized analysis parameters.');
disp('________________________________');
disp(' ');
disp(' ');

%% ALL SESSIONS >>> CONVERT TO PREFERRED DATA ORGANIZATION - DEPRECATED

totCount = 0;
numMice = size(PData.mouse,2);

for a=1:numMice

    numSess = size(PData.mouse(a).session,2);

    for b=1:numSess

        disp(['MOUSE: ' num2str(a) ' | Current session: ' num2str(b) ' of ' num2str(numSess)]);
        count       = 1;
        tmp         = find(PData.mouse(a).session(b).Nev.Data.Spikes.Unit>0);
        allTrodes   = unique(PData.mouse(a).session(b).Nev.Data.Spikes.Electrode)

        for i=1:numel(allTrodes)

            currTrode = allTrodes(i)
            
            for j=1:6

                spkInds = find(PData.mouse(a).session(b).Nev.Data.Spikes.Electrode==currTrode & PData.mouse(a).session(b).Nev.Data.Spikes.Unit==j);

                if numel(spkInds)>0

                    % index through channels looking for sorted units
                    PopData.mouse(a).session(b).unit(count).ts   = round(double(PData.mouse(a).session(b).Nev.Data.Spikes.Timestamps(spkInds))./30);
                    PopData.mouse(a).session(b).unit(count).el   = currTrode;
                    PopData.mouse(a).session(b).unit(count).un   = j;

                    % Get the electrode map
                    %[row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(currTrode,'NN_w64');
                    %PopData.session(currSes).unit(count).row = row;
                    %PopData.session(currSes).unit(count).col = col;

                    count = count+1;

                end
                
            end
            
        end

        count       = count - 1;
        totCount    = totCount + count;
        
        PopData.mouse(a).session(b).numUnits = count;
        clear count;        

    end
    
end

disp(['Total number of units: ' num2str(totCount)]);
%% ALL SESSIONS >>> CONVERT TO PREFERRED DATA ORGANIZATION - ACTIVE

totCount = [];
numMice = size(PData.mouse,2);

for a=1:numMice

    numSess = size(PData.mouse(a).session,2);

    for b=1:numSess

        disp(['MOUSE: ' num2str(a) ' | Current session: ' num2str(b) ' of ' num2str(numSess)]);

        PData.mouse(a).session(b).numUnits = numel(PData.mouse(a).session(b).unit);
        totCount = [totCount PData.mouse(a).session(b).numUnits]

        PopData.mouse(a).session(b).Behavior = PData.mouse(a).session(b).Behavior;
        
        for c=1:numUnits

            % index through channels looking for sorted units
            PopData.mouse(a).session(b).unit(c).ts   = PData.mouse(a).session(b).unit(c).ts;
            PopData.mouse(a).session(b).unit(c).el   = PData.mouse(a).session(b).unit(c).el;
            PopData.mouse(a).session(b).unit(c).un   = PData.mouse(a).session(b).unit(c).un;
            
            [PopData.mouse(a).session(b).unit(c).isi] = TNC_QuantISI(PData.mouse(a).session(b).unit(c).ts);
            
        end
            
    end
    
end

disp(['Total number of units: ' num2str(sum(totCount)) ' | ' num2str(median(totCount)) ' per session']);

% Tidy up...
% clear PData;
% PData = PopData;
% clear PopData;

%% ALL SESSIONS >>> APPROACH SEP BY SIDE

a       = 1;
currSes = 1;
ySpc    = 4;
count   = 1;

% ALIGN ON PRESS (right VS RIGHT)
for a=1:2

    for currSes = 1:size(PData.mouse(a).session,2)

        disp(['Extracting data from session ' num2str(currSes) ' of ' num2str(size(PData.mouse(a).session,2)) ' total.']);

        PData.mouse(a).session(currSes).numUnits = numel(PData.mouse(a).session(currSes).unit);
        
        for j = 1:PData.mouse(a).session(currSes).numUnits

            numStamps   = numel(PData.mouse(a).session(currSes).unit(j).ts);
            delta       = zeros(1,ceil(PData.mouse(a).session(currSes).unit(j).ts(numStamps)));
            delta(PData.mouse(a).session(currSes).unit(j).ts) = 1;
            tmpSmooth   = conv(delta,currParams.filter.kernel,'same');
            tmpSmoothZ  = ( tmpSmooth - mean(tmpSmooth) ) ./ std(tmpSmooth);

            leftInds = find(PData.mouse(a).session(currSes).Behavior.LeverSide==1 & PData.mouse(a).session(currSes).Behavior.TrialType==2);            
            [rsLevPress] = TNC_AlignRasters(tmpSmoothZ , PData.mouse(a).session(currSes).unit(j).ts , -1 ,round( PData.mouse(a).session(currSes).Behavior.Approach(leftInds)-currParams.offset ), [currParams.winParams.prior currParams.winParams.after],0,1);
            analysis.approach.left.pAvg(count,:)    = rsLevPress.image.psthAVG;
            analysis.approach.left.pZav(count,:)    = rsLevPress.image.psthZ;
            analysis.approach.left.sess(count)      = currSes;
            analysis.approach.left.unit(count)      = j;
            analysis.approach.left.mouse(count)     = a;

            rightInds = find(PData.mouse(a).session(currSes).Behavior.LeverSide==0 & PData.mouse(a).session(currSes).Behavior.TrialType==2);
            [rsLevPress] = TNC_AlignRasters(tmpSmoothZ , PData.mouse(a).session(currSes).unit(j).ts , -1 ,round( PData.mouse(a).session(currSes).Behavior.Approach(rightInds)-currParams.offset ), [currParams.winParams.prior currParams.winParams.after],0,1);
            analysis.approach.right.pAvg(count,:)   = rsLevPress.image.psthAVG;
            analysis.approach.right.pZav(count,:)   = rsLevPress.image.psthZ;
            analysis.approach.right.sess(count)     = currSes;
            analysis.approach.right.unit(count)     = j;
            analysis.approach.right.mouse(count)    = a;

%             numInds = numel(leftInds);
%             candidateInds = PData.mouse(a).session(currSes).Behavior.Approach(leftInds);
%             [sink] = TNC_ExtTrigWins(tmpSmoothZ,candidateInds,[currParams.winParams.prior currParams.winParams.after]);
%             analysis.approach.left.pZavg(count,:) = sink.avg;
%             analysis.approach.left.pZsem(count,:) = sink.err;
%             analysis.approach.left.sess(count) = currSes;
%             analysis.approach.left.unit(count) = j;
%             analysis.approach.left.mouse(count) = a;            
            
            PData.mouse(a).session(currSes).unit(j).meanRate = rsLevPress.meanRate;
            PData.mouse(a).session(currSes).unit(j).stdRate = rsLevPress.stdRate;

            count = count+1;

        end
    end
end

clear currData*

for j = 1:size(analysis.approach.right.pZav,1)
    currDataR(j) = trapz( analysis.approach.right.pZav(j,2000:6000) );
    currDataL(j) = trapz( analysis.approach.left.pZav(j,2000:6000) );
end

[vals,indsR] = sort(currDataR,'descend');
[vals,indsL] = sort(currDataL,'descend');

figure(10); clf;
[mapName] = TNC_CreateRBColormap(1024,'mbr');

subplot(4,3,[1 4 7])
imagesc(analysis.approach.left.pAvg(indsL,:),[-2 2]); hold on;
plot([2000 2000],[0 200],'k-')
colormap(mapName);

subplot(4,3,[1 4 7]+1)
imagesc(analysis.approach.right.pAvg(indsL,:),[-2 2]); hold on;
plot([2000 2000],[0 200],'k-')

subplot(4,3,[1 4 7]+2)
imagesc(abs( analysis.approach.left.pAvg(indsL,:)-analysis.approach.right.pZav(indsL,:) ) , [-2 2]) ; hold on;
plot([2000 2000],[0 200],'k-')


subplot(4,3,10+0)
plot(mean(analysis.approach.left.pAvg),'b'); hold on;
plot([2000 2000],[-2 2],'k-');
axis([0 7000 -0.3 0.2]);

subplot(4,3,10+1)
plot(mean(analysis.approach.right.pAvg),'b'); hold on;
plot([2000 2000],[-2 2],'k-');
axis([0 7000 -0.3 0.2]);

subplot(4,3,10+2); hold off;
plot(mean( abs( analysis.approach.left.pAvg(indsL,:)-analysis.approach.right.pAvg(indsL,:) ) ),'b'); hold on;
plot([2000 2000],[-2 2],'k-');
axis([0 7000 0 0.75]);


% Clock representation
[mappedA, mapping] = compute_mapping(analysis.approach.left.pAvg','PCA',3);
[mappedB, mappingB] = compute_mapping(analysis.approach.right.pAvg','PCA',3);
% 
% for i=1:7001
%     for mm = 1:3
%         mappedB(i,mm) = dot( analysis.approach.right.pAvg(:,i) , mapping.M(:,mm) );
%     end
% end
% 
% theta   = (mappedA(:,1) - min(mappedA(:,1))) ./ ( max(mappedA(:,1)) - min(mappedA(:,1)) ) .* 2 .* pi;
% rho     = (mappedA(:,2) - min(mappedA(:,2))) ./ ( max(mappedA(:,2)) - min(mappedA(:,2)) );
% 
% figure(1); clf; polar(theta,rho,'k-'); hold on; polar(theta(2000:5000),rho(2000:5000),'b-'); 

figure(2); clf; scatter(mappedA(:,2),mappedA(:,3),2,[1:7001]); hold on; scatter(mappedB(:,2),mappedB(:,3),10,[1:7001]);
xlabel('PC2');ylabel('PC3');
figure(3); 
for mm=1:3
    subplot(3,1,mm);
    [val,inds] = sort(mappingB.M(:,mm));
    plot(1:222 , mappingB.M(inds,mm) , 'ko' , 1:222 , mapping.M(inds,mm) , 'ro');
    title(num2str(dot(mappingB.M(:,mm),mapping.M(:,mm))));
end

figure(4);
plot(-2000:5000,mappedA(:,1));

%% ALL SESSIONS >>> LEAVE SEP BY SIDE

    tmpPrior = currParams.winParams.prior;
    currParams.winParams.prior = 4000;
    tmpAfter = currParams.winParams.after;
    currParams.winParams.after = 2000;

a = 1;
currSes = 1;
ySpc = 4;
count=1;

% ALIGN ON PRESS (right VS RIGHT)
for a=1:2

    for currSes = 1:size(PData.mouse(a).session,2)

        disp(['Extracting data from session ' num2str(currSes) ' of ' num2str(size(PData.mouse(a).session,2)) ' total.']);

        for j = 1:PData.mouse(a).session(currSes).numUnits

            numStamps   = numel(PData.mouse(a).session(currSes).unit(j).ts);
            delta       = zeros(1,ceil(PData.mouse(a).session(currSes).unit(j).ts(numStamps)));
            delta(PData.mouse(a).session(currSes).unit(j).ts) = 1;
            tmpSmooth   = conv(delta,currParams.filter.kernel,'same');

            leftInds = find(PData.mouse(a).session(currSes).Behavior.LeverSide==1 & PData.mouse(a).session(currSes).Behavior.TrialType==2);
            [rsLevPress] = TNC_AlignRasters(tmpSmooth , PData.mouse(a).session(currSes).unit(j).ts , -1 , round( PData.mouse(a).session(currSes).Behavior.Leave(leftInds)-currParams.offset ), [currParams.winParams.prior currParams.winParams.after],0,1);
            analysis.leave.left.pAvg(count,:) = rsLevPress.image.psthAVG;
            analysis.leave.left.pZav(count,:) = rsLevPress.image.psthZ;
            analysis.leave.left.sess(count) = currSes;
            analysis.leave.left.unit(count) = j;
            analysis.leave.left.mouse(count) = a;

            rightInds = find(PData.mouse(a).session(currSes).Behavior.LeverSide==0 & PData.mouse(a).session(currSes).Behavior.TrialType==2);
            [rsLevPress] = TNC_AlignRasters(tmpSmooth , PData.mouse(a).session(currSes).unit(j).ts , -1 , round( PData.mouse(a).session(currSes).Behavior.Leave(rightInds)-currParams.offset ), [currParams.winParams.prior currParams.winParams.after],0,1);
            analysis.leave.right.pAvg(count,:) = rsLevPress.image.psthAVG;
            analysis.leave.right.pZav(count,:) = rsLevPress.image.psthZ;
            analysis.leave.right.sess(count) = currSes;
            analysis.leave.right.unit(count) = j;
            analysis.leave.right.mouse(count) = a;

            count = count+1;

        end
    end
end


clear currData*

for j = 1:size(analysis.leave.right.pZav,1)
    currDataR = sgolayfilt( cumsum( analysis.leave.right.pZav(j,2500:5000) ),   3,  21);
    currDataL = sgolayfilt( cumsum( analysis.leave.left.pZav(j,2500:5000) ),   3,  21);
    
    peakTimeR(j) = find(diff(currDataR)==max(diff(currDataR)),1,'first');
    peakTimeL(j) = find(diff(currDataL)==max(diff(currDataL)),1,'first');
end

[vals,indsR] = sort(peakTimeR);
[vals,indsL] = sort(peakTimeL);

figure(12); clf;
[mapName] = TNC_CreateRBColormap(1024,'rb');

subplot(4,3,[1 4 7])
imagesc(-currParams.winParams.prior:currParams.winParams.after,1:200,analysis.leave.left.pZav(indsL,:),[-2 2]); hold on;
plot([0 0],[0 200],'k-')
set(gca,'TickDir','out');
ylabel('Cell Index');

subplot(4,3,[1 4 7]+1)
imagesc(-currParams.winParams.prior:currParams.winParams.after,1:200,analysis.leave.right.pZav(indsR,:),[-2 2]); hold on;
plot([0 0],[0 200],'k-')
set(gca,'TickDir','out');

subplot(4,3,[1 4 7]+2)
imagesc(-currParams.winParams.prior:currParams.winParams.after,1:200,abs(analysis.leave.left.pZav(indsL,:)-analysis.leave.right.pZav(indsL,:)),[0 3]); hold on;
plot([0 0],[0 200],'k-')
set(gca,'TickDir','out');

colormap(mapName);

subplot(4,3,10+0)
plot(-currParams.winParams.prior:currParams.winParams.after, mean(analysis.leave.left.pZav),'k','LineWidth',2);
hold on; plot([0 0],[-5 5],'k-');
axis([-currParams.winParams.prior currParams.winParams.after -0.15 0.4]);
set(gca,'TickDir','out');
box off;
ylabel('Z-score');
xlabel('Time relative to leaving (ms)');

subplot(4,3,10+1)
plot(-currParams.winParams.prior:currParams.winParams.after, mean(analysis.leave.right.pZav),'k','LineWidth',2);
hold on; plot([0 0],[-5 5],'k-');
axis([-currParams.winParams.prior currParams.winParams.after -0.15 0.4]);
set(gca,'TickDir','out');
box off;

subplot(4,3,10+2)
plot(-currParams.winParams.prior:currParams.winParams.after, mean(abs(analysis.leave.left.pZav-analysis.leave.right.pZav)),'k','LineWidth',2); 
hold on; plot([0 0],[-5 5],'k-');
axis([-currParams.winParams.prior currParams.winParams.after 0 0.5]);
set(gca,'TickDir','out');
box off;

    
    currParams.winParams.prior = tmpPrior;
    currParams.winParams.after = tmpAfter;

%% ALL SESSIONS >>> approach SEP BY SIGMA AND SIDE

currSes = 1;
ySpc = 4;
count=1;

% ALIGN ON PRESS
for a=1:2

    for currSes = 1:size(PData.mouse(a).session,2)

        disp(['Extracting data from session ' num2str(currSes) ' of ' num2str(size(PData.mouse(a).session,2)) ' total.']);

        for j = 1:PData.mouse(a).session(currSes).numUnits

            numStamps   = numel(PData.mouse(a).session(currSes).unit(j).ts);
            delta       = zeros(1,ceil(PData.mouse(a).session(currSes).unit(j).ts(numStamps)));
            delta(PData.mouse(a).session(currSes).unit(j).ts) = 1;
            tmpSmooth   = conv(delta,currParams.filter.kernel,'same');

            leftInds = find(PData.mouse(a).session(currSes).Behavior.LeverSide==1 & PData.mouse(a).session(currSes).Behavior.TrialType==2);
            [rsLevPress] = TNC_AlignRasters(tmpSmooth , PData.mouse(a).session(currSes).unit(j).ts , -1 ,round( PData.mouse(a).session(currSes).Behavior.Approach(leftInds)-currParams.offset ), [currParams.winParams.prior currParams.winParams.after],0,1);
            analysis.approach.left.pAvg(count,:) = rsLevPress.image.psthAVG;
            analysis.approach.left.pZav(count,:) = rsLevPress.image.psthZ;
            analysis.approach.left.sess(count) = currSes;
            analysis.approach.left.unit(count) = j;
            analysis.approach.left.mouse(count) = a;

            rightInds = find(PData.mouse(a).session(currSes).Behavior.LeverSide==0 & PData.mouse(a).session(currSes).Behavior.TrialType==2);
            [rsLevPress] = TNC_AlignRasters(tmpSmooth , PData.mouse(a).session(currSes).unit(j).ts , -1 ,round( PData.mouse(a).session(currSes).Behavior.Approach(rightInds)-currParams.offset ), [currParams.winParams.prior currParams.winParams.after],0,1);
            analysis.approach.right.pAvg(count,:) = rsLevPress.image.psthAVG;
            analysis.approach.right.pZav(count,:) = rsLevPress.image.psthZ;
            analysis.approach.right.sess(count) = currSes;
            analysis.approach.right.unit(count) = j;
            analysis.approach.right.mouse(count) = a;

            PData.mouse(a).session(currSes).unit(j).meanRate = rsLevPress.meanRate;
            PData.mouse(a).session(currSes).unit(j).stdRate = rsLevPress.stdRate;

            count = count+1;

        end
    end
end

clear currData*

for j = 1:size(analysis.approach.right.pZav,1)
    currDataR(j) = trapz( analysis.approach.right.pZav(j,2000:6000) );
    currDataL(j) = trapz( analysis.approach.left.pZav(j,2000:6000) );
end

[vals,indsR] = sort(currDataR,'descend');
[vals,indsL] = sort(currDataL,'descend');

% indsR = 1:size(analysis.approach.right.pZav,1);
% indsL = 1:size(analysis.approach.right.pZav,1);

figure(9); clf;
[mapName] = TNC_CreateRBColormap(1024,'rb');

subplot(4,3,[1 4 7])
imagesc(-currParams.winParams.prior:currParams.winParams.after,1:200,analysis.approach.left.pZav(indsL,:),[-2 2]); hold on;
plot([0 0],[0 200],'k-')
colormap(mapName);
set(gca,'TickDir','out');
ylabel('Cell Index');

subplot(4,3,[1 4 7]+1)
imagesc(-currParams.winParams.prior:currParams.winParams.after,1:200,analysis.approach.right.pZav(indsR,:),[-2 2]); hold on;
plot([0 0],[0 200],'k-')
set(gca,'TickDir','out');

subplot(4,3,[1 4 7]+2)
imagesc(-currParams.winParams.prior:currParams.winParams.after,1:200,abs(analysis.approach.left.pZav(indsL,:)-analysis.approach.right.pZav(indsL,:))); hold on;
plot([0 0],[0 200],'k-')
set(gca,'TickDir','out');

subplot(4,3,10+0)
plot(-currParams.winParams.prior:currParams.winParams.after,mean(analysis.approach.left.pZav),'k', 'LineWidth', 2); hold on;
plot([0 0],[-2 2],'k--');
axis([-currParams.winParams.prior currParams.winParams.after -0.3 0.2]);
set(gca,'TickDir','out'); box off;
ylabel('Z-score');
xlabel('Time relative to entry (ms)');

subplot(4,3,10+1)
plot(-currParams.winParams.prior:currParams.winParams.after,mean(analysis.approach.right.pZav),'k', 'LineWidth', 2); hold on;
plot([0 0],[-2 2],'k--');
set(gca,'TickDir','out'); box off;
axis([-currParams.winParams.prior currParams.winParams.after -0.3 0.2]);

subplot(4,3,10+2)
plot(-currParams.winParams.prior:currParams.winParams.after,mean(abs(analysis.approach.left.pZav-analysis.approach.right.pZav)),'k', 'LineWidth', 2); hold on;
plot([0 0],[-2 2],'k--');
set(gca,'TickDir','out'); box off;
axis([-currParams.winParams.prior currParams.winParams.after 0 0.75]);

%% ALL SESSIONS >>> GET APPROACH PROJECTION

figure(20); clf; subplot(121); imagesc(corr(analysis.approach.left.pZav(indsL,:)'), [-1 1]); colormap(mapName);
    
covMatForTimeCourse = cov(analysis.approach.left.pZav');
[v,d] = eig(covMatForTimeCourse);
analysis.approach.pca.vecs = zeros(size(analysis.approach.left.pZav,1),3);
analysis.approach.pca.vals = diag(d)./sum(diag(d));
analysis.approach.pca.vecs(:,1) = v(:,numel(analysis.approach.pca.vals));
analysis.approach.pca.vecs(:,2) = v(:,numel(analysis.approach.pca.vals)-1);
analysis.approach.pca.vecs(:,3) = v(:,numel(analysis.approach.pca.vals)-2);
analysis.approach.pca.varExp = sum(analysis.approach.pca.vals(numel(analysis.approach.pca.vals)-2:numel(analysis.approach.pca.vals)));
disp(['Variance explained: ' num2str(sum(analysis.approach.pca.vals(numel(analysis.approach.pca.vals)-2:numel(analysis.approach.pca.vals)))) ' | ' num2str(analysis.approach.pca.vals(numel(analysis.approach.pca.vals))) ' ' num2str(analysis.approach.pca.vals(numel(analysis.approach.pca.vals)-1)) ' ' num2str(analysis.approach.pca.vals(numel(analysis.approach.pca.vals)-2)) ]);
disp(' ');

for i=1:size(analysis.approach.left.pZav,2)
    
    % project onto leading dimensions
    analysis.approach.pca.lTraj(1,i) = dot(analysis.approach.left.pZav(:,i),analysis.approach.pca.vecs(:,1));
    analysis.approach.pca.lTraj(2,i) = dot(analysis.approach.left.pZav(:,i),analysis.approach.pca.vecs(:,2));
    analysis.approach.pca.lTraj(3,i) = dot(analysis.approach.left.pZav(:,i),analysis.approach.pca.vecs(:,3));

    % project onto leading dimensions
    analysis.approach.pca.rTraj(1,i) = dot(analysis.approach.right.pZav(:,i),analysis.approach.pca.vecs(:,1));
    analysis.approach.pca.rTraj(2,i) = dot(analysis.approach.right.pZav(:,i),analysis.approach.pca.vecs(:,2));
    analysis.approach.pca.rTraj(3,i) = dot(analysis.approach.right.pZav(:,i),analysis.approach.pca.vecs(:,3));

end

figure(20); clf;
plot3(analysis.approach.pca.lTraj(1,:),analysis.approach.pca.lTraj(2,:),analysis.approach.pca.lTraj(3,:),'Color',[0 0.67 1],'LineWidth',2); hold on;
plot3(analysis.approach.pca.rTraj(1,:),analysis.approach.pca.rTraj(2,:),analysis.approach.pca.rTraj(3,:),'r','LineWidth',2);
plot3(analysis.approach.pca.lTraj(1,1),analysis.approach.pca.lTraj(2,1),analysis.approach.pca.lTraj(3,1),'s','LineWidth',5,'Color',[0 0.67 1]);
plot3(analysis.approach.pca.rTraj(1,1),analysis.approach.pca.rTraj(2,1),analysis.approach.pca.rTraj(3,1),'s','LineWidth',5,'Color',[1 0 0]);
grid on;
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

% c=1:numel(analysis.approach.pca.rTraj(1,:));
% scatter3(analysis.approach.pca.lTraj(1,:),analysis.approach.pca.lTraj(2,:),analysis.approach.pca.lTraj(3,:),5,c); hold on;
% scatter3(analysis.approach.pca.rTraj(1,:),analysis.approach.pca.rTraj(2,:),analysis.approach.pca.rTraj(3,:),5,c);
% colormap(jet);

%% CURR SESSION >>> CONVERT DATASTRUCTURE TO THE PopData.session(N).unit(N) formalism
count       = 1;
tmp         = find(opSigPop.Nev.Data.Spikes.Unit>0);
allTrodes   = unique(opSigPop.Nev.Data.Spikes.Electrode(tmp));

for i=1:numel(allTrodes)

    currTrode = allTrodes(i);

    for j=1:6

        spkInds = find(opSigPop.Nev.Data.Spikes.Electrode==currTrode & opSigPop.Nev.Data.Spikes.Unit==j);

        if numel(spkInds)>0
            
            % index through channels looking for sorted units
            opSigPop.session(currSes).unit(count).ts   = round(double(opSigPop.Nev.Data.Spikes.Timestamps(spkInds))./30);
            opSigPop.session(currSes).unit(count).el   = currTrode;
            opSigPop.session(currSes).unit(count).un   = j;

            % Get the electrode map
            %[row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(currTrode,'NN_w64');
            %PopData.session(currSes).unit(count).row = row;
            %PopData.session(currSes).unit(count).col = col;

            count = count+1;

        end
    end
end

count = count - 1;

opSigPop.session(currSes).numUnits = count;
clear count;

%% CURR SESSION >>> LOAD CONTINUOUS LEVER DATA

    sLeverData(1,:) = sgolayfilt(decimate(opSigPop.continuousData.Data(1,:),10),9,101);
    sLeverData(2,:) = sgolayfilt(decimate(opSigPop.continuousData.Data(2,:),10),9,101);

    difLick = diff(decimate(opSigPop.continuousData.Data(3,:),10));
    opSigPop.Behavior.eLick   = find(difLick>1000);

    % EXTRACT VELOCITY DATA FROM CONT LEVER DATA
    numSamples = size(sLeverData,2);
    tmpLeverData(1,:) = sgolayfilt(sLeverData(1,:),3,201);
    tmpLeverData(2,:) = sgolayfilt(sLeverData(2,:),3,201);

    disp(' '); disp(' '); disp('Extracting velocity...');
    dLeft = diff(tmpLeverData(1,:));
    dRight = diff(tmpLeverData(2,:));
    disp(' '); disp(' Complete. '); disp(' ');
    
%% CURR SESSION >>> FIND ALL LEVER PRESSES

    clear progSt* reach

    pre         = 10;
    post        = 10;       
    minSpace    = 250;
    countL      = 1;
    countR      = 1;

    % threshold the velocities
    vLeft       = find(dLeft(1:numel(dLeft)-(minSpace.*2))<-1);
    vRight      = find(dRight(1:numel(dLeft)-(minSpace.*2))<-1);

    progStart.left(countL)   = vLeft(1);
    progStart.right(countR)  = vRight(1);

    for j=2:numel(vLeft)

        if vLeft(j)>vLeft(j-1)+minSpace

            postN = find(dLeft(vLeft(j-1):vLeft(j))>0,1,'first');
            if numel(postN)<1
                disp('Cannot find stop');                        
                postN=post;
            end

            progStop.left(countL)   = vLeft(j-1)+postN;
            countL                  = countL+1;

            preN = find(dLeft(vLeft(j-1) : vLeft(j))<1,1,'last');

            if numel(preN)<1
                disp('Cannot find start');
                progStart.left(countL)  = vLeft(j)-pre;                    
            else
                progStart.left(countL)  = vLeft(j-1)+preN;
            end

        end

        if j==numel(vLeft)
            progStop.left(countL)   = vLeft(j)+post;
        end

    end
    
    for j=2:numel(vRight)

        if vRight(j)>vRight(j-1)+minSpace

            postN = find(dRight(vRight(j-1):vRight(j))>0,1,'first');
            if numel(postN)<1
                disp('Cannot find stop');                        
                postN=post;
            end

            progStop.right(countR)  = vRight(j-1)+postN;
            countR                  = countR+1;

            preN = find(dRight(vRight(j-1) : vRight(j))<1,1,'last');

            if numel(preN)<1
                disp('Cannot find start');
                progStart.right(countR)  = vRight(j)-pre;                    
            else
                progStart.right(countR)  = vRight(j-1)+preN;
            end

        end

        if j==numel(vRight)
            progStop.right(countR)   = vRight(j)+post;
        end

    end
    
    figure(4); clf;
    
%     subplot(211)
    plot(dLeft,'r'); hold on;
    plot(progStart.left,ones(1,numel(progStart.left)).*-2,'go');
    plot(progStop.left,zeros(1,numel(progStop.left)),'ms');

%     subplot(212)
    plot(dRight,'k'); hold on;
    plot(progStart.right,ones(1,numel(progStart.right)).*-2,'ro');
    plot(progStop.right,zeros(1,numel(progStop.right)),'bs');
    
%% CURR SESSION >>> EXTRACT PROPERTIES OF EACH LEVER PRESS
    clear reach;
    
    countL      = 1;
    countR      = 1;
    
    for k = 1:numel(progStart.left)

        if k==1
            reach.init = 1;
        end

        % reaches must be at least 60 ms long
        if progStop.left(k)-progStart.left(k)>60 & progStart.left(k)>minSpace

            velTraj     = dLeft(progStart.left(k) : progStop.left(k));
            xVals       = tmpLeverData(1,progStart.left(k) : progStop.left(k));

            reach.left.start(countL)  = progStart.left(k);
            reach.left.stop(countL)   = progStop.left(k);
            reach.left.dist(countL)   = trapz(velTraj);

            reach.left.dur(countL,1)    = progStop.left(k) - progStart.left(k);
            reach.left.vel(countL,1)    = max(velTraj);
            reach.left.vel(countL,2)    = trapz(velTraj) ./ reach.left.dur(countL,1);
            reach.left.vel(countL,3)    = var(velTraj);

            reach.left.acc(countL,1)   = max(diff(velTraj));
            reach.left.acc(countL,2)   = mean(diff(velTraj));
            reach.left.acc(countL,3)   = max(diff(velTraj(1:50))); % max in first 50 ms of movement

            countL                 = countL+1;

        end            
    end
    
    
    for k = 1:numel(progStart.right)

        if k==1
            reach.init = 1;
        end

        % reaches must be at least 60 ms long
        if progStop.right(k)-progStart.right(k)>60 & progStart.right(k)>minSpace

            velTraj     = dRight(progStart.right(k) : progStop.right(k));
            xVals       = tmpLeverData(1,progStart.right(k) : progStop.right(k));

            reach.right.start(countR)  = progStart.right(k);
            reach.right.stop(countR)   = progStop.right(k);
            reach.right.dist(countR)   = trapz(velTraj);

            reach.right.dur(countR,1)    = progStop.right(k) - progStart.right(k);
            reach.right.vel(countR,1)    = max(velTraj);
            reach.right.vel(countR,2)    = trapz(velTraj) ./ reach.right.dur(countR,1);
            reach.right.vel(countR,3)    = var(velTraj);

            reach.right.acc(countR,1)   = max(diff(velTraj));
            reach.right.acc(countR,2)   = mean(diff(velTraj));
            reach.right.acc(countR,3)   = max(diff(velTraj(1:50))); % max in first 50 ms of movement

            countR                 = countR+1;

        end            
    end
    
    reach.left.numReaches   = countL-1;
    reach.right.numReaches  = countR-1;
    
    opSigPop.Behavior.reach = reach;
    
%% CURR SESSION >>> FOR EVERY VALID TRIAL RETURN PROPERTIES OF THE LEVER PRESS (NUMBER, VELOCITY, DURATION)

trialLinds = find(opSigPop.Behavior.LeverSide==1);
trialRinds = find(opSigPop.Behavior.LeverSide==0);
clear TrigPress

for i = 1:numel(trialLinds)
    
    currTrialStart                                  = opSigPop.Behavior.CueLight(trialLinds(i));
    
    if trialLinds(i)>1
        lastTrialEnd    = opSigPop.Behavior.Leave(trialLinds(i)-1);
    else
        lastTrialEnd    = 0;
    end
    trigReachInd                                    = find(reach.left.start<currTrialStart,1,'last');
    firstReachInd                                   = find(reach.left.start<currTrialStart & reach.left.start>lastTrialEnd,1,'first');

    TrigPress(trialLinds(i))      = reach.left.start(trigReachInd);
    LeverVel(trialLinds(i),:)     = reach.left.vel(trigReachInd,:);
    % opSigPop.Behavior.FirstPress(trialLinds(i))     = reach.left.start(firstReachInd);
    % opSigPop.Behavior.NumPress(trialLinds(i))       = numel(find(reach.left.start<currTrialStart & reach.left.start>lastTrialEnd));

    if opSigPop.Behavior.CueLight(trialLinds(i))-TrigPress(trialLinds(i))>1000
        figure(10);
        plot(dLeft(1,opSigPop.Behavior.CueLight(trialLinds(i))-1000:opSigPop.Behavior.CueLight(trialLinds(i))+1000));
        title(num2str(opSigPop.Behavior.CueLight(trialLinds(i))-TrigPress(trialLinds(i))));
        pause();
    end
end

for i = 1:numel(trialRinds)
    
    currTrialStart                                  = opSigPop.Behavior.CueLight(trialRinds(i));
    if trialRinds(i)>1
        lastTrialEnd    = opSigPop.Behavior.Leave(trialRinds(i)-1);
    else
        lastTrialEnd    = 0;
    end
    trigReachInd                                    = find(reach.right.start<currTrialStart,1,'last');
    firstReachInd                                   = find(reach.right.start<currTrialStart & reach.right.start>lastTrialEnd,1,'first');

    TrigPress(trialRinds(i))      = reach.right.start(trigReachInd);
    LeverVel(trialRinds(i),:)     = reach.right.vel(trigReachInd,:);
    % opSigPop.Behavior.FirstPress(trialLinds(i))     = reach.left.start(firstReachInd);
    % opSigPop.Behavior.NumPress(trialLinds(i))       = numel(find(reach.left.start<currTrialStart & reach.left.start>lastTrialEnd));

end

opSigPop.Behavior.TrigPress = TrigPress;
opSigPop.Behavior.LeverVel = LeverVel;
    
%% CURR SESSION >>> DECISION TO PRESS LEVER
currSes = 1;
ySpc = 4;
clear analysis;

% ALIGN ON PRESS (right VS RIGHT)
    for j = 1:opSigPop.session(currSes).numUnits

        numStamps   = numel(opSigPop.session(currSes).unit(j).ts);
        delta       = zeros(1,ceil(opSigPop.session(currSes).unit(j).ts(numStamps)));
        delta(opSigPop.session(currSes).unit(j).ts) = 1;
        tmpSmooth   = conv(delta,currParams.filter.kernel,'same');
        
        leftInds = find(opSigPop.Behavior.LeverSide==1 & opSigPop.Behavior.TrialType==2);
        [rsLevPress] = TNC_AlignRasters(tmpSmooth , opSigPop.session(currSes).unit(j).ts , -1 , opSigPop.Behavior.Approach(leftInds)-currParams.offset, [currParams.winParams.prior currParams.winParams.after],0,1);
        analysis.alLeverPress.left.pAvg(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthAVG;
        analysis.alLeverPress.left.pZav(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthZ;

        rightInds = find(opSigPop.Behavior.LeverSide==0 & opSigPop.Behavior.TrialType==2);
        [rsLevPress] = TNC_AlignRasters(tmpSmooth , opSigPop.session(currSes).unit(j).ts , -1 , opSigPop.Behavior.Approach(rightInds)-currParams.offset, [currParams.winParams.prior currParams.winParams.after],0,1);
        analysis.alLeverPress.right.pAvg(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthAVG;
        analysis.alLeverPress.right.pZav(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthZ;
        
        opSigPop.session(currSes).unit(j).meanRate = rsLevPress.meanRate;
        opSigPop.session(currSes).unit(j).stdRate = rsLevPress.stdRate;
        
    end

    figure(1); clf;
    % subplot(121);
    % imagesc(analysis.alLeverPress.pZav,[-2 2]); hold on;
    % set(gca,'YDir','normal');
    % [mapName] = TNC_CreateRBColormap(1024,'rb');
    % colormap(mapName);
    % plot([currParams.winParams.prior currParams.winParams.prior],[0 opSigPop.session(currSes).numUnits+1],'k--');


    % subplot(122);
    for j = 1:opSigPop.session(currSes).numUnits
        subplot(5,3,j); ySpc=0;
        plot(-currParams.winParams.prior:currParams.winParams.after,(analysis.alLeverPress.left.pAvg(j,:).*1000) + (ySpc.*j),'Color',[1 0 0],'LineWidth',2); hold on;
        plot(-currParams.winParams.prior:currParams.winParams.after,(analysis.alLeverPress.right.pAvg(j,:).*1000) + (ySpc.*j),'Color',[0 0.67 1],'LineWidth',2); hold on;
%         text(50-currParams.winParams.prior,(ySpc.*j)+2,num2str(opSigPop.session(currSes).unit(j).meanRate.*1000));
        plot([0 0],[min(analysis.alLeverPress.left.pAvg(j,:).*1000) max(analysis.alLeverPress.left.pAvg(j,:).*1000)],'k-');
        plot([0 0],[min(analysis.alLeverPress.right.pAvg(j,:).*1000) max(analysis.alLeverPress.right.pAvg(j,:).*1000)],'k-');
        set(gca,'TickDir','out');
    end

%     plot([0 0],[0 ySpc.*(opSigPop.session(currSes).numUnits+1)],'k--');

% WARP ALIGN BETWEEN LEAVE AND PRESS

%% CURR SESSION >>> DECISION TO LEAVE PORT

% COMPARE LAST LICK > LEAVE INTERVAL with the LEAVE ON PROBES

%% CURR SESSION >>> ALIGN TO LEAVE, SORTED BY SIGMArew

    tmpPrior = currParams.winParams.prior;
    currParams.winParams.prior = 3000;
    tmpAfter = currParams.winParams.after;
    currParams.winParams.after = 2000;
    
    figure(2); clf;
    ySpc = 3;

    for j = 1:opSigPop.session(currSes).numUnits

        numStamps   = numel(opSigPop.session(currSes).unit(j).ts);
        delta       = zeros(1,ceil(opSigPop.session(currSes).unit(j).ts(numStamps)));
        delta(opSigPop.session(currSes).unit(j).ts) = 1;
        tmpSmooth   = conv(delta,currParams.filter.kernel,'same');

        sigsUsed = unique(opSigPop.Behavior.SigmaR);
        numSigs = numel(sigsUsed);
        
        for p=1:numSigs
            currInds = find(opSigPop.Behavior.SigmaR==sigsUsed(p) & opSigPop.Behavior.TrialType==1 & opSigPop.Behavior.EarlyLeaveInd==0);
            [rsLevPress] = TNC_AlignRasters(tmpSmooth , opSigPop.session(currSes).unit(j).ts , -1 , opSigPop.Behavior.Leave(currInds)-currParams.offset, [currParams.winParams.prior currParams.winParams.after],0,1);
            if p==1
                plot(-currParams.winParams.prior:currParams.winParams.after,zeros(1,numel(rsLevPress.image.psthAVG))+(ySpc.*j),'-','Color',[0.5 0.5 0.5],'LineWidth',1); hold on;
%                 text(50-currParams.winParams.prior,(ySpc.*j)+(ySpc./2),num2str(opSigPop.session(currSes).unit(j).meanRate.*1000),'FontWeight','bold');
                if j==1
                    analysis.leave.sigr(p).pAvg = zeros(opSigPop.session(currSes).numUnits,numel(rsLevPress.image.psthAVG));
                    analysis.leave.sigr(p).pZav = zeros(opSigPop.session(currSes).numUnits,numel(rsLevPress.image.psthAVG));
                end
            end
            analysis.leave.sigr(p).pAvg(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthAVG;
            analysis.leave.sigr(p).pZav(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthZ;
            analysis.leave.sigr(p).sigma                                     = sigsUsed(p);
            analysis.leave.sigr(p).side                                      = mode(opSigPop.Behavior.LeverSide(currInds));

            subplot(5,3,j); ySpc=0;
            plot(-currParams.winParams.prior:currParams.winParams.after,(analysis.leave.sigr(p).pZav(j,:).*1) + (ySpc.*j),'Color',[(p-1)./(numSigs-1) 0.67-(0.67.*((p-1)./(numSigs-1))) 1-((p-1)./(numSigs-1))],'LineWidth',2); hold on;                        
            plot([0 0],[0 ySpc.*(opSigPop.session(currSes).numUnits+1)],'k--'); box off;
            set(gca,'TickDir','out');
%             text((1000.*p)-currParams.winParams.prior,(ySpc./2),num2str(analysis.leave.sigr(p).sigma),'Color',[(p-1)./(numSigs-1) 0.67-(0.67.*((p-1)./(numSigs-1))) 1-((p-1)./(numSigs-1))],'FontWeight','bold');
            
        end
    end

    set(gca,'TickDir','out');
    plot([0 0],[0 ySpc.*(opSigPop.session(currSes).numUnits+1)],'k--');
    
    currParams.winParams.prior = tmpPrior;
    currParams.winParams.after = tmpAfter;
    
%     figure(3);
%     for p=1:3
%         subplot(1,3,p);
%         imagesc(analysis.leave.sigr(p).pZav,[-1 1]);
%         [mapName] = TNC_CreateRBColormap(1024,'rb');
%         colormap(mapName);
%     end
%     

%% CURR SESSION >>> COMPARE PROBE LEAVE AND INTRABLOCK ERROR LEAVE

    tmpPrior = currParams.winParams.prior;
    currParams.winParams.prior = 5000;
    tmpAfter = currParams.winParams.after;
    currParams.winParams.after = 5000;
    
    figure(3); clf;
    ySpc = 3;

    for j = 1:opSigPop.session(currSes).numUnits

        numStamps   = numel(opSigPop.session(currSes).unit(j).ts);
        delta       = zeros(1,ceil(opSigPop.session(currSes).unit(j).ts(numStamps)));
        delta(opSigPop.session(currSes).unit(j).ts) = 1;
        tmpSmooth   = conv(delta,currParams.filter.kernel,'same');

        sigsUsed = unique(opSigPop.Behavior.SigmaR);
        numSigs = numel(sigsUsed);
        
        currIndsEtmp = find(opSigPop.Behavior.SigmaR==sigsUsed(1) & opSigPop.Behavior.IntraBlockErr==1);
%             currIndsE = currIndsEtmp(find(diff(currIndsEtmp)==1));
            currIndsE = currIndsEtmp;
        currIndsP = find(opSigPop.Behavior.SigmaR==sigsUsed(3) & opSigPop.Behavior.TrialType==2);
        currIndsR = find(opSigPop.Behavior.SigmaR==sigsUsed(3) & opSigPop.Behavior.TrialType==1);
            
        [rsLevPress] = TNC_AlignRasters(tmpSmooth , opSigPop.session(currSes).unit(j).ts , -1 , opSigPop.Behavior.Leave(currIndsE)-currParams.offset, [currParams.winParams.prior currParams.winParams.after],0,1);

        analysis.leave.PrErr(1).pAvg(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthAVG;
        analysis.leave.PrErr(1).pZav(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthZ;

        [rsLevPress] = TNC_AlignRasters(tmpSmooth , opSigPop.session(currSes).unit(j).ts , -1 , opSigPop.Behavior.Leave(currIndsP)-currParams.offset, [currParams.winParams.prior currParams.winParams.after],0,1);

        analysis.leave.PrErr(2).pAvg(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthAVG;
        analysis.leave.PrErr(2).pZav(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthZ;

        [rsLevPress] = TNC_AlignRasters(tmpSmooth , opSigPop.session(currSes).unit(j).ts , -1 , opSigPop.Behavior.Leave(currIndsR)-currParams.offset, [currParams.winParams.prior currParams.winParams.after],0,1);

        analysis.leave.PrErr(3).pAvg(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthAVG;
        analysis.leave.PrErr(3).pZav(j,1:numel(rsLevPress.image.psthAVG)) = rsLevPress.image.psthZ;
        
        plot(-currParams.winParams.prior:currParams.winParams.after,zeros(1,numel(rsLevPress.image.psthAVG))+(ySpc.*j),'-','Color',[0.5 0.5 0.5],'LineWidth',1); hold on;
        plot(-currParams.winParams.prior:currParams.winParams.after,(analysis.leave.PrErr(1).pZav(j,:).*1) + (ySpc.*j),'Color',[0 0.67 1],'LineWidth',2);
        plot(-currParams.winParams.prior:currParams.winParams.after,(analysis.leave.PrErr(2).pZav(j,:).*1) + (ySpc.*j),'Color',[1 0 0],'LineWidth',2);                  
        plot(-currParams.winParams.prior:currParams.winParams.after,(analysis.leave.PrErr(3).pZav(j,:).*1) + (ySpc.*j),'Color',[0.5 0.33 0.5],'LineWidth',2);                  
            
    end

    set(gca,'TickDir','out');
    plot([0 0],[0 ySpc.*(opSigPop.session(currSes).numUnits+1)],'k--');
        
    currParams.winParams.prior = tmpPrior;
    currParams.winParams.after = tmpAfter;
    
% EXAMINE THE LEAVE TRANSITIONS DURING INTERBLOCK ERRORS

%% CURR SESSION >>> TRANSITION TO EXPLORATION

%% COMPILED DATA >>> IDENTIFY NUMBER OF CLUSTERS, CALCULATE MEAN PSTHs FOR EACH CLUSTER

% Find the number of clusters

Y = 3:0.2:5;

for i = 1:numel(Y);

   for k = 2:15; %define possible sizes of k

       Color = {'m','b','g','r','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colors.

       k
       
       %Kmeans requires rows to contain observations (so cells), columns to contain variables
       FullChar_2 = Probe.Block2000.pZav;
%        FullChar_2 = Probe.Block50.pZav;
%        FullChar_2 = Probe.Block750.pZav;
% FullChar_2 = analysis.approach.right.pAvg;

       [idx,C,sumd,D] = kmeans(FullChar_2,k,'replicates',3);


       %find distance from closest cluster centroid

       for j = 1:size(FullChar_2);

           ClostestToCentroid(j,1) = D(j,idx(j,:));

       end

       %compute the distortion - i.e. the normalised square distance between each
       %observation and its clostest cluster center.

       normD = ClostestToCentroid/max(ClostestToCentroid);

       meanDistortion = mean(normD.^2);

       valk(k,1) = meanDistortion;

       YmeanDistortion(i,k) = meanDistortion.^Y(1,i);

   end

end

figure1 = figure('Color',[1 1 1]);
set(gca, 'box','off');


for i = 1:numel(Y);

   [peak,idx] = max(YmeanDistortion(i,:));

   subplot(3,4,i);
   plot(YmeanDistortion(i,:));
   hold on
   plot(idx,peak,'ro');

   xlabel ('Number of Clusters')
   ylabel ('Jump');

   str = sprintf('Y = %f',Y(1,i));

   title(str);
   hold on

   peakKVal(i,1) = idx;

end

%% Based upon the median peakKVal display the clustered correlation matrix
% and the mean PSTHs
% numClusts = mode(peakKVal)

% FullChar_2 = Probe.Block2000.pZav; figOff = 0;
FullChar_2 = Probe.Block50.pZav; figOff = 10;
% FullChar_2 = Probe.Block750.pZav; figOff = 20;
% FullChar_2 = analysis.approach.right.pAvg;

clear PopAnalysis

numClusts = 4;
[idx,C,sumd,D] = kmeans(FullChar_2,numClusts,'replicates',3);
[vals,inds] = sort(idx);
[mapName] = TNC_CreateRBColormap(1024,'mbr');
figure(22+figOff); subplot(121); imagesc(FullChar_2(inds,:),[-2 2]);
colormap(mapName);
figure(22+figOff); subplot(122); imagesc(corr(FullChar_2(inds,:)'),[-1 1]);

cMap = TNC_CreateRBColormap(numClusts,'bo');
figure(20+figOff); clf;
for k=1:numClusts
    
    tInds = find(idx==k);
    numel(tInds)
    if numel(tInds)>1
        PopAnalysis.clust(k).avg = mean(FullChar_2(tInds,:));
        PopAnalysis.clust(k).sem = std(FullChar_2(tInds,:),[],1) ./ sqrt(numel(tInds)-1);
        PopAnalysis.clust(k).cor = corr(FullChar_2(tInds,:));
    else 
        PopAnalysis.clust(k).avg = FullChar_2(tInds,:);
        PopAnalysis.clust(k).sem = zeros(1,numel(FullChar_2(tInds,:)));
        PopAnalysis.clust(k).cor = [FullChar_2(tInds,:) ; FullChar_2(tInds,:)];
    end
    PopAnalysis.clust(k).num = numel(tInds);
    
    PopAnalysis.timeOfPeak(k) = find(abs(PopAnalysis.clust(k).avg) == max(abs(PopAnalysis.clust(k).avg)),1,'first');
%     if max(PopAnalysis.clust(k).avg) - min(PopAnalysis.clust(k).avg) > 0
%         PopAnalysis.signOfPeak(k) = 1;
%     else
%         PopAnalysis.signOfPeak(k) = -1;
%     end
    
    PopAnalysis.signOfPeak(k) = mean(PopAnalysis.clust(k).avg(4000:4500));
    
end

% [valsP,indsP] = sort(PopAnalysis.timeOfPeak.*PopAnalysis.signOfPeak);
[valsP,indsP] = sort(PopAnalysis.signOfPeak);

for k=1:numClusts
    
    currInd = indsP(k);
    shadedErrorBar(-2000:5000 , PopAnalysis.clust(currInd).avg + (-k.*2), PopAnalysis.clust(currInd).sem , {'-k','color',cMap(k,:)}); hold on;
    plot([-2000 5000] , [(-k.*2) (-k.*2)] , 'k-' , 'Color' , [0.5 0.5 0.5]);
    text( 3000 , (-k.*2)+0.1 , num2str(PopAnalysis.clust(currInd).num) );
    
end

plot([0 0] , [0 (-k.*2)-2] , 'k-');
plot([2100 2100] , [0 (-k.*2)-2] , 'k--');
plot([2500 2500] , [0 (-k.*2)-2] , 'k--');
axis tight;

figure(21); clf;
for k=1:numClusts
    
    currInd = indsP(k);
    figure(21); subplot(3,3,k);
    imagesc(PopAnalysis.clust(currInd).cor,[-1 1]); colormap(mapName);
    hold on; plot([0 7000],[2000 2000],'k-'); plot([2000 2000],[0 7000],'k-');

end

%% PLOT the summary correlation with sigma for sorted population data

figure(500); clf;

shadedErrorBar(-2000:5000 , SummaryData.type(1).avg , SummaryData.type(1).sem , {'-k','color',cMap2(6,:)}); hold on;
shadedErrorBar(-2000:5000 , SummaryData.type(2).avg , SummaryData.type(2).sem , {'-k','color',cMap2(2,:)});
shadedErrorBar(-2000:5000 , SummaryData.type(3).avg , SummaryData.type(3).sem , {'-k','color',cMap2(10,:)});

plot([-2000 5000] , [0 0] , 'k--' , 'Color' , [0.5 0.5 0.5]);
plot([0 0] , [0.5 -1.5] , 'k-');
plot([2100 2100] , [0.5 -1.5] , 'k--' , 'Color' , [0.5 0.5 0.5]);
plot([2500 2500] , [0.5 -1.5] , 'k--' , 'Color' , [0.5 0.5 0.5]);

%% CURR SESSION >>> DISPLAY THE CHANGE IN APPROACH ALIGNED RESPONSES AROUND A TRANSITION

    figure(106);clf;    
    clear pre post
    
for i=2:14
    
    response = mean( [ mean(AllProbeTrials.PreSwitch.trial(i-1).pZav(:,2000:4500)) ; mean(AllProbeTrials.PreSwitch.trial(i).pZav(:,2000:4500)) ; mean(AllProbeTrials.PreSwitch.trial(i+1).pZav(:,2000:4500)) ] );
    p = polyfit([0:2500],response-response(1),1);
    
    pre.slope(i-1) = p(1);
    
    plot([0:2500] + ((i-1).*3500) , response - response(1) , 'Color' , 'k'); hold on;
    plot([0:2500] + ((i-1).*3500)  , polyval(p,[0:2500]) , 'k-' , 'LineWidth' , 2);
    plot([0 0] + ((i-1).*3500), [-1 1] , 'k--');
    
    
end

for i=2:12
    responseA = mean( [ mean(AllProbeTrials.PostSwitch.trial(i-1).pZav(:,2000:4500)) ; mean(AllProbeTrials.PostSwitch.trial(i).pZav(:,2000:4500)) ; mean(AllProbeTrials.PostSwitch.trial(i+1).pZav(:,2000:4500)) ] ) ;
    p = polyfit([0:2500],responseA-responseA(1),1);
    
    post.slope(i-1) = p(1);
    
    plot([0:2500] + ((i+13).*3500), responseA - responseA(1) , 'Color' , 'r');
    plot([0:2500] + ((i+13).*3500)  , polyval(p,[0:2500]) , 'r-' , 'LineWidth' , 2);
    plot([0 0] + ((i+13).*3500), [-1 1] , 'k--');

end

figure(107); clf; 
p1 = polyfit(-12:0,pre.slope,1);
p2 = polyfit(1:11,post.slope,1);

plot(-12:0,pre.slope,'ko'); hold on;
plot(-12:0,polyval(p1,-12:0),'k-','LineWidth',2);
plot(1:11,post.slope,'ro');
plot(1:11,polyval(p2,1:11),'r-','LineWidth',2);
axis([-12 12 -3.2e-4 0.2e-4]);



%% ALTERNATIVE APPROACH TO DEFINE RESPONSE TYPES IS TO USE PCA

numDims = 4;
[mappedB, mappingB] = compute_mapping(analysis.approach.right.pAvg,'PCA',numDims);

figure(21); clf
for i=1:numDims
    
    plot([-2000:5000] , mappingB.M(:,i) - (i.*0.05) ); hold on;
    
end

% figure(22);
% [mapName] = TNC_CreateRBColormap(2056,'bo');
% scatter(mappingB.M(:,1) , mappingB.M(:,2), 5 , [-2000:5000]); colormap(mapName);


%% >>> IF I JUST USED AN ABSOLUTE VALUE SORT OF COMPUTATION HOW WOULD I SCALE IT TO GROW WITH SIGMA BEST?

% dist(1).sigma = 50;
% dist(1).values = randn(1,1000) .* dist(1).sigma;
% 
% dist(2).sigma = 750;
% dist(2).values = randn(1,1000) .* dist(2).sigma;
% 
% dist(3).sigma = 2000;
% dist(3).values = randn(1,1000) .* dist(3).sigma;

sigmas = [50, 250, 500, 750, 1000, 1200, 1375, 1750, 2000]

for i=1:9
    
    dist(i).values = randn(1,20000) .* sigmas(i);
    sd(i) = std(dist(i).values);
    AbsAvg(i) = mean(abs(dist(i).values));
  
end

cMap = TNC_CreateRBColormap(21,'bo');
figure(100); clf;

for j=1:21
    alpha = ((j-1)./20);
    plot(sd,(1+alpha).*AbsAvg,'k.-','Color',cMap(j,:)); hold on;
    text(2050,(1+alpha).*AbsAvg(9),num2str(alpha));
end

plot(sigmas,sigmas,'ko-','LineWidth',1);
% plot(sigmas,AbsAvg,'k.','LineWidth',2);
axis([0 2500 0 3500]);

%% >>> MODEL of instantaneous probability that reward will come in the next delta
tot_time = 10000;
clear pRewInOff;
rew_pdf_750 = TNC_CreateGaussian(3000,750,tot_time,1);
rew_pdf_50 = TNC_CreateGaussian(3000,(3000.*0.16),tot_time,1);
rew_pdf_2000 = TNC_CreateGaussian(3000,2000,tot_time,1);
figure(1); subplot(4,1,1); plot(1:tot_time,rew_pdf_50,'r',1:tot_time,rew_pdf_750,'g',1:tot_time,rew_pdf_2000,'b');
offset = 250;

% rew_pdf_50=rew_pdf_50./max(rew_pdf_50);
% rew_pdf_750=rew_pdf_750./max(rew_pdf_750);
% rew_pdf_2000=rew_pdf_2000./max(rew_pdf_2000);

% plot(1:tot_time,(log(cumsum(rew_pdf_50))+log(rew_pdf_50)),'r',...
% 1:tot_time,(log(cumsum(rew_pdf_750))+log(rew_pdf_750)),'g',...
% 1:tot_time,(log(cumsum(rew_pdf_2000))+log(rew_pdf_2000)),'b');
% plot(1:tot_time,(1-(cumsum(rew_pdf_50)).*(rew_pdf_50)),'r',...
% 1:tot_time,(1-(cumsum(rew_pdf_750)).*(rew_pdf_750)),'g',...
for i=1:tot_time-offset
    pRewInOff_50(i) = trapz(rew_pdf_50(i:i+offset));
    pRewInOff_750(i) = trapz(rew_pdf_750(i:i+offset));
    pRewInOff_2000(i) = trapz(rew_pdf_2000(i:i+offset));
end
figure(1); subplot(4,1,2); 
plot(1:tot_time-offset,-log(pRewInOff_50),'r',1:tot_time-offset,-log(pRewInOff_750),'g',1:tot_time-offset,-log(pRewInOff_2000),'b'); axis([0 5000 0 12]);
figure(1); subplot(4,1,3); 
% plot(1:tot_time-offset,log(cumsum(pRewInOff_50)),'r',1:tot_time-offset,log(cumsum(pRewInOff_750)),'g',1:tot_time-offset,log(cumsum(pRewInOff_2000)),'b'); axis([0 5000 -6 6]);
plot(1:tot_time-offset,log(0.85.*(1:tot_time-offset))-log(cumsum(pRewInOff_50)),'r',1:tot_time-offset,log(0.85.*(1:tot_time-offset))-log(cumsum(pRewInOff_750)),'g',1:tot_time-offset,log(0.85.*(1:tot_time-offset))-log(cumsum(pRewInOff_2000)),'b'); axis([0 5000 0 12]);
% axis([0 6000 -20 -5]);
figure(1); subplot(4,1,4); 
figure(3);
plot(1:tot_time-offset,-log(pRewInOff_50)-(log(0.85.*(1:tot_time-offset))-log(cumsum(pRewInOff_50))),'r',1:tot_time-offset,-log(pRewInOff_750)-(log(0.85.*(1:tot_time-offset))-log(cumsum(pRewInOff_750))),'g',1:tot_time-offset,-log(pRewInOff_2000)-(log(0.85.*(1:tot_time-offset))-log(cumsum(pRewInOff_2000))),'b');%,1:tot_time-offset,log(0.85.*(1:tot_time-offset))-log(cumsum(pRewInOff_2000)),'b'); 
axis([0 5000 -2 2]);


%% >>> HOW IS Treward COMPUTED FOR EACH TRIAL
% One idea is that the dissimilarity of the activity trace evoked by
% left presses and right presses grows linearly with time. Should be
% revealed by alignment on approach.

%% >>> HOW IS Twait CHOSEN/READOUT FROM THE MODEL