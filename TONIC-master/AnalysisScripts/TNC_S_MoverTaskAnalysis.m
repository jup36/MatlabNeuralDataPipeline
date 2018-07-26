% script for analysis of mover, sigmamover, and mitomover

%% ANALYSIS PARAMETERS
% Define the smoothing to apply to the timeseries data
currParams.smthParams.rise     = 1;
% currParams.smthParams.decay    = 100;
% currParams.smthParams.decay    = 5;
% currParams.smthParams.decay    = 10;
currParams.smthParams.decay    = 25;
% currParams.smthParams.decay    = 60;
currParams.filter.causal = 0;

if currParams.filter.causal
    [currParams.filter.kernel]  = TNC_CreateCausalKernel(currParams.smthParams.rise,currParams.smthParams.decay,1);
else
    [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
end

currParams.winParams.prior     = 2e3;
currParams.winParams.after     = 3e3; 

currParams.popVec.winSize       = 1;

currParams.stableTrials         = 5; %derived from summary behavior data

PopData.currParams = currParams;

halfVal = max(currParams.filter.kernel) ./ 2;

% figure(1); clf; plot(PopData.currParams.filter.kernel,'ko-','LineWidth',2,'Color',[0 0 0]);

disp(' ');
disp(' ');
disp('________________________________');
disp(' ');
disp('Initialized analysis parameters.');
disp(['Kernel FWHM: ' num2str( find(currParams.filter.kernel>halfVal,1,'last') - find(currParams.filter.kernel>halfVal,1) ) ' ms']);
disp('________________________________');
disp(' ');
disp(' ');

% need 5 colors for displaying sorted units
colorPalette(1,:)   = [1 0 0];
colorPalette(2,:)   = [0 0.67 1];
colorPalette(3,:)   = [0 0 0];
colorPalette(4,:)   = [1 0.67 0];
colorPalette(5,:)   = [0.75 0.75 0.75];

%% ALIGN RASTERS FOR LEVER MVMTS

disp(['___________________________________________________'])
disp(['STARTED aligning all raster plots and updating the PopData structure...']);

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
psthWin = [0.5e3,1.5e3];
k=1;

[mapName] = TNC_CreateRBColormap(500,'rb');


for i=1:NumSessions

    numUnits = size(PopData.session(i).unit,2);
    plotDim = ceil(sqrt(numUnits));
    countUnits = countUnits+numUnits;
    PopData.session(i).reaches = size(ContData.behavior.reach.start,1);

    for j = 1:numUnits

        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,ceil(PopData.session(i).unit(j).ts)) = 1;

        [respCS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,ContData.behavior.reach.start(:,4)+ContData.behavior.reach.vel(:,4),psthWin,1,1);
%         [respCS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,ContData.behavior.reach.stop(:,4),psthWin,1,1);
%         [respCS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,ContData.behavior.reach.start(:,4),psthWin,1,1);
        PopData.session(i).unit(j).reach.raster      = respCS.raster;

        tmpSmooth = conv(delta,currParams.filter.kernel,'same');
        [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,ContData.behavior.reach.start(:,4)+ContData.behavior.reach.vel(:,4),psthWin,0,1);
%         [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,ContData.behavior.reach.stop(:,4),psthWin,0,1);
%         [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,ContData.behavior.reach.start(:,4),psthWin,0,1);
        PopData.session(i).unit(j).reachS            = respCSS;  

        % Display mean PSTHs
        figure(3); 
        if j==1 & i==1
            clf;
        end
        subplot(plotDim-1,plotDim,j);
%         patch((PopData.session(i).unit(j).respCSS.image.psthZ+PopData.session(i).unit(j).respCSS.image.psthZe)',(PopData.session(i).unit(j).respCSS.image.psthZ-PopData.session(i).unit(j).respCSS.image.psthZe)',colorPalette(2,:));
%         hold on;
        plot(-psthWin(1):psthWin(2),PopData.session(i).unit(j).reachS.image.psthZ,'Color',colorPalette(2,:)./2,'LineWidth',2);
                hold on;
        plot(-psthWin(1):psthWin(2),PopData.session(i).unit(j).reachS.image.psthZ-PopData.session(i).unit(j).reachS.image.psthZe,'-','Color',colorPalette(2,:),'LineWidth',1);
        plot(-psthWin(1):psthWin(2),PopData.session(i).unit(j).reachS.image.psthZ+PopData.session(i).unit(j).reachS.image.psthZe,'-','Color',colorPalette(2,:),'LineWidth',1);
        plot([-psthWin(1),psthWin(2)],[0 0],'--','Color',[0.2 0.2 0.2]);
        plot([0 0],[-0.5 2],'--','Color',[0.2 0.2 0.2]);
        axis([-psthWin(1) psthWin(2) -0.5 2]);
        title([num2str(1000./PopData.session(i).unit(j).isi.stats.mean) ' Hz ']);
        k=k+1;

        PopPSTH(j,:) = PopData.session(i).unit(j).reachS.image.psthZ ./ max(PopData.session(i).unit(j).reachS.image.psthZ);
        
    end

    figure(3); 
    subplot(plotDim-1,plotDim,j+1);
    imagesc(corr(PopPSTH(:,window)'),[-1 1]);
    colormap(mapName);
    
    disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' of ' num2str(NumSessions)])

end

disp(['___________________________________________________'])
disp(['COMPLETED ...']);

%% COMPUTE DIMENSIONS THAT CAPTURE VARIANCE OF POPULATION AROUND REACHES

window = 1:1000;
widthSig = currParams.smthParams.decay;
thisSession = 83;
sortedByLat = 1;


i = 1;
thisSession = i;

    [PopVec.kernel]  = TNC_CreateGaussian(widthSig.*15,widthSig,widthSig.*30,1);
        
%     PhaseMat = PopData.session(thisSession).PhaseMat;
%     numUnits = size(PhaseMat,1);

disp(' ');
disp(' ');
disp('________________________________');
disp('Computing principal components...');
disp(' ');

    % First we consider the basic features of temporal evolution of activity (an analogy to the continuous behavioral data we intend to extract)
    covMatForTimeCourse = cov(PopPSTH(:,window)');
    [v,d] = eig(covMatForTimeCourse);
    TimeCourse.pca.vecs = zeros(size(PopPSTH,1),3);
    TimeCourse.pca.vals = diag(d)./sum(diag(d));
    TimeCourse.pca.vecs(:,1) = v(:,numel(TimeCourse.pca.vals));
    TimeCourse.pca.vecs(:,2) = v(:,numel(TimeCourse.pca.vals)-1);
    TimeCourse.pca.vecs(:,3) = v(:,numel(TimeCourse.pca.vals)-2);
    TimeCourse.pca.varExp = sum(TimeCourse.pca.vals(numel(TimeCourse.pca.vals)-2:numel(TimeCourse.pca.vals)));
    disp(['Variance explained (time): ' num2str(sum(TimeCourse.pca.vals(numel(TimeCourse.pca.vals)-2:numel(TimeCourse.pca.vals)))) ' | ' num2str(TimeCourse.pca.vals(numel(TimeCourse.pca.vals))) ' ' num2str(TimeCourse.pca.vals(numel(TimeCourse.pca.vals)-1)) ' ' num2str(TimeCourse.pca.vals(numel(TimeCourse.pca.vals)-2)) ]);
    disp(' ');
 
%     % Second we consider the evolution of the population vector using the warped response patterns
%     covMatForPopVec = cov(PhaseMat(:,window)');
%     [v,d] = eig(covMatForPopVec);
%     PopVec.pca.vecs = zeros(size(PhaseMat,1),3);
%     PopVec.pca.vals = diag(d)./sum(diag(d));
%     PopVec.pca.vecs(:,1) = v(:,numel(PopVec.pca.vals));
%     PopVec.pca.vecs(:,2) = v(:,numel(PopVec.pca.vals)-1);
%     PopVec.pca.vecs(:,3) = v(:,numel(PopVec.pca.vals)-2);
%     PopVec.pca.varExp = sum(PopVec.pca.vals(numel(PopVec.pca.vals)-2:numel(PopVec.pca.vals)));
%     disp(['Variance explained (population): ' num2str(sum(PopVec.pca.vals(numel(PopVec.pca.vals)-2:numel(PopVec.pca.vals)))) ' | ' num2str(PopVec.pca.vals(numel(PopVec.pca.vals))) ' ' num2str(PopVec.pca.vals(numel(PopVec.pca.vals)-1)) ' ' num2str(PopVec.pca.vals(numel(PopVec.pca.vals)-2))]);

disp(' ');
disp('Visualizing structure of covariance matrices...');
disp(' ');

    figure(21);
    subplot(2,2,1);
    imagesc(covMatForTimeCourse);
    title('Covariance of PSTH time courses');

    subplot(2,2,2); hold off;
    plot(TimeCourse.pca.vecs(:,1),'Color',[0 0 0],'LineWidth',2); hold on;
    plot(TimeCourse.pca.vecs(:,2),'Color',[1 0 0],'LineWidth',2);
    plot(TimeCourse.pca.vecs(:,3),'Color',[0 0.67 1],'LineWidth',2);
    title('Leading eigenvectors for TimeCourse');
% 
%     subplot(2,3,3); hold off;
%     plot3(v(1,:),v(2,:),v(3,:),'LineWidth',2);
%     grid on; 
%    
%     title('Evolution of TimeCourse dimensions');

%     subplot(2,3,4);
%     imagesc(covMatForPopVec);
%     title('Covariance of individual PSTHs');
% 
%     subplot(2,3,5); hold off;
%     plot(PopVec.pca.vecs(:,1),'Color',[0 0 0],'LineWidth',2); hold on;
%     plot(PopVec.pca.vecs(:,2),'Color',[1 0 0],'LineWidth',2);
%     plot(PopVec.pca.vecs(:,3),'Color',[0 0.67 1],'LineWidth',2);
%     title('Leading eigenvectors for PopVec');

%     subplot(2,3,6); hold off;
%     % project all units for display
%     for i=1:size(PhaseMat,1)
%         TimeCourse.proj(1,i) = dot(PhaseMat(i,window),TimeCourse.pca.vecs(:,1));
%         TimeCourse.proj(2,i) = dot(PhaseMat(i,window),TimeCourse.pca.vecs(:,2));
%         TimeCourse.proj(3,i) = dot(PhaseMat(i,window),TimeCourse.pca.vecs(:,3));
% 
%         subplot(2,3,6);
% 
%         plot3(TimeCourse.proj(1,i),TimeCourse.proj(2,i),TimeCourse.proj(3,i),'o','Color',[1-(i./numUnits) (0.67.*(i./numUnits)) i./numUnits],'LineWidth',2,'MarkerSize',8); hold on;
%     end
%     title('Projection onto TimeCourse dimensions');
%     grid on;

disp(' ');
disp('Projecting individual trial trajectories in low dimensions...');
disp(' ');
    
    numTrials = size(ContData.behavior.reach.start,1);
    duration = size(PopData.session(thisSession).unit(1).reachS.image.aligned,2); % dims are [trials x time]

%     if sortedByLat==1
%         disp('Sorting by latency...');
%         yLabelStr = '(sorted)';
%         latencies = PopData.session(thisSession).behavior.flCSon(validTrials);
%         [vals,trialIndsTmp] = sort(latencies);
%         trialInds = validTrials(trialIndsTmp);
%     else
%         disp('No sorting was applied...');
%         trialInds = validTrials;
%         yLabelStr = '(not sorted)';
%     end

    trialIndsTmp = 1:numTrials;

    for k=1:numTrials
        
        numUnits = size(PopData.session(thisSession).unit,2);
        
        VectorMat = zeros(numUnits,duration); 
        TimeCourse.trial(k).proj = zeros(3,duration);
        
        
        for l=1:numUnits
                        
            trialMat(l,:) = PopData.session(thisSession).unit(l).reachS.image.aligned(k,:);
            
        end
        
        for m = 1:duration
            PopVec.trial(k).proj(1,m) = dot(trialMat(:,m),TimeCourse.pca.vecs(:,1));
            PopVec.trial(k).proj(2,m) = dot(trialMat(:,m),TimeCourse.pca.vecs(:,2));
            PopVec.trial(k).proj(3,m) = dot(trialMat(:,m),TimeCourse.pca.vecs(:,3));
        end
        
    end
    
% disp(' ');
% disp('Visualizing individual trial trajectories...');
% disp(' ');
% 
%     figure(22); clf; hold on;
%     bins = 5;
%     binUnits = floor(numTrials./bins)
%     PopVec.trial(k).proj = zeros(3,7300); 
%     
%     for k=1:bins
% 
%         tempMat1 = zeros(binUnits,7300);
%         tempMat2 = zeros(binUnits,7300);
%         tempMat3 = zeros(binUnits,7300);
% 
%         for m=1:binUnits
%             thisTrial = m + ((k-1)*binUnits);
%             tempMat1(m,:) = PopVec.trial(thisTrial).proj(1,:);
%             tempMat2(m,:) = PopVec.trial(thisTrial).proj(2,:);
%             tempMat3(m,:) = PopVec.trial(thisTrial).proj(3,:);
%         end
%         
%         PopVec.trialBin(k).proj         = zeros(3,7300);
%         PopVec.trialBin(k).proj(1,:)    = mean(tempMat1);
%         PopVec.trialBin(k).proj(2,:)    = mean(tempMat2);
%         PopVec.trialBin(k).proj(3,:)    = mean(tempMat3);
%         
%         % colorize by fl latency
%         figure(22);        
%         thisColor = [1-(k./bins) (0.67.*(k./bins)) k./bins];
%         window = 1:7000;
%         plot3(PopVec.trialBin(k).proj(1,window),PopVec.trialBin(k).proj(2,window),PopVec.trialBin(k).proj(3,window),'-','Color',thisColor,'LineWidth',4); 
%         plot3(PopVec.trialBin(k).proj(1,1000),PopVec.trialBin(k).proj(2,1000),PopVec.trialBin(k).proj(3,1000),'o','Color',thisColor,'LineWidth',2,'MarkerSize',10);
%         plot3(PopVec.trialBin(k).proj(1,1718),PopVec.trialBin(k).proj(2,1718),PopVec.trialBin(k).proj(3,1718),'s','Color',thisColor,'LineWidth',2,'MarkerSize',10);
%         plot3(PopVec.trialBin(k).proj(1,2795),PopVec.trialBin(k).proj(2,2795),PopVec.trialBin(k).proj(3,2795),'d','Color',thisColor,'LineWidth',2,'MarkerSize',10);
%         plot3(PopVec.trialBin(k).proj(1,1055:1150),PopVec.trialBin(k).proj(2,1055:1150),PopVec.trialBin(k).proj(3,1055:1150),'k--','LineWidth',4);
%         grid on; view([144 62]);
%         xlabel('PCA 1');
%         ylabel('PCA 2');
%         zlabel('PCA 3');
%     end
%             
%     PopData.session(thisSession).PopVec = PopVec;
%     PopData.session(thisSession).TimeCourse = TimeCourse;
% 
%     clear PopVec TimeCourse;
% 
% end    
%     
%     
%     
%     
% disp(' ');
% disp('________________________________');
% disp(' ');

%% VIUALIZING THE LOW DIMENSIONAL PROJECTION OF THE DATA

startT = 250;
endT = 750;

disp('Sorting by max velocity...');
[vals,sortInds] = sort(ContData.behavior.reach.dur(:,1),'descend');
figure(5); clf; plot(ContData.behavior.reach.dur(sortInds,1)); pause(1); drawnow;

figure(5); clf;

for p = [1:50 201:250]
    
    i = sortInds(p);
    
    if p<100
        colorSpec = [1 0 0];
    else
        colorSpec = [0 0.67 1];
    end        
    
    figure(5); plot3(PopVec.trial(i).proj(1,500),PopVec.trial(i).proj(2,500),PopVec.trial(i).proj(3,500),'o','MarkerSize',5,'Color',colorSpec); grid on;
    hold on;
    figure(5); plot3(PopVec.trial(i).proj(1,endT),PopVec.trial(i).proj(2,endT),PopVec.trial(i).proj(3,endT),'s','MarkerSize',20,'Color',colorSpec); grid on;
%     figure(5); plot3(PopVec.trial(i).proj(1,startT:500),PopVec.trial(i).proj(2,startT:500),PopVec.trial(i).proj(3,startT:500),'Color',[0.7 0.7 0.7]); grid on;
%     figure(5); plot3(PopVec.trial(i).proj(1,500:endT),PopVec.trial(i).proj(2,500:endT),PopVec.trial(i).proj(3,500:endT),'Color',colorSpec); grid on;

%     figure(5); plot(PopVec.trial(i).proj(3,500),PopVec.trial(i).proj(2,500),'o','MarkerSize',4,'Color',colorSpec); grid on;
%     hold on;
%     figure(5); plot(PopVec.trial(i).proj(3,endT),PopVec.trial(i).proj(2,endT),'s','MarkerSize',20,'Color',colorSpec); grid on;
%     figure(5); plot(PopVec.trial(i).proj(3,startT:500),PopVec.trial(i).proj(2,startT:500),'Color',[0.7 0.7 0.7]); grid on;
%     figure(5); plot(PopVec.trial(i).proj(3,500:endT),PopVec.trial(i).proj(2,500:endT),'Color',colorSpec); grid on;

end

%% COMPILE ALL PSTHS

%% FIND PEAKS

%% TRIAL BY TRIAL COMPARISONS - For reaches display the psth of all units

%% DEPRECATED CODE
% %% LOAD NEX DATA
% % WALK THROUGH ALL FILES IN A DIRECTORY AND LOAD DATA INTO MEMORY
% 
% clear
% 
% preCheck        = 0;
% error           = 0;
% stableTrials    = 10;
% winPsth         = 10;
% 
% RunTonicInfo.activeSet = 1;
% RunTonicInfo.loadedFiles = 0;
% totalFiles = 0;
% extLogic =0;
% 
% % create a list of all the folders in the directory
% allDirectories  = dir; % creates a structure, use the "name" field
% dirs            = size(allDirectories,1);
% currSes         = 0;
% unitCount       = 0;
% 
% disp(' ');
% disp(' ');
% disp(' ');
% disp(' ');
% disp(['________________LOADING ALL NEX FILES IN DIRECTORY_________________________']);
% disp(['Directory: ' pwd]);
% disp(['Timestamp: ' datestr(now)]);
% 
% fid = fopen([date '_' num2str(round(rand(1).*100)) '.log'],'w');
% 
% fprintf(fid,['\n________________Running WalkAllDir... Script_________________________\n']);
% fprintf(fid,['Directory: ' pwd '\n']);
% fprintf(fid,['Timestamp: ' datestr(now) '\n']);
% 
% disp(' ');
% 
% % loop through all the directories
% for i = 1:dirs
%     
%     % test for a valid directory
%     if allDirectories(i).isdir == 1
%         
%         dirPath = allDirectories(i).name;
%         if strncmp('.',dirPath,1)==0
%             disp(['________________' allDirectories(i).name '_________________________']);
%             fprintf(fid,['\n________________' allDirectories(i).name '_________________________\n']);
%             
%             % create a list of all the files
%             allFiles = dir(sprintf('%s/*.nex',dirPath));
%                         
%             for j = 1:size(allFiles,1)
%                 disp(' ');
%                 disp(allFiles(j).name);   
% 
%                 currSes = currSes+1;
%                 disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
%                 disp(['Current file index: ' num2str(currSes)]);
%                 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %_____________APPLY ANALYSIS TO INDIVIDUAL FILES HERE_____________         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                     clear propUnits; clear propEvents; p = 1; propUnits.unitInds = []; 
% 
%                     RunTonicInfo.activeSet = RunTonicInfo.loadedFiles;
%                     RunTonicInfo.loadedFiles = RunTonicInfo.loadedFiles+totalFiles;
% 
%                     fileNameStr = [dirPath '/' allFiles(j).name];
%                     tmpName = allFiles(j).name(1:strfind(allFiles(j).name,'.')-1);
% 
%                     % remove '-'s from the name
%                     nonDashInds = strfind(tmpName,'-');
%                     tmpName(nonDashInds) = '_';
%                     % remove ' ' from the name
%                     spaceInds = strfind(tmpName,' ');
%                     tmpName(spaceInds) = '_';
%                     nameInMem = tmpName;
% 
%                     RunTonicInfo.activeSet = RunTonicInfo.activeSet+1;
%                     RunTonicInfo.fileNames(RunTonicInfo.activeSet).name     = nameInMem;
%                     RunTonicInfo.fileNames(RunTonicInfo.activeSet).loadTime = datestr(now);
% 
%                     RD_Active = TNC_LoadData(0, 0, fileNameStr);
% 
%                     for m = 1:size(RD_Active.neurons,1)    
%                         if findstr(RD_Active.neurons{m}.name,'elec') == 1% electrode
% 
%                             testi = findstr(RD_Active.neurons{m}.name,'U');
%                             if size(testi,1) == 0
%                                 
%                                 propUnits.unitInds = [propUnits.unitInds,m];
%                                 propUnits.names(p).str = RD_Active.neurons{m}.name;
%                                 p = p+1;
%                                 
%                                 fprintf(fid,[RD_Active.neurons{m}.name '\n']);                                
% 
%                             end
%                                                        
%                         end
%                         
%                         if findstr(RD_Active.neurons{m}.name,'ainp9a')
%                             disp(RD_Active.neurons{m}.name);
%                         end
%                     end
%                     
%                     % Digital signals
%                         % early data files are -0006, -0007, and -0008
%                         
%                     
%                     % create a unique id for the session
%                     PopData.session(currSes).sessId = [nameInMem '>>' num2str(currSes)];
% 
%                     disp(       ['Session id: ' PopData.session(currSes).sessId ' | units: ' num2str(size(propUnits.unitInds,2))]);
%                     fprintf(fid,['Session id: ' PopData.session(currSes).sessId ' | units: ' num2str(size(propUnits.unitInds,2)) '\n']);
%                     
%                     % FLAG TO PREVENT COMPLETE DATA EXTRACTION
% 
%                     if preCheck == 0
%                         
%                         % for reference: [tonicDataStructure] = TNC_PhotoStimToStruct(dataStructure,unitArrayToLoad,toneID,liteID)
%                         [TD_Active] = TNC_MoverToStruct(RD_Active,propUnits.unitInds);
% 
%                         % Store the general session event data:
%                         PopData.session(currSes).events         = TD_Active.events;
%                         
%                         kMax = size(TD_Active.unit,2);
%                         disp(['Total units: ' num2str(kMax)]);
%                         
%                         for k=1:kMax
% 
%                             % create a unique id for the cell/session pair
%                             PopData.session(currSes).unit(k).uID  = [num2str(currSes) '_' num2str(k)];    
%                             PopData.session(currSes).unit(k).name = TD_Active.unit(k).name;    
%                             PopData.session(currSes).unit(k).ts   = TD_Active.unit(k).ts;   
% 
%                             % extract and store waveforms
%                             PopData.session(currSes).unit(k).WfMean = mean(TD_Active.unit(k).wf,2);
%                             PopData.session(currSes).unit(k).WfStd  = std(TD_Active.unit(k).wf,0,2);
% 
%                             % calculate and store isis
%                             [isi] = TNC_QuantISI(TD_Active.unit(k).ts);
%                             PopData.session(currSes).unit(k).isi = isi;    
% 
%                         end
% 
%                     end
% 
%                     disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
%                     disp('  ');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %_____________APPLY ANALYSIS TO INDIVIDUAL FILES HERE_____________         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%             end
% 
%             disp('  ');
%         end
%         
%     end
%     
% end
% 
% disp(['Unit count: ' num2str(unitCount)]);
% clear DA_Ma*
% 
% fprintf(fid,'_____________________________________________________________________________________\n');
% fprintf(fid,['Unit count: ' num2str(unitCount) '\n']);
% % close the file
% fclose(fid);
% 
% disp(['Saving data as ' nameInMem(1:8) 'PD.mat']);
% eval(['save ' nameInMem(1:8) 'PD PopData']);
% 
% %% EXTRACT NEX DATA AND CONVERT TO SORTING FORMAT
% 
% %% LOAD NEV DATA (or can be done manually into a structure named NEV)
% [filename,pathname] = uigetfile('.nev');
% fileNameStr = [pathname filename]; 
% 
% disp(' ');disp(' ');
% disp(['Loading the file>>> ' filename]);
% disp(['...from the path>>> ' pathname]);
% disp(' ');disp(' ');
% 
% [NEV] = TNC_LoadData(0, 0, fileNameStr);
% 
% %% CONVERT NEV DATA TO SORTING FORMAT
% 
% % Pack into the standard TONIC formatted data structure
% 
% 
% % for j = 1:128 % only check ephys channels at the moment
% %     
% %     evOnCurrChan = find(NEV.Data.Spikes.Electrode==j);
% % 
% %     if numel(evOnCurrChan)>0
% %         PreSortData.session(i).electrode(j).ons.events.ts   = double(NEV.Data.Spikes.Timestamps(evOnCurrChan));
% %         PreSortData.session(i).electrode(j).ons.events.wfs  = double(NEV.Data.Spikes.Waveform(evOnCurrChan,:));
% %         PreSortData.session(i).electrode(j).ons.events.uid  = double(NEV.Data.Spikes.Unit(evOnCurrChan));
% %         PreSortData.session(i).electrode(j).ons.events.tsRes= NEV.MetaTags.SampleRes;
% %         PreSortData.session(i).electrode(j).ons.events.thr  = NEV.ElectrodesInfo{j}.LowThreshold;
% %         PreSortData.session(i).electrode(j).ons.events.winL = 10;
% %         PreSortData.session(i).electrode(j).ons.events.winR = size(NEV.Data.Spikes.Waveform,2)-10+1;    
% %     else
% %         PreSortData.session(i).electrode(j).ons.events.ts    = [];
% %     end
% %     
% % end
% % 
% % for j = 129:144 % now get event data
% %     
% %     evOnCurrChan = find(NEV.Data.Spikes.Electrode==j);
% % 
% %     if numel(evOnCurrChan)>0
% %         PreSortData.session(i).events.analog(j-128).ts      = double(NEV.Data.Spikes.Timestamps(evOnCurrChan));
% %         PreSortData.session(i).events.analog(j-128).uid     = double(NEV.Data.Spikes.Unit(evOnCurrChan));
% %     else
% %         PreSortData.session(i).events.analog(j-128).ts      = [];
% %         PreSortData.session(i).events.analog(j-128).uid     = [];
% %     end
% % end
% % 
% % % Align the events
% % for j = 1:128 % only check ephys channels at the moment
% %     if numel(PreSortData.session(i).electrode(j).ons.events.ts)>0
% % 
% %         if align==1
% %             [PreSortData.session(i).electrode(j).ons.events] = TNC_EventAlign(PreSortData.session(i).electrode(j).ons.events,clean,upRatio);
% %         else
% %             disp(['Upsampling channel ' num2str(j) ' ...']);
% %             indSize = size(PreSortData.session(i).electrode(j).ons.events.wfs);  
% %             upRatio = 10;
% %             x   = 1:1:indSize(1,2);
% %             xx  = 1:(1/upRatio):indSize(1,2);
% %             PreSortData.session(i).electrode(j).ons.events.Sh_x = xx;
% %             PreSortData.session(i).electrode(j).ons.events.upRatio = upRatio;
% %             
% %             PreSortData.session(i).electrode(j).ons.events.Sh_wfs = zeros(indSize(1,1),numel(xx));
% %             for k = 1:indSize(1,1)
% %                 % use spline interpolation to get smoother spikes
% %                 splineSpk(1,:) = spline(x,PreSortData.session(i).electrode(j).ons.events.wfs(k,:),xx);
% %                 PreSortData.session(i).electrode(j).ons.events.Sh_wfs(k,:) = splineSpk;
% %             end
% %         end
% %     end
% % end
% 
% %% POPUP A DISPLAY OF THE CURRENT CHANNELS FOR QUICK EVALUATION
% i=1;
% plotMean = 1;
% 
% % plot in a 4x8 grid
% 
% figure(1); clf;
% 
% for j=1:32
%     
%     disp(['Channel ' num2str(j) '  |  ' num2str(numel(PreSortData.session(i).electrode(j).ons.events.ts)) ' events.']);
%     
%     if numel(PreSortData.session(i).electrode(j).ons.events.ts)>0
%         
%         clear meanSpike errSpike
%         
%         numSamps = size(PreSortData.session(i).electrode(j).ons.events.Sh_wfs,2);
%         
%         disp(['Each event has ' num2str(numSamps) ' points.']);
%         uUIDsRaw = unique(PreSortData.session(i).electrode(j).ons.events.uid);
%         
%         tmp = find(uUIDsRaw>0 & uUIDsRaw<10);
%         uUIDs = uUIDsRaw(tmp);
%         
%         figure(1);
% 
%         if plotMean
%             
%             subplot(4,8,j); hold on;
% 
%             for p = 1:numel(uUIDs)
%                 currUID = uUIDs(p);
%                 clear allspikes
%                 allWfs = find(PreSortData.session(i).electrode(j).ons.events.uid==currUID);
%                 allspikes = PreSortData.session(i).electrode(j).ons.events.Sh_wfs(allWfs,:);
%                 meanSpike(p,:)  = mean(allspikes,1);
%                 errSpike(p,:)   = std(allspikes,[],1);
%                 
%                 currColor = colorPalette(currUID,:);
%                 plot(1:numSamps, meanSpike(p,:), '-', 'Color', currColor, 'LineWidth', 2);
%                 plot(1:numSamps, meanSpike(p,:)+errSpike(p,:), '--', 'Color', currColor, 'LineWidth', 1);
%                 plot(1:numSamps, meanSpike(p,:)-errSpike(p,:), '--', 'Color', currColor, 'LineWidth', 1);
% 
%                 if p==numel(uUIDs)
%                     PreSortData.session(i).electrode(j).ons.meanSpike = meanSpike;
%                     PreSortData.session(i).electrode(j).ons.errSpike = errSpike;
%                 end
% 
%             end
% 
% 
%         else
%             for k = 1:numel(PreSortData.session(i).electrode(j).ons.events.ts)
% 
%                 currWf = PreSortData.session(i).electrode(j).ons.events.wfs(k,:);
% 
%                 if k==1
%                     subplot(4,8,i); hold on;
%                 end
%                 if PreSortData.session(i).electrode(j).ons.events.uid(k)==0
%                     currColor = [0.2 0.2 0.2];
%                 elseif PreSortData.session(i).electrode(j).ons.events.uid(k)>5
%                     currColor = [0.2 0.2 0.2];
%                     plot(1:numSamps, currWf, 'Color', currColor);                
%                 else
%                     currColor = colorPalette(PreSortData.session(i).electrode(j).ons.events.uid(k),:);
%                     plot(1:numSamps, currWf, 'Color', currColor);
%                 end
%             end
%         end
%         
%         title(['Channel ' num2str(j)]);
% %         pause();
%         
%     end
%     
% end
% 
% %% CHECK the ONLINE SORTING AND ADD IN UNSORTED SPIKES
% % TO DO: Compile all of this data into a single report about the quality of the sorting
% 
% i=1;
% plotMean = 1; plotTemplates=0;
% 
% ChanSelector = 27;
% figure(2); clf;
% 
% j = ChanSelector;
% disp(['Channel ' num2str(j) '  |  ' num2str(numel(PreSortData.session(i).electrode(j).ons.events.ts)) ' events.']);
%     
%     if numel(PreSortData.session(i).electrode(j).ons.events.ts)>0
%         eventsLocal = PreSortData.session(i).electrode(j).ons.events;
% 
%         upRatio = PreSortData.session(i).electrode(j).ons.events.upRatio;
%         
% %         if isfield(eventsLocal,'scalar')==0
%             disp('Calculating scalar values for spikes.');
%             [eventsLocal] = TNC_EventQuant(eventsLocal,'scalar','dot',upRatio);
% %         end
% %         if isfield(eventsLocal,'gmono')==0
%             disp('Calculating gmono values for spikes.');
%             [eventsLocal] = TNC_EventQuant(eventsLocal,'gmono','dot',upRatio);
% %         end
% %         if isfield(eventsLocal,'frequency')==0
%             disp('Calculating frequency values for spikes.');
%             [eventsLocal] = TNC_EventQuant(eventsLocal,'frequency','dot',upRatio);
% %         end
% %         if isfield(eventsLocal,'pca')==0
%             disp('Calculating pca values for spikes.');
%             [eventsLocal] = TNC_EventQuant(eventsLocal,'pca','dot',upRatio);
% %         end
% 
%         uUIDsRaw = unique(PreSortData.session(i).electrode(j).ons.events.uid);
% 
%         tmp = find(uUIDsRaw>0 & uUIDsRaw<10);
%         uUIDs = uUIDsRaw(tmp);
% 
%         figure(2); clf;
%         subplot(442); hold on;
% %         plot3(eventsLocal.pca.corrSh.values(:,3),eventsLocal.pca.corrSh.values(:,2),eventsLocal.gmono.dot.values(:,2),'o','MarkerSize',6,'Color',[0.5 0.5 0.5]);
%         plot3(eventsLocal.gmono.dot.values(:,1),eventsLocal.scalar.scl.values(:,6),eventsLocal.scalar.scl.values(:,7),'o','MarkerSize',3,'Color',[0.5 0.5 0.5]);
%         view([65 40]); grid on;
% 
%         for q = 1:numel(uUIDs)
%             currUID = uUIDs(q);
%             indices = find(PreSortData.session(i).electrode(j).ons.events.uid == currUID);
%             currColor = colorPalette(currUID,:);
%             
%             figure(2);
%             subplot(442); hold on;
% %             plot3(eventsLocal.pca.corrSh.values(indices,3),eventsLocal.pca.corrSh.values(indices,2),eventsLocal.gmono.dot.values(indices,2),'.','MarkerSize',6,'Color',currColor);            
%             plot3(eventsLocal.gmono.dot.values(indices,1),eventsLocal.scalar.scl.values(indices,6),eventsLocal.scalar.scl.values(indices,7),'.','MarkerSize',2,'Color',currColor);            
%             xlabel('Gmono'); ylabel('Neg Width'); zlabel('Positive Width')
%             
%             % check stability over time
%             figure(2); subplot(4,4,[14:16]);view([0 90]); hold on;
%             plot(eventsLocal.ts(indices),eventsLocal.scalar.scl.values(indices,4),'.','Color',currColor);
%             
%             % check out the correlation matrix for all values
%             forCorr  = [eventsLocal.scalar.scl.values(indices,[4,6:7])';eventsLocal.gmono.dot.values(indices,:)';eventsLocal.pca.dot.values(indices,:)';eventsLocal.frequency.dot.values(indices,:)'];
%             forCorr2 = [eventsLocal.scalar.scl.values(indices,3:5)'];
% %             forCorr = [eventsLocal.gmono.corrSh.values(indices,:)'];
% %             forCorr = [eventsLocal.scalar.scl.values(indices,1:3)'];
% %             forCorr2 = [eventsLocal.pca.corrSh.values(indices,:)'];
% %             forCorr2 = [eventsLocal.Sh_wfs(indices,:)'];
% 
%             if q==1
%                 finalCorr  = forCorr;
%                 finalCorr2 = forCorr2;
%             else
%                 finalCorr  = [finalCorr,forCorr];                
%                 finalCorr2 = [finalCorr2,forCorr2];                
%             end
%             
%             meanSpikeDisp   = PreSortData.session(i).electrode(j).ons.meanSpike;
%             errSpikeDisp    = PreSortData.session(i).electrode(j).ons.errSpike;
%             figure(2);
%             subplot(441); hold on;
%             plot(meanSpikeDisp(q,:), '-o', 'Color', currColor, 'LineWidth', 2);
%             plot(meanSpikeDisp(q,:)+errSpikeDisp(q,:), '--', 'Color', currColor, 'LineWidth', 1);
%             plot(meanSpikeDisp(q,:)-errSpikeDisp(q,:), '--', 'Color', currColor, 'LineWidth', 1);
%         end
%         
%         disp('Compile event values for blind clustering.');
%         forClust = [eventsLocal.scalar.scl.values(:,1:4)';eventsLocal.gmono.dot.values(:,:)';eventsLocal.pca.dot.values(:,:)';eventsLocal.frequency.dot.values(:,:)'];
% %         forClust = [eventsLocal.scalar.scl.values(:,8:11)'];
%         
%         eventsLocal.forClust = forClust;
%         
%         if plotTemplates==1
%             figure(2);
%             subplot(441);
%             plot(1:size(eventsLocal.pca.template,1),eventsLocal.pca.template.*5000,'g-', 'LineWidth', 2)
%             plot(1:size(eventsLocal.frequency.template,1),eventsLocal.frequency.template.*500,'b-', 'LineWidth', 2)
% 
%             figure(50); clf; hold on;
%             for q=1:5
%                 plot(1:size(eventsLocal.pca.template,1),eventsLocal.pca.template(:,q),'Color',colorPalette(q,:),'LineWidth',2)
%             end        
%         end
%         
%         PreSortData.session(i).electrode(j).ons.events = eventsLocal;
% 
%     end
% 
%     blindClusts = 3;
%     recurClusts = 5;
%     [clustIds, centroids, dispersion, distances] = kmeans(forClust',blindClusts);
%     [values,sinds] = sort(clustIds,'ascend');
%         
%     figure(2);
%     subplot(445);
%     imagesc(corr(finalCorr));        
%     title('Online sorting based clustering');
% %     subplot(132);
% %     imagesc(corr(finalCorr2));        
% %     title('Online sorting based clustering');
%     subplot(446);
%     imagesc(corr(forClust(:,sinds)));
%     title('K-means clustering on quantified data');
% 
%     for i=1:blindClusts
%         indsBlClust = find(clustIds==i);
%         blindClustWfs(i,:) = mean(eventsLocal.Sh_wfs(indsBlClust,:),1);
%         figure(2);subplot(443);hold on;
%         plot(1:size(eventsLocal.Sh_wfs,2),blindClustWfs(i,:),'Color',colorPalette(i,:), 'LineWidth', 4);
%     end
%     
%     % recalculate distances using the first round of clustering and compare
%     clear round2vals;
%     for i=1:blindClusts
%         for j=1:size(eventsLocal.Sh_wfs,1)
%             round2vals(i,j) = pdist2(eventsLocal.Sh_wfs(j,:),blindClustWfs(i,:));
% %             round2vals(i,j) = max(xcorr(eventsLocal.Sh_wfs(j,:),blindClustWfs(i,:),50));
%         end
%     end    
%     
%     forClust2 = round2vals;
%     size(forClust)
%     size(forClust2)
%     [clustIds2, centroids2, dispersion2, distances2] = kmeans(forClust2',5);
%     [values2,sinds2] = sort(clustIds2,'ascend');
%     size(sinds2)
%     figure(2); subplot(447);
%     imagesc(corr(forClust2(:,sinds2)));
%     title('K-means clustering (recursion 1)');
%     
%     for i=1:recurClusts
%         indsBlClust2 = find(clustIds2==i);
%         blindClustWfs2(i,:) = mean(eventsLocal.Sh_wfs(indsBlClust2,:),1);
%         figure(2);subplot(444);hold on;
%         plot(1:size(eventsLocal.Sh_wfs,2),blindClustWfs2(i,:),'Color',colorPalette(i,:), 'LineWidth', 4);
%         figure(2); subplot(449); hold on;
%         plot3(forClust2(1,indsBlClust2),forClust2(2,indsBlClust2),forClust2(3,indsBlClust2),'.','Color',colorPalette(i,:), 'LineWidth', 4);
%         plot3(centroids2(i,1),centroids2(i,2),centroids2(i,3),'bo','MarkerSize', 10,'LineWidth',3);
%         view([35 25]);
%         figure(2); subplot(4,4,[10:12]); hold on;
%         plot3(eventsLocal.ts(indsBlClust2),eventsLocal.scalar.scl.values(indsBlClust2,4),eventsLocal.scalar.scl.values(indsBlClust2,1),'.','Color',colorPalette(i,:));
%     end
%     
% %% PROMOTE SORTED UNIT ARRAY TO ANALYZABLE DATA IN THE DATA STRUCTURE
% 
% %% LOAD NS5 DATA
% 
% %% FILTER AND GENERATE LF AND HF
% 
% %% EXTRACT BEHAVIOR / EVENT DATA
% 
% %% COMPRESS SEQ VIDEO DATA
% 
% %% MAKE TIMESTAMP ARRAY FOR VIDEO FRAMES
% 
% %% T0 >>> PROCESS LEVER DATA
% 
% % Load the behavioral data from analog channels
% leverX = 'ainp9';
% leverY = 'ainp10';
% 
% licking= 'ainp11';
% 
% eventA = 'ainp14';
% eventB = 'ainp15';
% 
% vidFrame= 'ainp13';
% vidFrame= 'ainp16';
% 
% % filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'ns4']
% filenamestrE = '20100429-m15-d1-006.ns4';
% Ns4DATA = openNSx('report','read',filenamestrE);
% 
% clear leverData sLeverData tmpLeverData 
% 
% rawXi = find(Ns4DATA.MetaTags.ChannelID==137);
% rawYi = find(Ns4DATA.MetaTags.ChannelID==138);
% rawLi = find(Ns4DATA.MetaTags.ChannelID==139);
% 
% leverData(1,:) = decimate(Ns4DATA.Data(rawXi,:),10);
% leverData(2,:) = decimate(Ns4DATA.Data(rawYi,:),10);
% 
% sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);
% sLeverData(2,:) = sgolayfilt(leverData(2,:),9,101);
% 
% tmpLick         = diff(sgolayfilt(Ns4DATA.Data(rawLi,:),9,101));
% lickData        = decimate(tmpLick,10);
% 
% behavior.sLeverData = sLeverData;
% behavior.lickData   = lickData;
% 
% rawData = decimate(Ns4DATA.Data(5,:),10);
% 
% % design the bandpass filter to use for field potentials
% dBT = fdesign.bandpass('n,fc1,fc2', 300, 4, 55, 1000);
% HdBT = design(dBT);   
% egLFP = filtfilt(HdBT.Numerator,1,rawData); % zero-phase filtering
% 
% behavior.egLFP      = egLFP;
% 
% %% T0 >>> LOAD the associated behavior file for the session
% % codeType = 2;
% % 
% % switch codeType
% %     
% %     case 2
% %         clear displace* trial* 
% % 
% %         fileNameStr = '_p.csv'
% %         trialParams = dlmread(fileNameStr,',',1,0);
% % 
% %         fileNameStr = '.csv'
% %         trialTimes = dlmread(fileNameStr,',',2,0);
% % 
% %         fileNameStr = 'tX.csv'
% %         tmpX = dlmread(fileNameStr,',',1,0);
% %         displaceX = tmpX(:,1:100);
% % 
% %         fileNameStr = 'tY.csv'
% %         tmpY = dlmread(fileNameStr,',',1,0);
% %         displaceY = tmpY(:,1:100);
% %         
% %     case 3
% %         
% % end
% 
% %% T1 >>> EXTRACT SPECTRAL DATA
% 
% %% T1 >>> POWER AS A FUNCTION OF TIME, FREQ CONSTANT
% 
% %% T1 >>> LFP STATE AS A FUNCTION OF TIME
% 
% % take spectral time slices
% 
% % find dimensions
% 
% % create time series data
% 
% %% T1 >>> ALIGNED PETH FOR LSTATE AND BPOW
% 
% %% T2 >>> EXTRACT EVENTS FROM HF DATA
% 
%     filenamestr = '20100429-m15-d1-006-01-sorted.nev';
%     dataEvents = openNEV(filenamestr,'read','nosave');
% 
%     MMdataSet.nevRes = dataEvents.MetaTags.SampleRes;
% 
%     rewardInds = find(dataEvents.Data.Spikes.Electrode==142);       % reward
%     MMdataSet.behavior.rewardInds = round(dataEvents.Data.Spikes.Timestamps(rewardInds)./30);
% 
%     threshInds = find(dataEvents.Data.Spikes.Electrode==143);       % movement
%     MMdataSet.behavior.threshInds = round(dataEvents.Data.Spikes.Timestamps(threshInds)./30);
% 
%     lickInds = find(dataEvents.Data.Spikes.Electrode==139);         % licking
%     MMdataSet.behavior.lickInds = round(dataEvents.Data.Spikes.Timestamps(lickInds)./30);
% 
%     numSpks = numel(dataEvents.Data.Spikes.Timestamps);
% 
%     [phys] = TNC_ExtractNEV2IndChan(dataEvents);
% 
%     causal = 0;
%     [popVec] = TNC_ExtractPopVector(phys,causal);
% 
% %% find popVec event times
% threshold = 3;
% popVecChan = 1;
% tmpPCA = abs(popVec.proj(popVecChan,:));
% suThr = find(tmpPCA>threshold);
% evInds = suThr(find(diff(suThr)>1)+1);
% 
% %% T2 >>> VISUALIZE RAW DATA
% 
% TNC_DisplayExampleWaveforms(dataEvents,numSpks);
% TNC_AnimatePopVector(popVec,[300 1400],5,3);
% 
% %% T2 >>> RECONCILE EVENTS FROM HF DATA
% 
% %% T2 >>> ALIGN EVENTS FROM HF DATA
% 
% %% T2 >>> EXTRACT EVENTS IMAGES FROM HF DATA
% 
% %% T2 >>> QUANTIFY IMAGE DISTANCES
% 
% %% T2 >>> CLASSIFY IMAGES
% 
% %% T3 >>> EVALUATION OF ONLINE SORTING FROM NEV
% 
% %% T3 >>> AUTOMATED, RECURSIVE SORTING OF NS5 DATA

