function [S,pc,sortedBinSpkCell] = pcaPSTHtwoPopulations_js2p0( filePath, fileName, varName, varargin  )
%pcaPSTH takes filePath, fileName, varName and other parameters to run
% a pca on trial-averaged z-scored PSTHs of units whose mean FR greater
% than the FRcut, and the units comprising top PCs selected by the user.
% pcaPSTH returns results of PCA and binned spike count structures/cells
% with the units sorted based on their PC scores of the top PCs.
% Modified on 9/6/18 to read in/create spikeCountMat temporarily (without
% saving) from the spikeTimes

cd(filePath)

%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919';
%fileName = {'binSpkCountWR40_081919.mat'}; 
%varName = 'rStartToPull'; 

p = parse_input_psthPCA( filePath, fileName, varName, varargin );
%p = parse_input_psthPCA( filePath, fileName, varName, {'binSize',100,'stepSize',100,'binSizeZ',50,'useAllUnits',false,'PCsLogic',true,'PCs',5,'expVarLogic',false,'expVarCut',80,'FRcut',.5,'nanTrialCut',.1,'cmap','hot','cAxis',[-3 10], 'rmvNaNtrials', true , 'rmvNaNunits', false} );

% imagesc of the PSTHs (potentially modify to use TNC_CreateRBColormap.m)
cmap = TNC_CreateRBColormapJP(100,p.Results.cmap); % generate a colormap for imagesc psth

%% preprocessing and PCA on trial-averaged z-scored PSTHs
numbTrialEachFile = zeros(length(fileName),1); % take the number of file in each file to ensure that files have the common number of trials
integratedNaNtrials = []; % collect NaN trials index of all the files
for f = 1:length(fileName) % increment files
    
    S = load(fileName{f},varName); % load the psth structure renamed as S
    S = S.(varName); % rename the structure as S (get rid of the original varName)
    
    numbTrialEachFile(f,1) = unique(cellfun(@length,S.currEvt)); % take the number of trials in each file
    
    if f == length(fileName)
        if length(unique(numbTrialEachFile)) ~= 1
            error('Files have different number of events!!!')
        end
    end
    
    % preprocess spike count matrices
    S.binSpk   = []; % trial x bin x units
    S.binSpkZ  = []; % trial x bin (50 ms, i.e., decimation factor 50)
    S.nanTrials  = []; % trial x units
    S.valFrCellId  = []; % valid cell id (cells passing the FR cut)
    
    
    valFrCellCnt = 0; % valid cell count
    for u = 1:length(S.SpkCountMatZ) % increment units
        S.binSpkZ(u,:)  = binAvg1msSpkCountMat( S.SpkCountMatZ{u}, p.Results.binSizeZ, p.Results.binSizeZ ); % bin/average the SpkCountMatZ - unit x bin (50 ms) z score
        %S.binSpkZ(u,:)  = decimate(S.SpkCountMatZ{u},p.Results.dcFactor); % decimate the SpkCountMatZ - unit x bin (50 ms) z score
        
        tempSpkCountMat = cell2mat(getSpkCntMatFromSpkTimes( S.SpkTimes{u}, S.params )); % get the current unit's spikeCountMat (trial-by-1msBin)
        % tempSpkCountMat = full(cell2mat(S.SpkCountMat{u}));     % temp spike count mat STR
        % imagesc(tempSpkCountMat) % image the unit's spikecount mat
        
        if nanmean(tempSpkCountMat(:))*1000 > p.Results.FRcut % in case the mean FR greater than 1Hz
            valFrCellCnt = valFrCellCnt + 1;
            S.valFrCellId(valFrCellCnt,1) = u; % valid unit ID (surpassing the FR cut)
            tempBinSpk      = bin1msSpkCountMat(tempSpkCountMat,p.Results.binSize,p.Results.stepSize); % 1 ms binSize, 1 ms stepSize
            S.binSpk(:,:,valFrCellCnt) = tempBinSpk; % trial x bin x unit (raw spike counts)
            S.nanTrials(:,valFrCellCnt) = isnan(sum(tempBinSpk,2)); % trial x units (detect nan trials - nan due to an event(s) out of a unit's spike range)
        else
            
        end
        
    end
    clearvars u
    
    %% run PCA on the trial-averaged z-scored PSTHs and sort units based on their PC scores
    [pc.loadings, pc.PCscore, pc.eigValues, ~, pc.expVar] = pca(S.binSpkZ(S.valFrCellId(:,1),:)); % run pca on trial-averaged z scores (include all units passing the FR cut)
    [~,pc.maxDim] = max(pc.PCscore,[],2); % dimension for each unit whose PC score is greatest
    
    pc.pcUnitSort = [];
    pcMaxDims = unique(pc.maxDim);  % pc max dimensions
    S.valFrCellId(:,2) = [1:1:size(S.valFrCellId,1)]; % valid cell id before sort
    S.valFrCellId(:,3) = pc.maxDim; % put pc max dimensions
    
    % sort cells by their PC scores
    for i = 1:length(pcMaxDims) % increment pc max dims
        tempUnitList = S.valFrCellId(S.valFrCellId(:,3)==pcMaxDims(i),:); % valid unit id
        tempUnitList(:,4) = pc.PCscore(tempUnitList(:,2),pcMaxDims(i)); % their max PC score (on that dimension)
        srtTempUnitList = sortrows(tempUnitList,-4); % sortRows by the PC score in a descending manner
        pc.pcUnitSort = [pc.pcUnitSort; srtTempUnitList];  % concatenate the unit list
    end
    clearvars i temp*
    
    % select units with max PCs corresponding to the user-defined top PCs
    if p.Results.useAllUnits
        pc.pcUnitCutLogic = ones(size(pc.pcUnitSort,1),1); % use all units
        
    elseif p.Results.PCsLogic
        pc.pcUnitCutLogic = pc.pcUnitSort(:,3)<=p.Results.PCs;
        pc.pcUnitSort = pc.pcUnitSort(pc.pcUnitCutLogic,:); % select units whose max PC scores are of the set top PCs
    elseif p.Results.expVarLogic
        sumExpVar = zeros(length(pc.expVar),1);
        for i = 1:length(pc.expVar)
            sumExpVar(i) = sum(pc.expVar(1:i)); % cumulative explained variance
        end
        pc.pcUnitCutLogic = pc.pcUnitSort(:,3)<=find(sumExpVar>=p.Results.expVarCut,1);
        pc.pcUnitSort = pc.pcUnitSort(pc.pcUnitCutLogic,:); % select units whose max PC scores are of the top PCs whose sum explained variance greater than the set expVarCut
    else
        error('Determine whether to select units based on their PCs for DCA/GPFA!!')
    end
    
    %% put the relevant structures S and pc in to the structure named file
    file(f).S  = S;
    file(f).pc = pc;
    integratedNaNtrials = [integratedNaNtrials, file(f).S.nanTrials];
    clearvars S pc pcMaxDims
    
end
clearvars f

%% eliminate nan trials either by eliminating trials or units so that the resulting S.binSpkVal has no NaN trials
if p.Results.rmvNaNtrials % in case dealing NaN trials by excluding trials
    missingTrialLimit = floor(unique(numbTrialEachFile)*p.Results.nanTrialCut); % missing trial limit (default: 10 % of total trials)
    if sum(sum(integratedNaNtrials,2)>0)<missingTrialLimit % in case there are not so many missing trials (across units)
        nonNaNtrialId = find(sum(integratedNaNtrials,2)==0);  % find trials at which no unit had a NaN trial
        for f = 1:length(fileName)
            file(f).S.nonNaNtrialId = nonNaNtrialId; % valid trials at which none of the units had a NaN trial
            file(f).S.binSpkVal = file(f).S.binSpk(nonNaNtrialId,:,:); % eliminate NaN trials binSpkVal (trial x timebin x unit) must be free of any NaN values
            % organize the binSpkVal into a cell array
            permbinSpk = permute(file(f).S.binSpkVal,[3 1 2]); % permute dimensions so that binSpk to be unit x trials x timebin
            sortedBinSpkCell = cell(1,size(permbinSpk,3));     % each cell contains data for each time bin
            for b = 1:length(sortedBinSpkCell) % increment bins
                file(f).sortedBinSpkCell{1,b} = permbinSpk(file(f).pc.pcUnitSort(:,2),:,b); % each cell contains data for each time bin (unit x trials)
            end
            clearvars b
        end
        clearvars f
    else
        error('There might be a unit(s) with too many NaN trials!!')
    end
    
elseif p.Results.rmvNaNunits % in case dealing NaN trials remove NaN units
    for f = 1:length(fileName)
        if sum(sum(file(f).S.nanTrials,1)>0)< floor(size(file(f).S.nanTrials,2)*p.Results.nanUnitCut) % if there are units with NaN trials fewer than the allowed proportion (p.Results.nanUnitCut, e.g. 10%)
            file(f).S.unitNaNidx = sum(file(f).S.nanTrials,1)>0;        % idx for units with 1 or more nan trials, to be deselected below
            file(f).S.binSpkVal  = file(f).S.binSpk(:,:,~file(f).S.unitNaNidx & file(f).pc.pcUnitCutLogic'); % select units satisfying the PC criterion and w/o NaN trials binSpkVal must be free of any NaN values
            % organize the selected cells into a cell array
            permbinSpk = permute(file(f).S.binSpkVal,[3 1 2]); % permute dimensions so that binSpk to be unit x trials x timebin
            sortedBinSpkCell = cell(1,size(permbinSpk,3));     % each cell contains data for each time bin
            
            for b = 1:length(sortedBinSpkCell) % increment bins
                file(f).sortedBinSpkCell{1,b} = permbinSpk(:,:,b); % each cell contains data for each time bin (unit x trials)
            end
            clearvars b
            
        else % if there are too many units with NaN trials to be removed
            error('prog:input','File %d has as many as %d units with NaN trials to be removed!!!',f,sum(sum(file(f).S.nanTrials,1)>0));
        end
    end
else
    error('Determine how to deal with NaN trials for DCA/GPFA!!')
end

%% Generate plots and Save Data
for f = 1:length(fileName)
    % imagesc the trial-averaged z-score psths of valid units (subjected to PCA)
    Xs = file(f).S.params.binEdges1ms(1):p.Results.binSizeZ:file(f).S.params.binEdges1ms(end); 
    figure; 
    imagescJP(file(f).S.binSpkZ(file(f).pc.pcUnitSort(:,1),:),cmap,p.Results.cAxis,'X',Xs,'cBar',true); 
    pbaspect([1 1 1]); % plot the trial-averaged z-scored PSTHs
    xlim(p.Results.imagescXlim)
    print( fullfile(filePath, 'Figure', [fileName{f},'_trialAvgZscore',varName]),'-dpdf','-bestfit','-painters')
    
    % plot relevant pc loadings
    pcLoadings = unique(file(f).pc.pcUnitSort(:,3)); % pcLoadins of the top PCs
    figure; hold on;
    for i = 1:length(pcLoadings)
        plot(file(f).pc.loadings(:,pcLoadings(i)));
        pcLoadingLabel{i} = strcat('pc',num2str(pcLoadings(i)));
    end
    hold off; pbaspect([1 1 1]); clearvars i
    xlim([41, 80])
    legend(pcLoadingLabel);
    print( fullfile(filePath, 'Figure', [fileName{f},'_trialAvgPsthPCs',varName]),'-dpdf','-bestfit','-painters')
    
    %% Save variables
    cd(filePath)
    saveName = strcat(fileName{f},'pcaPSTH',varName,num2str(p.Results.binSize),'ms');
    S = file(f).S;
    pc = file(f).pc;
    sortedBinSpkCell = file(f).sortedBinSpkCell;
    
    save(saveName,'sortedBinSpkCell','pc','S','p');
    %clearvars S pc sortedBinSpkCell
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = parse_input_psthPCA( filePath, fileName, varName, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
%parse input, and extract name-value pairs for the main function
% 'pcaPSTH.m'

default_binSize  = 200; % 200 ms bins
default_stepSize = default_binSize; % by default, use the stepSize same as the binSize
default_binSizeZ = 50;  % binSize for SpkCountMatZ
default_PCs      = 3;   % use top 3 PCs
default_expVarCut = 80; % to include top PCs whose summed explained variance greater than 80% of total variance
default_FRcut = 1; % exclude units with mean FR lower than the FRcut from PCA
default_nanTrialCut = .1; % to exclude NaN trials, if the proportion of them relative to total trials is smaller than this proportion
default_nanUnitCut = .1;  % to exclude units with NaN trials, if the proportion of them relative to total units is smaller than this proportion
default_useAllUnits = true;  % logical to include all units in the sortedBinSpkCell prepared for DCA
default_PCsLogic = false;    % logical to choose to include only the top PCs units
default_expVarLogic = false; % logical to choose to include units with PCs of which cumulative explained variance surpassing the set expVarCut (e.g. 80 %)
default_cmap = 'cb';    % default colormap
default_cAxis = [-3 3]; % default colorAxis
default_rmvNaNunits = false; % exclude units with NaN trials from PCA and DCA preprocessing
default_rmvNaNtrials = true; % remove the trials in which one or more units have NaN trials
default_imagescXlim = [-2000 2000]; % default xlim to be used for imagesc plot 1ms 

p = inputParser; % create parser object
addRequired(p,'filePath');
addRequired(p,'fileName');
addRequired(p,'varName');

addParameter(p,'binSize', default_binSize)
addParameter(p,'stepSize',default_stepSize)
addParameter(p,'binSizeZ', default_binSizeZ)
addParameter(p,'PCs', default_PCs)
addParameter(p,'expVarCut', default_expVarCut)
addParameter(p,'FRcut', default_FRcut)
addParameter(p,'nanTrialCut', default_nanTrialCut)
addParameter(p,'nanUnitCut', default_nanUnitCut)
addParameter(p,'useAllUnits',default_useAllUnits)
addParameter(p,'PCsLogic',default_PCsLogic)
addParameter(p,'expVarLogic',default_expVarLogic)
addParameter(p,'cmap',default_cmap)
addParameter(p,'cAxis',default_cAxis)
addParameter(p,'rmvNaNunits',default_rmvNaNunits)
addParameter(p,'rmvNaNtrials',default_rmvNaNtrials)
addParameter(p,'imagescXlim',default_imagescXlim)

parse(p, filePath, fileName, varName, vargs{:})

end

end % end of the function
