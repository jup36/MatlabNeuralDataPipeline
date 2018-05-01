function [S,pc,sortedBinSpkCell] = pcaPSTH( filePath, fileName, varName, varargin  )
%pcaPSTH takes filePath, fileName, varName and other parameters to run 
% a pca on trial-averaged z-scored PSTHs of units whose mean FR greater
% than the FRcut, and the units comprising top PCs selected by the user. 
% pcaPSTH returns results of PCA and binned spike count structures/cells 
% with the units sorted based on their PC scores of the top PCs. 

cd(filePath)

p = parse_input_psthPCA( filePath, fileName, varName, varargin );

S = load(fileName,varName); % load the psth structure renamed as S
S = S.(varName); % rename the structure simply as S (get rid of the original varName)

% imagesc of the PSTHs (potentially modify to use TNC_CreateRBColormap.m)
cmap = TNC_CreateRBColormap(100,p.Results.cmap); % generate a colormap for imagesc psth

%% preprocessing
Xs = cell(1,2); % Xs comprises two neural datasets: M1 and DMS 

%   Xs: (1 x M cell), data, M datasets where Xs{iset} is num_variables x num_samples

% STR cell pre-processing
S.binSpk   = [];  % trial x bin x units
S.binSpkZ  = [];  % trial x bin (50 ms, i.e., decimation factor 50)
nanTrials  = []; % trial x units
valCellId  = []; % valid cell id

valCellCnt = 0; % valid cell count

for u = 1:length(S.SpkCountMatZ) % increment units

    S.binSpkZ(u,:)    = decimate(S.SpkCountMatZ{u},p.Results.dcFactor); % decimate the SpkCountMatZ - unit x bin (50 ms) z score
    tempSpkCountMat = full(cell2mat(S.SpkCountMat{u}));     % temp spike count mat STR    
    
    if nanmean(tempSpkCountMat(:))*1000 > p.Results.FRcut % in case no missing trial & mean FR greater than 1Hz    
        valCellCnt = valCellCnt + 1;
        valCellId(valCellCnt,1) = u; % valid unit ID (surpassing the FR cut)
        tempBinSpk      = bin1msSpkCountMat(tempSpkCountMat,p.Results.binSize,p.Results.stepSize); % 200 ms binSize, 200 ms stepSize
        S.binSpk(:,:,valCellCnt) = tempBinSpk; % trial x bin x unit (raw spike counts)
        nanTrials(:,valCellCnt) = isnan(sum(tempBinSpk,2)); % trial x units (detect nan trials - nan due to an event(s) out of a unit's spike range)
    else
        
    end

end
clearvars u 

% eliminate nan trials
missingTrialLimit = floor(size(nanTrials,1)*p.Results.nanTrialCut);  % missing trial limit (default: 10 % of total trials)
unitNaNidx1  = sum(nanTrials,1)>missingTrialLimit;          % idx for units with too many nan trials, to be deselected below
unitNaNidx2  = sum(nanTrials,1)~=mode(sum(nanTrials,1));    % idx for units with uncommon number of nan trials, to be deselected below
unitNaNidx = unitNaNidx1 | unitNaNidx2;

binSpkValunits = S.binSpk(:,:,~unitNaNidx);  % eliminate units with too many nan trials
valCellId = valCellId(~unitNaNidx,1);        % update valCelIdx accordingly
trialNaNidx = sum(nanTrials(:,~unitNaNidx),2)>=1; % idx for NaN trials (delete the trial if there's any unit with NaN value on that trial) 
S.binSpkVal = binSpkValunits(~trialNaNidx,:,:);   % eliminate NaN trials binSpkVal must be free of any NaN values

%imagescJP(S.binSpkZ,hotcolormap,[-2 4]) % colormap for sanity check

%% run PCA and sort units based on their PC scores
[pc.loadings, pc.PCscore, pc.eigValues, ~, pc.expVar] = pca(S.binSpkZ(valCellId,:)); % run pca on trial-averaged z scores
[~,pc.maxDim] = max(pc.PCscore,[],2); % dimension for each unit whose PC score is greatest

pc.pcUnitSort = [];
pcMaxDims = unique(pc.maxDim); % pc max dimensions 
valCellId(:,2) = [1:1:size(valCellId,1)]; % valid cell id before sort
valCellId(:,3) = pc.maxDim;    % put pc max dimensions

for i = 1:length(pcMaxDims) % increment pc max dims
    tempUnitList = valCellId(find(valCellId(:,3)==pcMaxDims(i)),:); % valid unit id 
    tempUnitList(:,4) = pc.PCscore(tempUnitList(:,2),pcMaxDims(i)); % their max PC score (on that dimension)
    srtTempUnitList = sortrows(tempUnitList,-4); % sortRows by the PC score in a descending manner
    pc.pcUnitSort = [pc.pcUnitSort; srtTempUnitList];  % concatenate the unit list 
end
clearvars i temp*

% select units with max PCs corresponding to the user-defined top PCs 
if p.Results.useAllUnits
   
elseif p.Results.PCsLogic
   pc.pcUnitSort = pc.pcUnitSort(pc.pcUnitSort(:,3)<=p.Results.PCs,:); % select units whose max PC scores are of the set top PCs
elseif p.Results.expVarLogic
   sumExpVar = zeros(length(pc.expVar),1);
   for i = 1:length(pc.expVar) 
        sumExpVar(i) = sum(pc.expVar(1:i)); % cumulative explained variance
   end
   pc.pcUnitSort = pc.pcUnitSort(pc.pcUnitSort(:,3)<=find(sumExpVar>=p.Results.expVarCut,1),:); % select units whose max PC scores are of the top PCs whose sum explained variance greater than the set expVarCut 
else
    
end


%% Organize the selected cells into a cell array
permbinSpk = permute(S.binSpkVal,[3 1 2]);        % permute dimensions so that binSpk to be unit x trials x timebin 
sortedBinSpkCell = cell(1,size(permbinSpk,3));  % each cell contains data for each time bin  

for b = 1:length(sortedBinSpkCell) % increment bins
    sortedBinSpkCell{1,b} = permbinSpk(pc.pcUnitSort(:,2),:,b); % each cell contains data for each time bin (unit x trials)
end
clearvars b

%% Generate plots
% figure; imagescJP(squeeze(nanmean(permbinSpk(pc.pcaCellIdx,:,:),2)),hotcolormap,[-2 5]) % colormap for sanity check
% imagesc the trial-averaged z-score psths of valid units (subjected to PCA)
imagescJP(S.binSpkZ(pc.pcUnitSort(:,1),:),cmap,p.Results.cAxis); pbaspect([1 1 1]); 
% plot relevant pc loadings
pcLoadings = unique(pc.pcUnitSort(:,3)); % pcLoadins of the top PCs
figure; hold on; 
for i = 1:length(pcLoadings)
    plot(pc.loadings(:,pcLoadings(i))); 
    pcLoadingLabel{i} = strcat('pc',num2str(pcLoadings(i)));
end
hold off; pbaspect([1 1 1]); clearvars i
legend(pcLoadingLabel);

%% Save variables
cd(filePath)
saveName = strcat(fileName,'pcaPSTH',varName,num2str(p.Results.binSize),'ms');
save(saveName,'sortedBinSpkCell','pc','S','p');

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_psthPCA( filePath, fileName, varName, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
        %parse input, and extract name-value pairs for the main function
        % 'pcaPSTH.m'
        
        default_binSize  = 200; % 200 ms bins
        default_stepSize = default_binSize; % by default, use the stepSize same as the binSize
        default_dcFactor = 50;  % decimation factor
        default_PCs      = 3;   % use top 3 PCs
        default_expVarCut = 80; % to include top PCs whose summed explained variance greater than 80% of total variance
        default_FRcut = 1; % exclude units with mean FR lower than the FRcut from PCA
        default_nanTrialCut = .1; % to exclude units of which NaN trial proportion greater than this  
        default_useAllUnits = true;  % logical to include all units in the sortedBinSpkCell prepared for DCA
        default_PCsLogic = false;    % logical to choose to include only the top PCs units
        default_expVarLogic = false; % logical to choose to include units with PCs of which cumulative explained variance surpassing the set expVarCut (e.g. 80 %)
        default_cmap = 'cb'; % default colormap
        default_cAxis = [-3 3];
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addRequired(p,'fileName');
        addRequired(p,'varName');
       
        addParameter(p,'binSize', default_binSize)
        addParameter(p,'stepSize',default_stepSize)
        addParameter(p,'dcFactor', default_dcFactor)
        addParameter(p,'PCs', default_PCs)
        addParameter(p,'expVarCut', default_expVarCut)
        addParameter(p,'FRcut', default_FRcut)
        addParameter(p,'nanTrialCut', default_nanTrialCut)
        addParameter(p,'useAllUnits',default_useAllUnits)
        addParameter(p,'PCsLogic',default_PCsLogic)
        addParameter(p,'expVarLogic',default_expVarLogic)
        addParameter(p,'cmap',default_cmap)       
        addParameter(p,'cAxis',default_cAxis)
        
        parse(p, filePath, fileName, varName, vargs{:})
        
    end

end

