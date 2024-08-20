
filePath = '/Volumes/Extreme SSD/js2p0/gpfa_model';
load(fullfile(filePath, 'simDataGpfa.mat'), 'trialgroup')

%% preprocessing
trialgroup_spec_C = cell(4, 1);
for i = 1:length(trialgroup.spec)
    trialgroup_spec_C{1}{end+1}= bin1msSpkCountMat(max(trialgroup.spec(i).psth(:, :, 1), 0), 50, 50); % 50-ms bin
    trialgroup_spec_C{2}{end+1}= bin1msSpkCountMat(max(trialgroup.spec(i).psth(:, :, 2), 0), 50, 50);
    trialgroup_spec_C{3}{end+1}= bin1msSpkCountMat(max(trialgroup.spec(i).psth(:, :, 3), 0), 50, 50);
    trialgroup_spec_C{4}{end+1}= bin1msSpkCountMat(max(trialgroup.spec(i).psth(:, :, 4), 0), 50, 50);
end

trialgroup_selt_C = cell(4, 1);
for i = 1:length(trialgroup.selt)
    trialgroup_selt_C{1}{end+1}= bin1msSpkCountMat(max(trialgroup.selt(i).psth(:, :, 1), 0), 50, 50); % 50-ms bin
    trialgroup_selt_C{2}{end+1}= bin1msSpkCountMat(max(trialgroup.selt(i).psth(:, :, 2), 0), 50, 50);
    trialgroup_selt_C{3}{end+1}= bin1msSpkCountMat(max(trialgroup.selt(i).psth(:, :, 3), 0), 50, 50);
    trialgroup_selt_C{4}{end+1}= bin1msSpkCountMat(max(trialgroup.selt(i).psth(:, :, 4), 0), 50, 50);
end

% get datSpec
count = 0;
for j = 1:4
    for jj = 1:10
        count = count + 1;
        datSpec(count).trialId = count;
        datSpec(count).spikes = max(trialgroup_spec_C{j}{jj}, 0);
    end
end

% get datSel
count = 0;
for j = 1:4
    for jj = 1:10
        count = count + 1;
        datSelt(count).trialId = count;
        datSelt(count).spikes = max(trialgroup_selt_C{j}{jj}, 0);
    end
end

%% GPFA specification model
% specification model
saveNameSpec = fullfile(filePath, 'gpfa_spec');
if exist(saveNameSpec, "dir")~=7
    mkdir(saveNameSpec)
end

if exist(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')), 'file')==2
    gpfaRezSpec = load(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')));
else
    gpfaRezSpec = neuralTrajAlreadyBinned(saveNameSpec, datSpec, 'method', 'gpfa', ...
    'xDim', 5, 'kernSDList', 5, 'binWidth', 1); %
end

% gpfa postprocess
gpfaRezSpec.method = 'gpfa';
% Orthonormalize neural trajectories
[gpfaRezSpec.estParamsOrtho, gpfaRezSpec.seqTrainOrtho] = postprocess(gpfaRezSpec, 'kernSD', 5);

gpfaRezSpec_tt1 = projYtoGetXsm(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSpec(1:10));  % trial type 1
gpfaRezSpec_tt2 = projYtoGetXsm(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSpec(11:20)); % trial type 2
gpfaRezSpec_tt3 = projYtoGetXsm(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSpec(21:30)); % trial type 3
gpfaRezSpec_tt4 = projYtoGetXsm(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSpec(31:40)); % trial type 4

plot3D_seqCells_color({gpfaRezSpec_tt1, gpfaRezSpec_tt2, gpfaRezSpec_tt3, gpfaRezSpec_tt4}, {[0 0 1], [0 1 1], [1 0 1], [1 0 0]}, 10);

%% GPFA selection model
saveNameSelt = fullfile(filePath, 'gpfa_selt');
if exist(saveNameSelt, "dir")~=7
    mkdir(saveNameSelt)
end

if exist(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')), 'file')==2
    gpfaRezSelt = load(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')));
else
    gpfaRezSelt = neuralTrajAlreadyBinned(saveNameSelt, datSelt, 'method', 'gpfa', ...
    'xDim', 5, 'kernSDList', 5, 'binWidth', 1); %
end

% gpfa postprocess
gpfaRezSelt.method = 'gpfa';
% Orthonormalize neural trajectories
[gpfaRezSelt.estParamsOrtho, gpfaRezSelt.seqTrainOrtho] = postprocess(gpfaRezSelt, 'kernSD', 5);

gpfaRezSelt_tt1 = projYtoGetXsm(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSelt(1:10));  % trial type 1
gpfaRezSelt_tt2 = projYtoGetXsm(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSelt(11:20)); % trial type 2
gpfaRezSelt_tt3 = projYtoGetXsm(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSelt(21:30)); % trial type 3
gpfaRezSelt_tt4 = projYtoGetXsm(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSelt(31:40)); % trial type 4

plot3D_seqCells_color({gpfaRezSelt_tt1, gpfaRezSelt_tt2, gpfaRezSelt_tt3, gpfaRezSelt_tt4}, {[0 0 1], [0 1 1], [1 0 1], [1 0 0]}, 10);


%% PCA selection model
saveNameSeltPca = fullfile(filePath, 'pca_selt');
if exist(saveNameSeltPca, "dir")~=7
    mkdir(saveNameSeltPca)
end

if exist(fullfile(saveNameSeltPca, strcat(sprintf('%s_xDim%02d', 'pca', 5), '.mat')), 'file')==2
    pcaRezSelt = load(fullfile(saveNameSeltPca, strcat(sprintf('%s_xDim%02d', 'pca', 5), '.mat')));
else
    pcaRezSelt = neuralTrajAlreadyBinned(saveNameSeltPca, datSelt, 'method', 'pca', ...
    'xDim', 5, 'kernSDList', 5, 'binWidth', 1); %
end

pcaXpost1 = zeros(3, 50, 10); 
count = 0; 
for j = 1:10
    count = count+1; 
    pcaXpost1(:, :, count) = pcaRezSelt.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost2 = zeros(3, 50, 10); 
count = 0; 
for j = 11:20
    count = count+1;
    pcaXpost2(:, :, count) = pcaRezSelt.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost3 = zeros(3, 50, 10); 
count = 0; 
for j = 21:30
    count = count+1;
    pcaXpost3(:, :, count) = pcaRezSelt.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost4 = zeros(3, 50, 10); 
count = 0; 
for j = 31:40
    count = count+1;
    pcaXpost4(:, :, count) = pcaRezSelt.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpostC = {pcaXpost1, pcaXpost2, pcaXpost3, pcaXpost4}; 

plot3D_pcaSeqCells_color(pcaXpostC, {[0 0 1], [0 1 1], [1 0 1], [1 0 0]}, 10);

%% PCA specification model
saveNameSpecPca = fullfile(filePath, 'pca_Spec');
if exist(saveNameSpecPca, "dir")~=7
    mkdir(saveNameSpecPca)
end

if exist(fullfile(saveNameSpecPca, strcat(sprintf('%s_xDim%02d', 'pca', 5), '.mat')), 'file')==2
    pcaRezSpec = load(fullfile(saveNameSpecPca, strcat(sprintf('%s_xDim%02d', 'pca', 5), '.mat')));
else
    pcaRezSpec = neuralTrajAlreadyBinned(saveNameSpecPca, datSpec, 'method', 'pca', ...
    'xDim', 5, 'kernSDList', 5, 'binWidth', 1); %
end

pcaXpost1 = zeros(3, 50, 10); 
count = 0; 
for j = 1:10
    count = count+1; 
    pcaXpost1(:, :, count) = pcaRezSpec.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost2 = zeros(3, 50, 10); 
count = 0; 
for j = 11:20
    count = count+1;
    pcaXpost2(:, :, count) = pcaRezSpec.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost3 = zeros(3, 50, 10); 
count = 0; 
for j = 21:30
    count = count+1;
    pcaXpost3(:, :, count) = pcaRezSpec.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost4 = zeros(3, 50, 10); 
count = 0; 
for j = 31:40
    count = count+1;
    pcaXpost4(:, :, count) = pcaRezSpec.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpostC = {pcaXpost1, pcaXpost2, pcaXpost3, pcaXpost4}; 

plot3D_pcaSeqCells_color(pcaXpostC, {[0 0 1], [0 1 1], [1 0 1], [1 0 0]}, 10);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq = projYtoGetXsm(pathToGpfaRez, dat)
%This utility function handles projection of test trials (each trial dimension: NxT)
% onto the dimensionality-reduction subspace (e.g., gpfa).
% "exactInferenceWithLL.m" is the original code that handles the projection.
%       dif      = bsxfun(@minus, [seq(nList).y], estParams.d); % yDim x sum(T)
%       term1Mat = reshape(CRinv * dif, xDim*T, []); % (xDim*T) x length(nList)

% pathToGpfaRez = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/gpfa_ctx/gpfa_xDim05';
% estParams = gpfaRezCtx.estParams;
% dat: 1xTrials struct array with two fields - trialId, spikes.

% load parameters
load(fullfile(pathToGpfaRez), 'estParams', 'hasSpikesBool');

% get seq and ensure the y-dim matches with C-dim
seq = getSeq(dat, 1);
[yDim, xDim] = size(estParams.C);

if size(estParams.C, 1)==sum(hasSpikesBool)
    if any(hasSpikesBool==false)
        for t = 1:length(seq)
            seq(t).y = seq(t).y(hasSpikesBool, :);
        end
    end
else
    error("Dimensions of C and dat mismatch!")
end

%% Main computations to get the xsm
if estParams.notes.RforceDiagonal
    Rinv     = diag(1./diag(estParams.R));
    logdet_R = sum(log(diag(estParams.R)));
else
    Rinv     = inv(estParams.R);
    Rinv     = (Rinv+Rinv') / 2; % ensure symmetry
    logdet_R = logdet(estParams.R);
end
CRinv  = estParams.C' * Rinv;
CRinvC = CRinv * estParams.C;

Tall = [seq.T];
Tu   = unique(Tall);

T = Tu;

[K_big, K_big_inv, logdet_K_big] = make_K_big(estParams, T);

% There are three sparse matrices here: K_big, K_big_inv, and CRinvC_inv.
% Choosing which one(s) to make sparse is tricky.  If all are sparse,
% code slows down significantly.  Empirically, code runs fastest if
% only K_big is made sparse.
%
% There are two problems with calling both K_big_inv and CRCinvC_big
% sparse:
% 1) their sum is represented by Matlab as a sparse matrix and taking
%    its inverse is more costly than taking the inverse of the
%    corresponding full matrix.
% 2) term2 has very few zero entries, but Matlab will represent it as a
%    sparse matrix.  This makes later computations with term2 ineffficient.

K_big = sparse(K_big);

blah        = cell(1, T);
[blah{:}]   = deal(CRinvC);
%CRinvC_big = blkdiag(blah{:});     % (xDim*T) x (xDim*T)
[invM, logdet_M] = invPerSymm(K_big_inv + blkdiag(blah{:}), xDim,...
    'offDiagSparse', true);

% Note that posterior covariance does not depend on observations,
% so can compute once for all trials with same T.
% xDim x xDim posterior covariance for each timepoint
Vsm = nan(xDim, xDim, T);
idx = 1: xDim : (xDim*T + 1);
for t = 1:T
    cIdx       = idx(t):idx(t+1)-1;
    Vsm(:,:,t) = invM(cIdx, cIdx);
end

% T x T posterior covariance for each GP
VsmGP = nan(T, T, xDim);
idx   = 0 : xDim : (xDim*(T-1));
for i = 1:xDim
    VsmGP(:,:,i) = invM(idx+i,idx+i);
end

% Process all trials with length T
nList    = find(Tall == T);
dif      = bsxfun(@minus, [seq(nList).y], estParams.d); % yDim x sum(T)
term1Mat = reshape(CRinv * dif, xDim*T, []); % (xDim*T) x length(nList)

% Compute blkProd = CRinvC_big * invM efficiently
% blkProd is block persymmetric, so just compute top half
Thalf   = ceil(T/2);
blkProd = zeros(xDim*Thalf, xDim*T);
idx     = 1: xDim : (xDim*Thalf + 1);
for t = 1:Thalf
    bIdx            = idx(t):idx(t+1)-1;
    blkProd(bIdx,:) = CRinvC * invM(bIdx,:);
end
blkProd = K_big(1:(xDim*Thalf), :) *...
    fillPerSymm(speye(xDim*Thalf, xDim*T) - blkProd, xDim, T);
xsmMat  = fillPerSymm(blkProd, xDim, T) * term1Mat; % (xDim*T) x length(nList)

ctr = 1;
for n = nList
    seq(n).xsm   = reshape(xsmMat(:,ctr), xDim, T);
    seq(n).Vsm   = Vsm;
    seq(n).VsmGP = VsmGP;

    ctr = ctr + 1;
end

%% orthonormalize
[Xorth, Corth] = orthogonalize([seq.xsm], estParams.C);
seq = segmentByTrial(seq, Xorth, 'xorth');
estParams.Corth = Corth;

end
