function runDCA( filePath, X1fileName, X2fileName, numbTrialShuffle, saveName )
%runDCA performs distance covariance analysis using 'dca.m' from Byron's
% group. 
% Trial shuffle control data is generated using only the matching time
% bins. 

cd(filePath)

X1 = load(X1fileName); % X1 structure 
X2 = load(X2fileName); % X2 structure 
    
% Sanity check
if isequal(length(X1.sortedBinSpkCell), length(X2.sortedBinSpkCell))
    numbBin = length(X1.sortedBinSpkCell); % the number of spike count time bins
else
    error('Error: The number of spike count bins of the two variables must match!')
end

if isequal(size(X1.sortedBinSpkCell{1},2), size(X2.sortedBinSpkCell{1},2))
    numbTrial = size(X1.sortedBinSpkCell{1},2); % the number of trials across which the population 
else
    error('Error: The number of trials the two variables must match!')
end

dcaResults.dCovs = cell(numbBin,numbBin); % timeBin x timeBin distance covariance matrix
dcaResults.Us1 = cell(numbBin,numbBin);   % linear projections for population 1 
dcaResults.Us2 = cell(numbBin,numbBin);   % linear projections for population 2
dcaResults.Ps  = cell(numbBin,numbBin);   % parameters

%% run dca 
formatSpec = 'Processed timeBin %d\n of timeBin %d\n'; 
for i = 1:numbBin % increment time bins
    Xs{1} = X1.sortedBinSpkCell{1,i}; % take the neuron x trial spkc mat of the current time bin
    for ii = 1:numbBin % increment time bins    
        Xs{2} = X2.sortedBinSpkCell{1,ii}; % take the neuron x trial spkc mat of the current time bin
        [U, dcaResults.dCovs{i,ii}, p] = dca(Xs, 'num_dca_dimensions', 5, 'percent_increase_criterion', 0.01);
        dcaResults.Us1{i,ii} = U{1,1}; % linear projs of current time bin combination of population 1 
        dcaResults.Us2{i,ii} = U{1,2}; % linear projs of current time bin combination of population 2
        dcaResults.Ps{i,ii} = p; % store parameters 
        clearvars U p    
        fprintf(formatSpec, ii, i) % report progress
    end
end
clearvars i ii

for i = 1:numbBin % increment time bins
    dcaResults.dCovsDiag(1,i) = dcaResults.dCovs{i,i}(1); % dCovs from the matching time bins (diagonal dCovs)
end

%% run dca with trial shuffling (only run using matching time bins)
formatSpec2 = 'Processed timeBin %d\n of shuffle run %d\n'; 
for ts = 1:numbTrialShuffle % increment the number of trial shuffling to be done
    shufTrial1 = randsample(numbTrial,numbTrial); % random sample for X1
    shufTrial2 = randsample(numbTrial,numbTrial); % random sample for X2
    
    c{1,ts} = 'shufdCovs'; % just to label trial-shuffled dCovs to plot with gramm
    
    for i = 1:numbBin
        XsShuf{1} = X1.sortedBinSpkCell{1,i}(:,shufTrial1); % take the trial-shuffled sc mat of the current bin
        XsShuf{2} = X2.sortedBinSpkCell{1,i}(:,shufTrial2); % take the trial-shuffled sc mat of the current bin
        [Ushuf, dcaResults.dCovsShuf{ts,i}, p] = dca(XsShuf, 'num_dca_dimensions', 1, 'percent_increase_criterion', 0.01);
        dcaResults.Us1Shuf{ts,i} = Ushuf{1,1}; % linear projs of current time bin combination of population 1 
        dcaResults.Us2Shuf{ts,i} = Ushuf{1,2}; % linear projs of current time bin combination of population 2
        dcaResults.PsShuf{ts,i} = p; % store parameters 
        clearvars Ushuf p    
        fprintf(formatSpec2, i, ts) % report progress
    end
end
clearvars i ii

%% save dca Results
save(saveName,'dcaResults')

%% plot DCA dCovs colormap
dCovsMat = nan(size(dcaResults.dCovs,1),size(dcaResults.dCovs,2),5); % timeBin x timeBin x dimensions
xAxis = 200:200:4000; xAxis = xAxis-3000;
yAxis = 200:200:4000; yAxis = yAxis-3000;

for i = 1:size(dCovsMat,3) % increment dca dimensions, e.g. 5
    for ii = 1:size(dCovsMat,1)
        for iii = 1:size(dCovsMat,2)
            dCovsMat(ii,iii,i) = dcaResults.dCovs{ii,iii}(i);
        end
    end
    % plot
    if i == 1 % for the first DCA dimension
        cLimits = [0 floor(max(max(dCovsMat(:,:,i))))];
        %cLimits = [0 6];
    end
        fig = figure;     
        surf(xAxis,yAxis,dCovsMat(:,:,i));   
        colorbar
        view(0,90)
        axis tight
        set(gca,'TickDir','out')
        set(gca,'GridLineStyle','none')
        caxis(cLimits)
        pbaspect([1 1 1]); % plot box aspect (ratio)
        %print(fig, strcat('dCovs','Dim#',num2str(i)),'-dpdf'); % print figure as pdf
        clearvars fig
end
clearvars i ii iii

%% plot dCovs vs. trial-shuffled dCovs
clear g
g(1,1)=gramm('x',1:numbBin,'y',[dcaResults.dCovsDiag; cell2mat(dcaResults.dCovsShuf)],'color',[{'dCovs'},c]);
g(1,2)=copy(g(1));

g(1,1).geom_point();
g(1,1).set_title('geom_point()');

g(1,2).stat_smooth();
g(1,2).set_point_options('base_size',3);
g(1,2).set_title('stat_smooth()');

g.set_title('');

fig = figure('Position',[100 100 800 350]);
g.draw();
pbaspect([1 1 1]);

print(fig, strcat('dCovs_',saveName), '-dpdf'); % print figure as pdf

end