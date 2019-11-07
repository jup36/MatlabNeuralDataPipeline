% This code runs a simulation to compare PCA and LEM.

%% get a neuron-by-observation matrix
% define tuning of each 21 neurons
mu(1,:) = 0:2*pi/20:2*pi; % tuning curve center per neuron
mu(2,mu(1,:)<pi) = mu(1,mu(1,:)<pi)+2*pi; % tuning curve shifted by 2pi
mu(2,mu(1,:)>=pi) = mu(1,mu(1,:)>=pi)-2*pi; % tuning curve shifted by 2pi

sigma = [pi/3, pi/6, pi/9, pi/12, pi/15, pi/20]; % wide & narrow tuning
lambda = 15; % to draw from poisson
x = 0:2*pi/1000:2*pi-2*pi/1000; % 0 to 2pi
stim = 2*pi*rand(500,1)'; % random sample  
%stim = repmat(stim,[1,10]); 

figSaveDir = '/Users/parkj/Dropbox (HHMI)/Presentation/labMeeting_101019/Figure';
[cmap]=cbrewer('div', 'Spectral', 6);

for i=1:length(sigma) % increment sigmas
    % determine lower and upper bound for neuronal tuning 
    for ii = 1:size(mu,2) % increment neurons     
        % define tuning to sample vs. radians per neuron
        p1 = normpdf(x,mu(1,ii),sigma(i)); % centered tuning curve from stim
        p2 = normpdf(x,mu(2,ii),sigma(i)); % tuning from x (sample) curve shifted by 2pi
        rez(i).pdf(ii,:) = max([p1;p2],[],1); % tuning curves per neuron
        
        % define tuning to stim vs. radians per neuron
        p11 = normpdf(sort(stim),mu(1,ii),sigma(i)); % centered tuning curve from stim
        p22 = normpdf(sort(stim),mu(2,ii),sigma(i)); % tuning from stim curve shifted by 2pi
        rez(i).pp(ii,:) = max([p11;p22],[],1); % tuning curves per neuron
        rez(i).scaledTCstim(ii,:) = rez(i).pp(ii,:)*(lambda/max(rez(i).pp(ii,:))); % scale the tuning curves for the max rate/count to be at lambda
        rez(i).poissCntStim(ii,:) = poissrnd(rez(i).scaledTCstim(ii,:)); % draw poisson random counts from tuning curves     
    end
    figure; plot(rez(i).pdf','Color',cmap(i,:)); title(strcat('TuningCurves',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none'); 
    %print(fullfile(figSaveDir,sprintf('tuningCurvePerSigma#%d',i)),'-dpdf')
    close all; 
    
    %% Run PCA
    imagesc(cov(rez(i).scaledTCstim))
    pbaspect([1 1 1])
    title(strcat('Covariance matrix',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none')
    print(fullfile(figSaveDir,sprintf('covMatPerSigma#%d',i)),'-dpdf')
    
    rez(i).pcDim = pca(rez(i).scaledTCstim'); 
    rez(i).pcProj = rez(i).poissCntStim'*rez(i).pcDim(:,1:2); 
    
    rez(i).pcDimT = pca(rez(i).scaledTCstim); 
    rez(i).pcProjT = rez(i).poissCntStim*rez(i).pcDimT(:,1:2);  
    
    %% Run LEM with data points as nodes and data-point-pairs as edges
    [m,n] = size(rez(i).scaledTCstim);
    Dt = sqrt(sum((reshape(rez(i).scaledTCstim, [n,1,m])-reshape(rez(i).scaledTCstim, [1,n,m])).^2,3)); % pairwise distance bewteen time points
    upTriuDt = triu(Dt); % take the upper triangle 
    sortDt = sort(upTriuDt(upTriuDt~=0)); 
    nearestCut = sortDt(round(length(sortDt)*(10/100)),1); % distance threshold to define nearest neighbors
    adjDt = Dt<nearestCut; % the adjacency matrix
    imagesc(adjDt)
    title(strcat('Adjacency matrix',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none')
    pbaspect([1,1,1])
    print(fullfile(figSaveDir,sprintf('AdjMatPerSigmaEachDim#%d',i)),'-dpdf')
    
    deg = diag(sum(adjDt, 2)); % Compute the degree matrix
    L = deg - adjDt;    % Compute the laplacian matrix in terms of the degree and adjacency matrices
    imagesc(L)
    title(strcat('Laplacian matrix',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none')
    pbaspect([1,1,1])
    colormap('jet')
    caxis([-50 30])
    print(fullfile(figSaveDir,sprintf('LMatPerSigmaEachDim#%d',i)),'-dpdf')
    
    
    [rez(i).lemDim, rez(i).lemEigVals] = eig(L);    % Compute the eigenvalues/vectors of the laplacian matrix
    % lemdim has the topology among data points: drop the 1st vector and use 2nd and 3rd vectors. E.g. plot(rez.lemDim(:,2),lemDim(:,3)). 
    
    % plot to compare LEM vs. PCA-TimeDim vs. PCA-NeuronDim
    figure; 
    subplot(1,3,1)
    plot(rez(i).lemDim(:,2), rez(i).lemDim(:,3),'.','Color',cmap(i,:))       
    title(strcat('LEM',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none')
    pbaspect([1,1,1])
    
    subplot(1,3,2)
    plot(rez(i).pcDimT(:,1), rez(i).pcDimT(:,2),'.','Color',cmap(i,:))
    title(strcat('PCA-TimeDim',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none')
    pbaspect([1,1,1])
    
    subplot(1,3,3)
    plot(rez(i).pcProj(:,1), rez(i).pcProj(:,2),'.','Color',cmap(i,:))
    title(strcat('PCA-NeuronDim',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none')
    pbaspect([1,1,1])
    %print(fullfile(figSaveDir,sprintf('lemPcaPerSigma#%d',i)),'-dpdf')
    close all;

    % plot to compare LEM vs. PCA-TimeDim vs. PCA-NeuronDim
    figure; 
    subplot(1,2,1)
    hold on; 
    plot(rez(i).lemDim(:,2),'Color',cmap(i,:)) 
    plot(rez(i).lemDim(:,3),'Color',cmap(i,:))       
    hold off; 
    title(strcat('LEM',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none')
    pbaspect([1,1,1])
    
    subplot(1,2,2)
    hold on; 
    plot(rez(i).pcDimT(:,1),'Color',cmap(i,:))
    plot(rez(i).pcDimT(:,2),'Color',cmap(i,:))
    hold off; 
    title(strcat('PCA-TimeDim',sprintf('_w/Sigma_%.3g',sigma(i))),'Interpreter', 'none')
    pbaspect([1,1,1])
    print(fullfile(figSaveDir,sprintf('lemPcaPerSigmaEachDim#%d',i)),'-dpdf')
    
end
   


