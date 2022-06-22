function [features] = TNC_EventQuantSE(rawEvents,template)
% FUNCTION DETAILS: For a set of discrete events measure the "distance" between each event and a "target" waveform. Distances and targets are chosen to maximize the "clusterability" of the waveforms for subsequent steps. In many cases multiple targets and distance metrics should be used. See the help for details on supported targets and distance metrics.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% INPUTS:
% targets: 'template', 'frequency','pca', 'none'
% distance: 'pdist', 'corr', 'dot', 'ldist', 'none'
% 
% OUTPUTS:
% appear in the structure as .__target__.__distance__.values

                        
% NOTE: implement the theta + rho formalism in EventQuantSE

%% Set parameters
    numComponents = 3;          

%% Get scalar quantification
    
    numSamps = size(rawEvents.wfs(1).values,2);
    numEvents = numel(rawEvents.wfs);
    
    sclFeats= zeros(numEvents,10);
    theta   = zeros(numEvents,1);
    rho     = zeros(numEvents,1);
    
    for j=1:1:numEvents
        
        posInds = find(rawEvents.wfs(j).values(1,:)> 60);
        negInds = find(rawEvents.wfs(j).values(1,:)<-60);

        forPCA(j,:) = rawEvents.wfs(j).values(1,:);
        
        sclFeats(j,1) = min(rawEvents.wfs(j).values(1,:));                             % min value
        sclFeats(j,2) = max(rawEvents.wfs(j).values(1,:));                             % max value
        sclFeats(j,3) = trapz(rawEvents.wfs(j).values(1,posInds)) - trapz(rawEvents.wfs(j).values(1,negInds)); % peak-valley area
        sclFeats(j,4) = sclFeats(j,2) - sclFeats(j,1);                       % peak-valley amplitude
        sclFeats(j,5) = find(rawEvents.wfs(j).values(1,:)==sclFeats(j,1),1);           % time of peak
        sclFeats(j,6) = numel(negInds);                                      % valley width
        sclFeats(j,7) = numel(posInds);                                      % peak width
        sclFeats(j,8) = trapz(abs(rawEvents.wfs(j).values(1,:)));                      % energy
        sclFeats(j,9) = sum(rawEvents.wfs(j).values(1,:).^2);                          % nonlinear energy
        sclFeats(j,10) = median(abs(diff(rawEvents.wfs(j).values(1,:))));              % median derivative
        
        theta(j,1)  = atan2(sclFeats(j,1),sclFeats(j,2));
        rho(j,1)    = sqrt(sclFeats(j,1).^2 + sclFeats(j,2).^2);
        
        ts(j,1)     = rawEvents.inds(j);
        
    end

%% Calculate the components of feature space
    
    if isempty(template)
        [pca.u,pca.s,pca.v] = svd(sclFeats(:,[1 2 4 6 7 10]));
        features.template(:,1:numComponents) = pca.v(:,1:numComponents);
        d = diag(pca.s);    
        %disp(sprintf('     EVENT QUANT | Variance in first %g components of feature space: %g', numComponents, sum(d(1:numComponents,1))./sum(d) ));
    else
        features.template = template;        
    end

%% Project scalar quantities into lower dimensional feature space

    tmp = sclFeats(:,[1 2 4 6 7 10])*features.template;

    
%% Calculate PCA on the waveforms
    
    [mappedA, mapping] = compute_mapping(forPCA', 'PCA', 3);
    pcWaveform = forPCA*mappedA;

%% Update the output structure

    features.params = [ts sclFeats pcWaveform tmp theta rho];
    features.paramNames = {'ts' , 'min', 'max', 'pevaA', 'pevaV', 'pTime', 'widthV', 'widthP', 'energy', 'n energy', 'medDiff', 'PC1', 'PC2', 'PC3', 'fPC1', 'fPC2', 'fPC3', 'theta19', 'rho19'};
    
