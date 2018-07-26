function [features] = TNC_EventQuant(rawEvents,target,distance,pcaStruct)
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
% 

    getDistance = 1;
    disp(' ');
    
    switch lower(target)

        case 'exp' % create an exponential template
            [template] = TNC_CreateEXPTemplate(rawEvents.Sh_x);
            features.exp.template = template;
            features.template = template;    

        case 'frequency' % use a sinc function of multiple bandwidths
            [template] = TNC_CreateSINCTemplate(rawEvents.Sh_x);
            features.frequency.template = template;
            features.template = template;    
            
        case 'gmono' % use a gaussian monopulse as the template (appears quite good for interneuron waveforms)
            [template] = TNC_CreateGMONTemplate(rawEvents.Sh_x);
            features.gmono.template = template;
            features.template = template;    

        case 'pca' % use svd to find the components for the projection
            if isempty(pcaStruct)
                [template] = TNC_CreatePCATemplate(rawEvents);
                features.pca.template = template;
                features.template = template;
            else
                features.template = pcaStruct;
            end
            
        case 'scalar' % do not create a template
            getDistance = 0;
            numSamps = size(rawEvents.Sh_wfs,2);
            for j=1:1:size(rawEvents.Sh_wfs,1)
                tmp(j,1) = min(rawEvents.Sh_wfs(j,:));
                tmp(j,2) = max(rawEvents.Sh_wfs(j,:));
                tmp(j,3) = sum(rawEvents.Sh_wfs(j,:).^2);
                tmp(j,4) = abs(tmp(j,2) - tmp(j,1));
                tmp(j,5) = find(rawEvents.Sh_wfs(j,:)==tmp(j,1),1);
                tmp(j,6) = numel(find(rawEvents.Sh_wfs(j,:)<(tmp(j,1)./2)));
                tmp(j,7) = numel(find(rawEvents.Sh_wfs(j,:)>(tmp(j,2)./2)));
                tmp(j,8) = trapz(abs(rawEvents.Sh_wfs(j,:))); % quintile values
                tmp(j,9) = trapz( (rawEvents.Sh_wfs(j,:))); % quintile values
                tmp(j,10) = mean(rawEvents.Sh_wfs(j,numSamps-10:numSamps)); % quintile values                
            end

            features.template = [];
            
            eval(['features.' lower(target) '.scl.values = tmp;']);

    end

    disp(sprintf('Number of templates generated: %g',size(features.template,2)));
    
    if getDistance

% DISTANCE METRICS ARE Z-SCORED FOR COMPARISON AND CLUSTERING
        switch lower(distance)

            case 'pdist'
                for j=1:1:size(rawEvents.Sh_wfs,1)
                    tmp(j,:) = pdist2(rawEvents.Sh_wfs(j,:),features.template(:,:)','euclidean');
                end

                for p = 1:size(tmp,2)
                    tmp(:,p) = (tmp(:,p)-mean(tmp(:,p))) ./ std(tmp(:,p));
                end
                
                eval(['features.' lower(target) '.pdist.values = tmp;']);

            case 'corr'
                tmp = corr(rawEvents.Sh_wfs',features.template);

                for p = 1:size(tmp,2)
                    tmp(:,p) = (tmp(:,p)-mean(tmp(:,p))) ./ std(tmp(:,p));
                end
                
                eval(['features.' lower(target) '.corr.values = tmp;']);
                
            case 'corrsh'
                for j=1:1:size(rawEvents.Sh_wfs,1)
                    numTemplates = size(features.template,2);
                    for m = 1:numTemplates
                        tmpTmp = xcorr(rawEvents.Sh_wfs(j,:),features.template(:,m),30);
                        tmp(j,m) = max(tmpTmp);
                    end
                end
                
                for p = 1:size(tmp,2)
                    tmp(:,p) = (tmp(:,p)-mean(tmp(:,p))) ./ std(tmp(:,p));
                end
                
                eval(['features.' lower(target) '.corrSh.values = tmp;']);


            case 'dot'
                tmp = rawEvents.Sh_wfs*features.template;

                for p = 1:size(tmp,2)
                    tmp(:,p) = (tmp(:,p)-mean(tmp(:,p))) ./ std(tmp(:,p));
                end
                
                eval(['features.' lower(target) '.dot.values = tmp;']);
                
            case 'ldist'
                numTemplates = size(features.template,2);
                ones(size(features.template));
                for j=1:1:size(rawEvents.Sh_wfs,1)
                    tmpTmp = (rawEvents.Sh_wfs(j,:).*ones(size(features.template)))'-features.template(:,:);
                    tmp(j,:) = sum(tmpTmp,1);
                end

                for p = 1:size(tmp,2)
                    tmp(:,p) = (tmp(:,p)-mean(tmp(:,p))) ./ std(tmp(:,p));
                end
                
                eval(['features.' lower(target) '.ldist.values = tmp;']);
            
            case 'none'
                disp('No distance metrics were generated, only templates.');

        end

    end
    
    disp(' ');
    
end


%_________________________________________________________________________
%
% HELPER FUNCTIONS
%_________________________________________________________________________

function [template] = TNC_CreateEXPTemplate(x)
    
    samples = 0:1:length(x)-1;

    phaseOne = -6.*(1-exp(samples./-2)).*exp(-samples./3);
    phaseTwo = 4.*(1-exp(samples./-7)).*exp(-samples./7);

    phaseOneB = -8.*(1-exp(samples./-1)).*exp(-samples./2);
    phaseTwoB = 8.*(1-exp(samples./-2)).*exp(-samples./2);

    for i=14:length(x)
        if i>15
            fakeSpikeA(1,i) = phaseOne(i-10) + phaseTwo(i-14);
            fakeSpikeB(1,i) = phaseOneB(i-10) + phaseTwoB(i-14);
        else
            fakeSpikeA(1,i) = phaseOne(i-10);        
            fakeSpikeB(1,i) = phaseOneB(i-10);
        end
    end
    
    template(:,1) = -(fakeSpikeA ./ min(fakeSpikeA))';
    template(:,2) = -(fakeSpikeB ./ min(fakeSpikeB))';
    
end

function [template] = TNC_CreateSINCTemplate(x)
    template(:,1) = -sinc(2.*pi.*(x-5).*0.01);
    template(:,2) = -sinc(2.*pi.*(x-5).*0.02);
    template(:,3) = -sinc(2.*pi.*(x-5).*0.04);
    template(:,4) = -sinc(2.*pi.*(x-5).*0.2);
end

function [template] = TNC_CreateGMONTemplate(x)
    template(:,1) = gmonopuls((x-5)/300000,10000);
    template(:,2) = gmonopuls((x-5)/300000,20000);
    template(:,3) = gmonopuls((x-5)/300000,30000);
end

function [template] = TNC_CreatePCATemplate(rawEvents)

    numSpks = size(rawEvents.Sh_wfs(:,:),1);
    % only use 7500 spikes to generate PCA
    if numSpks<7500
        [pca.u,pca.s,pca.v] = svd(detrend(rawEvents.Sh_wfs(:,:)','constant')');
    else
        [pca.u,pca.s,pca.v] = svd(detrend(rawEvents.Sh_wfs(1:7500,:)','constant')');
    end
    numComponents = 5;
    template(:,1:numComponents) = pca.v(:,1:numComponents);
    d = diag(pca.s);
    
    disp(sprintf('Fraction of variance captured in first %g components: %g', numComponents, sum(d(1:numComponents,1))/sum(d) ));

end










