function [features] = TNC_EventQuantME(data, wfsStruct, pcaStruct, spikeInds, sampleRate,...
                                       winL, shankSize)
% FUNCTION DETAILS: 
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% INPUTS:


% FOR MULTISHANK PROBES JUST SCALAR CALCULATIONS SHOULD BE SUFFICIENT
    numSamps = size(wfsStruct(1).values,2);
    eventInds = spikeInds';
    timestamps = round(spikeInds(:,1)./sampleRate);
    for j=1:1:size(wfsStruct,2)
        features.minV(:,j)   = data(:,eventInds(j));
        features.maxV(:,j)   = max(wfsStruct(j).values(:,winL-10:winL+10),[],2);
        features.energy(:,j) = trapz(abs(wfsStruct(j).values(:,winL-10:winL+20)),2);
        features.peVa(:,j)   = features.maxV(:,j) - features.minV(:,j);
        features.minT(:,j)   = eventInds(:,j)./sampleRate;
    end
 
    if isempty(pcaStruct)
        [pca] = TNC_CreatePCATemplate(features);
    else
        pca = pcaStruct;
    end

    features.pca = pca;
    
% PROJECT IN TO LOW DIMENSIONS

    for j=1:1:size(wfsStruct,2)
        tmpProj(j,:) = features.minV(:,j)' * pca.template;
    end

    for p = 1:size(tmpProj,2)
        [h,pVal,ksstat,cv] = kstest(tmpProj(:,p),[],0.01,'unequal');
        mmStat(p) = ksstat;
    end

    features.lowd = tmpProj;
%   disp(['size(timestamps)=' num2str(size(timestamps)) ' size(features.lowd)=' num2str(size(features.lowd)) ...
%         ' size(features.energy)=' num2str(size(features.energy)) ' size(features.minV)=' num2str(size(features.minV))]);
 
    features.params = [timestamps, features.lowd, features.energy' features.minV' features.minT'];
    features.paramNames = cell(1,1+4*shankSize); 
    features.paramNames{1} = 'ts';
    for i=1:shankSize
        features.paramNames{1+i            } = strcat('cPC', num2str(i));
        features.paramNames{1+i+shankSize  } = strcat('eng', num2str(i));
        features.paramNames{1+i+shankSize*2} = strcat('min', num2str(i)); 
        features.paramNames{1+i+shankSize*3} = strcat('minT', num2str(i)); 
    end
end

%_________________________________________________________________________
%
% HELPER FUNCTION
%_________________________________________________________________________


    function [pca] = TNC_CreatePCATemplate(tmp)

        numSpks = size(tmp.minV,2);
        numChan = size(tmp.minV,1);

        [pca.u,pca.s,pca.v] = svd(tmp.minV');

        numComponents = numChan;
        template(:,1:numComponents) = pca.v(:,1:numComponents);
        d = diag(pca.s);

    %   disp(['size(d)=' num2str(size(d)) ' sum(d(1:3,1))=' num2str(sum(d(1:3,1))) ' sum(d)=' num2str(sum(d,1))]); 
        disp(sprintf('Fraction of variance captured in first %g components: %g', 3, sum(d(1:3,1))/sum(d(:,1)) ));

        pca.template = template;
        pca.numSpks = numSpks;
        pca.numChan = numChan;

    end
