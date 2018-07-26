%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

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
    numSamps = size(wfsStruct(1).values,2); % values    = matrix[K(num.compinents),length_of_waveform]
    eventInds = spikeInds';                 % eventInds = matrix[K(num.compinents),num_spikes] 
%   disp(['size(eventInds)=' num2str(size(eventInds)) ' size(wfsStruct,2)=' num2str(size(wfsStruct,2))]);

    % Compute features
    disp(['size(eventInds)=' num2str(size(eventInds))]);
    num_comp   = size(eventInds,1);
    num_spikes = size(wfsStruct,2);
    if num_comp > 0 && num_spikes> 0
        timestamps = round(median(spikeInds,2)./sampleRate);
        for k=1:num_comp                          
            X = [];
            for j=1:1:num_spikes                 
                features.minV(k,j)   = data(k,eventInds(k,j));
                features.maxV(k,j)   = max(wfsStruct(j).values(k,winL-10:winL+10),[],2);
                features.energy(k,j) = trapz(abs(wfsStruct(j).values(k,winL-10:winL+20)),2);
                features.peVa(k,j)   = features.maxV(k,j) - features.minV(k,j);
                features.minT(k,j)   = eventInds(k,j)./sampleRate;
                X = [X ; wfsStruct(j).values(k,:)];    
            end 
            [C,S,l]=princomp(X);
%       disp(['size(X)=' num2str(size(X)) ' size(C)=' num2str(size(C)) ' size(S)=' num2str(size(S)) ' size(l)=' num2str(size(l))]);
            for j=1:1:num_spikes
                features.pc1(k,j) = S(j,1);
                features.pc2(k,j) = S(j,2);
                features.pc3(k,j) = S(j,3); 
            end
        end
        if isempty(pcaStruct)
            [pca] = TNC_CreatePCATemplate(features);
        else
            pca = pcaStruct;
        end
        features.pca = pca;

        % PROJECT IN TO LOW DIMENSIONS
        %   disp(['size(features.minV)=' num2str(size(features.minV)) ' size(pca.template)=' num2str(size(pca.template))]);
        try 
            for j=1:1:num_spikes
                tmpProj(j,:) = features.minV(:,j)' * pca.template;
            end

            for p = 1:size(tmpProj,2)
                [h,pVal,ksstat,cv] = kstest(tmpProj(:,p),[],0.01,'unequal');
                mmStat(p) = ksstat;
            end
            features.lowd = tmpProj;
        catch
            features.lowd = [];
        end

        %   disp(['size(timestamps)=' num2str(size(timestamps)) ' size(features.lowd)=' num2str(size(features.lowd)) ...
    %         ' size(features.energy)=' num2str(size(features.energy)) ' size(features.minV)=' num2str(size(features.minV))]);

    else
        features.minV = []; features.maxV = []; features.energy = [];
        features.peVa = []; features.minT = []; 
        features.pc1  = []; features.pc2  = []; features.pc3    = [];
        features.pca  = []; features.lowd = [];
    end

    features.params = [timestamps, features.lowd, features.energy' ...
                       features.minV' features.minT' features.pc1' ...
                       features.pc2' features.pc3'];

    features.paramNames = cell(1,1+7*shankSize); 
    features.paramNames{1} = 'ts';
    for i=1:shankSize
        features.paramNames{1+i            } = strcat('cPC', num2str(i));
        features.paramNames{1+i+shankSize  } = strcat('eng', num2str(i));
        features.paramNames{1+i+shankSize*2} = strcat('min', num2str(i)); 
        features.paramNames{1+i+shankSize*3} = strcat('minT',num2str(i)); 
        features.paramNames{1+i+shankSize*4} = strcat('pc1', num2str(i));
        features.paramNames{1+i+shankSize*5} = strcat('pc2', num2str(i));
        features.paramNames{1+i+shankSize*6} = strcat('pc3', num2str(i));
    end
end

%_________________________________________________________________________
%
% HELPER FUNCTION
%_________________________________________________________________________

    function [pca] = TNC_CreatePCATemplate(tmp)

        numSpks = size(tmp.minV,2);
        numChan = size(tmp.minV,1);

        disp(['size(tmp.minV)=' num2str(size(tmp.minV))]);

        if numSpks > 0 && numChan > 0
            try 
                [pca.u,pca.s,pca.v] = svd(tmp.minV');

                numComponents = numChan;
                template(:,1:numComponents) = pca.v(:,1:numComponents);
                d = diag(pca.s);

            %   disp(['size(d)=' num2str(size(d)) ' sum(d(1:3,1))=' num2str(sum(d(1:3,1))) ' sum(d)=' num2str(sum(d,1))]); 
                disp(sprintf('Fraction of variance captured in first %g components: %g', 3, sum(d(1:3,1))/sum(d(:,1)) ));

                pca.template = template;
            catch
                pca.template = [];
            end
        else
            pca.template = []; 
        end
        pca.numSpks = numSpks;
        pca.numChan = numChan;
    end
