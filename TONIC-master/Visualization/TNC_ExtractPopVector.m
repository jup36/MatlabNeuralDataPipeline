function [popVec] = TNC_ExtractPopVector(phys,causal)

    currParams.smthParams.rise  = 25;
    currParams.smthParams.decay = 100;
    if causal
        [currParams.filter.kernel]  = TNC_CreateCausalKernel(currParams.smthParams.rise,currParams.smthParams.decay,1);
    else
        [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
    end
    popVec.currParams           = currParams;
    
    totTime = ceil(phys.maxTime./30);
    popVec.totTime = totTime;
    
% Create smoothed firing rates for each unit

    corrMat = zeros(totTime,phys.numUnits);

    for i=1:phys.numUnits
        delta = zeros(1,totTime);
        delta(round(phys.unit(i).ts./30)) = 1;
        popVec.unit(i).smth = conv(delta,currParams.filter.kernel,'same');
        if mean(popVec.unit(i).smth) > 0
            corrMat(:,i) = popVec.unit(i).smth;
        end
    end
    
    figure(2); 
    subplot(121); imagesc((corrMat));
    subplot(122); imagesc(corr(corrMat));
    drawnow;
    
% Compute the first three PCs

    [u,s,v] = svd(corrMat,'econ');

    popVec.pca.component(:,1:4) = v(:,1:4);
    d = diag(s);
    
    disp(sprintf('Fraction of variance captured in first %g components: %g', 3, sum(d(1:4,1))/sum(d) ));
  
% Project spiking into PC space at 1kHz resolution

    % estimate percent complete
    tenP = round(totTime./10);
    
    for i=1:totTime

        popVec.proj(1,i) = dot(corrMat(i,:)',popVec.pca.component(:,1));
        popVec.proj(2,i) = dot(corrMat(i,:)',popVec.pca.component(:,2));
        popVec.proj(3,i) = dot(corrMat(i,:)',popVec.pca.component(:,3));
        popVec.proj(4,i) = dot(corrMat(i,:)',popVec.pca.component(:,4));
    
        if rem(i,tenP)==0
            disp([num2str(round(i/totTime.*100)) '% complete.']);
        end
        
    end
