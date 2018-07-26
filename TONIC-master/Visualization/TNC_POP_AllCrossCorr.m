function [] = TNC_POP_AllCrossCorr(PopData,sessNum,unitArray,maxLag,figNum)

%% STANDARD ANALYSIS PARAMETERS
figure(figNum); clf;

kernel  = TNC_CreateGaussian(2.*15,2,2.*30,1);

%% Calculate all pairwise cross correlations

    if numel(unitArray)>1
        noUnitList = 1;
        numUnits = numel(unitArray);
    else
        noUnitList = 0;
        numUnits = numel(PopData.session(sessNum).unit);
    end
    
    for index=1:numUnits
        for jindex=1:index
            
            iStamps = PopData.session(sessNum).unit(index).ts;
            delta1 = zeros(1,ceil(iStamps(numel(iStamps))));
            delta1(iStamps) = 1;
            
            distance = (PopData.session(sessNum).unit(index).sh - PopData.session(sessNum).unit(jindex).sh) ./ 7;
            
            if index==jindex

                tmpCorr = xcorr(delta1,delta1,maxLag,'unbiased');
                diffCorrMat(jindex,index)   = 0;
                distMat(jindex,index)       = 0;
                corrMat(jindex,index)       = 0;
                
            else

                validStamps = find(PopData.session(sessNum).unit(jindex).ts <= iStamps(numel(iStamps)));            
                delta2 = zeros(1,ceil(iStamps(numel(iStamps))));
                delta2(PopData.session(sessNum).unit(jindex).ts(validStamps)) = 1;

%                 smth1 = zeros(size(delta1));
%                 smth2 = zeros(size(delta1));
%                 
%                 smth1 = conv(delta1,kernel,'same');
%                 smth2 = conv(delta2,kernel,'same');

                tmpCorr = xcorr(delta1,delta2,maxLag,'unbiased');
                diffCorrMat(jindex,index)   = ( trapz(tmpCorr(maxLag:2*maxLag)) - trapz(tmpCorr(1:maxLag)) ) ./ ( trapz(tmpCorr(maxLag:2*maxLag)) + trapz(tmpCorr(1:maxLag)) );
                distMat(jindex,index)       = distance;
                corrMat(jindex,index)       = tmpCorr(maxLag+1); 
                
            end

            figure(figNum);
            subplot(numUnits,numUnits,sub2ind([numUnits,numUnits],index,jindex));
            plot(-maxLag:maxLag,tmpCorr,'Color',[distance 0.5 1-distance]); 
            hold on; plot([0 0],[0 max(tmpCorr)],'k');
            axis tight; axis off;
            
        end
    end
    
    figure(figNum+1); clf;
    subplot(211);
    [mapName] = TNC_CreateRBColormap(1024,'bo');
    imagesc(diffCorrMat,[-0.7 0.7]);
    colormap(mapName);
    subplot(212);
    plot(distMat(1:numel(distMat)).*1400,diffCorrMat(1:numel(distMat)),'ko'); axis([-100 1300 -1 1]);
    
%% Display as a function of distance