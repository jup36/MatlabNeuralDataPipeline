function [mappedData] = TNC_POP_DimReduce(PopData,sessNum,smthTau,method,nDims)

%% STANDARD ANALYSIS PARAMETERS

    % Define the smoothing to apply to the timeseries data
    currParams.smthParams.decay    = smthTau;
    currParams.filter.causal       = 0;

    [kernel]  = TNC_CreateGaussian(smthTau.*15,smthTau,smthTau.*30,1);
    if currParams.filter.causal
        kernel(1:smthTau.*15) = 0;
    end
    
    [mapName] = TNC_CreateRBColormap(1024,'bo');

    disp(' ');
    disp(' ');
    disp('________________________________');
    disp(' ');
    disp('Initialized analysis parameters.');
    disp('________________________________');
    disp(' ');
    disp(' ');
    
%% COMPILE DATA IN TO A MATRIX SUITABLE FOR DIM REDUCTION
    
    clear dataMatrix

    for index = 1:numUnits        
        dataMatrix(:,index) = PopData.session(sessNum).unit(index).response.psthZ';
        shankLabels(1,index) = PopData.session(sessNum).unit(index).sh;
    end    
    
%% CALCULATE THE REDUCED DIMENSIONS

    [mappedA, mapping] = compute_mapping(dataMatrix, method, nDims);

%% PLOT THE PROJECTION

    figure(10); clf; scatter3(mappedA(:,1),mappedA(:,2),mappedA(:,3),10,[-2000:5000]); hold on; plot3(mappedA(2000,1),mappedA(2000,2),mappedA(2000,3),'ko','MarkerSize',10,'LineWidth',4);    

%% DUMP THE OUTPUT

    mappedData = mappedA;

