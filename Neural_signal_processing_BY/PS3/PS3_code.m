function PS3_3
cd('C:\Users\jup36\Dropbox\NSP\PS3')
load('ps3_simdata.mat','-mat');
[NumData NumClass]=size(trial);
for classIX=1:NumClass
    for dataIX=1:NumData
        dataArr(classIX,dataIX,:)=trial(dataIX,classIX).x;
    end
end
NumFea=size(dataArr,3);
%% %%%%%%%%%%%%%%%%%%%%%%%%% Part (b) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumModel=3;
% For model (1)
modParam{1}.mean=squeeze(mean(dataArr,2));  % squeeze removes all singleton dimensions
% Remove the mean of each class
dataArrRemMean=reshape(repmat(modParam{1}.mean,1,NumData),[NumClass NumFea NumData]);   % repmat replicate and tile array
dataArrRemMean=dataArr-permute(dataArrRemMean,[1 3 2]);     %  permute rearrange dimensions of N-D array
% Shared covariance matrix
modParam{1}.cov{1}=cov(reshape(dataArrRemMean,[],size(dataArr,3)));     % shared covariance
% For model (2)
modParam{2}.mean=squeeze(mean(dataArr,2));  % squeeze removes all singleton dimensions
for classIX=1:NumClass
    modParam{2}.cov{classIX}=cov(squeeze(dataArr(classIX,:,:)));
end
% For model (3)
modParam{3}.mean=squeeze(mean(dataArr,2));
for modelIX=1:NumModel      %
    %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MarkerPat={'rx','g+','bo'};
    figure(modelIX);
    for classIX=1:NumClass
        plot(squeeze(dataArr(classIX,:,1)),squeeze(dataArr(classIX,:,2)),...
            MarkerPat{classIX},'LineWidth',2,'MarkerSize',8);
        hold on;
    end
    axis([0 20 0 20]);
    xlabel('x_1');
    ylabel('x_2');
    %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MarkerCol={'r','g','b'};
    for classIX=1:NumClass
        plot(modParam{modelIX}.mean(classIX,1),modParam{modelIX}.mean(classIX,2),...
            'o','MarkerEdgeColor',MarkerCol{classIX},'MarkerFaceColor',...
            MarkerCol{classIX},'MarkerSize',10)
        hold on;
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (d) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Skip this part if using model (3)
    if modelIX<3
        NumCov=length(modParam{modelIX}.cov);   % shared covariance
        [X,Y] = meshgrid(0:.1:20,0:.1:20);
        feaVec=cat(3,X,Y);
        feaVec=reshape(feaVec,[],size(feaVec,3));
        for classIX=1:NumClass
            currMean=modParam{modelIX}.mean(classIX,:);
            covIX=min(classIX,NumCov);
            currCov=modParam{modelIX}.cov{covIX};
            % For each f=(x,y) calculate:
            % z=exp(-(x-f)*inv(cov)*(x f)/2)/sqrt(det(cov))
            Z=sum(((feaVec-repmat(currMean,size(feaVec,1),1))*inv(currCov)).*(feaVec-repmat(currMean,size(feaVec,1),1)),2);
            Z=exp(-Z/2)/sqrt(det(currCov));
            Z=reshape(Z,size(X));
            isoThr=0.02;
            contour(X,Y,Z,isoThr,MarkerCol{classIX},'LineWidth',2);
            hold on;
        end
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (e) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate dense samples
    [X,Y] = meshgrid(0:.1:20,0:.1:20);
    feaVec=cat(3,X,Y);
    feaVec=reshape(feaVec,[],size(feaVec,3));
    NumGrid=size(feaVec,1);
    switch modelIX
        % If using model (1)
        case {1}
            for classIX=1:NumClass
                tmp=feaVec-repmat(modParam{modelIX}.mean(classIX,:),NumGrid,1);
                logP(:,classIX)=sum(tmp*inv(modParam{modelIX}.cov{1}).*tmp,2);  % similar to the Z score
            end
            [minVal classLabel]=min(logP,[],2);
            % If using model (2)
        case {2}
            for classIX=1:NumClass
                tmp=feaVec-repmat(modParam{modelIX}.mean(classIX,:),NumGrid,1);
                logP(:,classIX)=sum(tmp*inv(modParam{modelIX}.cov{classIX}).*...
                    tmp,2)+log(det(modParam{modelIX}.cov{classIX}));
            end
            [minVal classLabel]=min(logP,[],2);
            % If using model (3)
        case {3}
            for classIX=1:NumClass
                currLambda=modParam{modelIX}.mean(classIX,:);
                logP(:,classIX)=-feaVec*log(currLambda')+sum(currLambda);
            end
            [minVal classLabel]=min(logP,[],2);
    end
    h=gscatter(X(:),Y(:),classLabel,'rgb','.',[],'off');
    set(h, 'Markersize', 1);
    hold on;
end
return;