function [unit] = TNC_EventCluster(featureData,clustNumber,cMethod,display)
% FUNCTION DETAILS: Performs the clustering (elsewhere called 'classification') step where individual events are given integer ids as members of particular clusters. Further a confidence value is also returned that provides a metric for the distance of each event from the center of the cluster of which it is a member.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% INPUTS:
% cMethod: 'gmix', 'km','hiearch', 'thresh', 'none'
% eMethod: 'euclidean', 'cosine', 'none'
% 
% OUTPUTS:
% appear in the structure as .clusters.ids
% appear in the structure as .clusters.distances
%

switch cMethod
    
    case 'gmix'
        options     = statset('Display','final');
        gm          = gmdistribution.fit(featureData,clustNumber,'Options',options);
        clustIds    = cluster(gm,featureData);
        unit.clustIds = clustIds;
        
    case 'km'
        if numel(clustNumber)~=1
            disp('Kmeans clustering requires the user to specify the number of clusters');
        else
            [clustIds, clustCenters] = kmeans(featureData,clustNumber,'Replicates',50,'Start','cluster');
            unit.clustIds = clustIds;
            unit.clustCenters = clustCenters;
        end
         
    case 'hierarch'
        distances   = pdist(featureData); % rows are observations, columns are variables
        links       = linkage(distances);
        clustIds    = cluster(links,'maxclust',clustNumber);
        unit.clustIds = clustIds;

    case 'thresh' 
        % a simple implementation of clustering for clear cases that uses
        % the properties of the histogram of projection values.
        
    case 'contour'
        % recently i have been thinking about how to implement this
        % algorithm. The implementation that I used for the grid cell
        % analysis seems quite relevant and is fairly simple to implement.
        % The algorithm is (1) Get contours of an image plot, (2) use a
        % heuristic to decide on the threshold level of the contour, (3)
        % find the nearest contour, (4) calculate the center of the given
        % contour, (5) for all chosen centers calculate the pairwise
        % distance between a given event and the center points, (6) define
        % a heuristic that chooses cluster membership based upon
        % distribution of distances, (7) for exach spike return the
        % distance metric and cluster id.
        
    case 'none'
        
end

if display==1
    % display the clustered output
    figure(200); clf;
    set(gcf, 'color', [0 0 0]);
    featureDims = size(featureData,2);

    for k=1:featureDims.^2
        set(gca, 'color', [0 0 0]);
        subplot(featureDims,featureDims,k);

        [xI,yI] = ind2sub([featureDims featureDims],k);
        if xI~=yI
            scatter(featureData(:,xI),featureData(:,yI),2,clustIds,'filled'); hold on;
            scatter(clustCenters(:,xI),clustCenters(:,yI),84,'k');
            whitebg('black');
            axis tight; axis off;
        else
            axis off;
        end
    end
end
