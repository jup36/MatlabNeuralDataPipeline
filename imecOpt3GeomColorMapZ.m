function [figHandle] = imecOpt3GeomColorMapZ( imecSites, Z, colorScheme, drawAllSites, colorAxis, probeDepth, varargin )
%This function generates a 'functional' probe map. Each probe site is drawn
% as a patch with its color coded to represent the quantity Z (a normalized quantity). 
% If drawAllSites boolean is true, all imecSites of the probe are patched,
% whereas false, imecSites within the relevant range of the probe are patched.
% Input: 
%   imecSites: imec sites
%   Z: normalized quantities across the imec sites

colorEdge = linspace(-max(abs(colorAxis)),max(abs(colorAxis)),100); % get a zero-centered color edges
[~, ~, colorIdx] = histcounts(Z, colorEdge);    % locate each Z score on the colorEdge
colorIdx(Z>max(colorEdge)) = length(colorEdge)-1; % this prevents assigning 0 to Z values greater than the upper limit of the colorEdge

[colorMap] = TNC_CreateRBColormap(length(colorEdge)+1, colorScheme); % generate a colormap

% imec option3 site location in micrometers (x and y)
geometry = zeros(384, 2); % col1: X, col2: Y
viHalf = 0:(384/2-1); % 192 rows
viHalf = -flip(viHalf,2); % flip the dimension and put negative to represent depth on y-axis
geometry(1:2:end,2) = viHalf * 20; % increment each row by 20 um
geometry(2:2:end,2) = geometry(1:2:end,2);
geometry(1:4:end,1) = 4;  %16; 
geometry(2:4:end,1) = 12; %48; 
geometry(3:4:end,1) = 0;  %0 
geometry(4:4:end,1) = 8;  %20;

xBounds = [0; 4; 4; 0];   % X bounds for each site
yBounds = [0; 0; 40; 40]; % Y bounds for each site

geometry(:,2) = geometry(:,2)-(probeDepth-(abs(geometry(1,2))))+40; % adjust to make the y-coordinates of the bottom site to be the probeDepth

% draw background imecSites first
if drawAllSites % draw all 384 imecSites
    xRep = repmat(geometry(:,1)', 4, 1) + repmat(xBounds,1,length(geometry)); % X coordinates
    yRep = repmat(geometry(:,2)', 4, 1) + repmat(yBounds,1,length(geometry)); % Y coordinates
    figHandle = figure('Position',[100 100 400 1200]); % figure('Position',[left bottom width height])
    patch(xRep,yRep,[0.7,0.7,0.7]) % patch all imecSites in grey
    axis tight
else
    xRep = repmat(geometry(min(imecSites):end,1)', 4, 1) + repmat(xBounds,1,length(geometry(min(imecSites):end,1))); % X coordinates
    yRep = repmat(geometry(min(imecSites):end,2)', 4, 1) + repmat(yBounds,1,length(geometry(min(imecSites):end,1))); % Y coordinates
    figHandle = figure('Position',[100 100 400 1200]); % figure('Position',[left bottom width height])
    patch(xRep,yRep,[0.7,0.7,0.7]) % patch imecSites in the range in grey
    axis tight
end

hold on;
xSitesRep = repmat(geometry(imecSites,1)',4,1) + repmat(xBounds,1,length(imecSites)); % x coordinates for imecSites
ySitesRep = repmat(geometry(imecSites,2)',4,1) + repmat(yBounds,1,length(imecSites)); % y coordinates for imecSites


for i = 1:size(imecSites,1)
    patch(xSitesRep(:,i),ySitesRep(:,i),colorMap(colorIdx(i)+1,:)) % patch all imecSites in grey
end
hold off;
caxis([-max(abs(colorAxis)) max(abs(colorAxis))]);
colormap(colorMap);
colorbar

end

