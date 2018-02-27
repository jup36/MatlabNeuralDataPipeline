function [figHandle] = imecOpt3GeomColorMap( sites, Z, drawAllSites )
%This function generates a 'functional' probe map. Each probe site is drawn
% as a patch with its color coded to represent the quantity 'Z'. 
% If drawAllSites boolean is true, all sites of the probe are patched,
% whereas false, sites within the relevant range of the probe are patched.
% Ensure that the order of data in sites and Z match!

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

%Xs = unique(geometry(:,1)); % X coordinates
%Ys = unique(geometry(:,2)); % Y coordinates

xBounds = [0; 4; 4; 0];   % X bounds for each site
yBounds = [0; 0; 40; 40]; % Y bounds for each site

% draw background sites first
if drawAllSites % draw all 384 sites
    xRep = repmat(geometry(:,1)', 4, 1) + repmat(xBounds,1,length(geometry)); % X coordinates
    yRep = repmat(geometry(:,2)', 4, 1) + repmat(yBounds,1,length(geometry)); % Y coordinates
    figHandle = figure('Position',[100 100 400 1200]); % figure('Position',[left bottom width height])
    patch(xRep,yRep,[0.7,0.7,0.7]) % patch all sites in grey
    axis tight
else
    xRep = repmat(geometry(min(sites):end,1)', 4, 1) + repmat(xBounds,1,length(geometry(min(sites):end,1))); % X coordinates
    yRep = repmat(geometry(min(sites):end,2)', 4, 1) + repmat(yBounds,1,length(geometry(min(sites):end,1))); % Y coordinates
    figHandle = figure('Position',[100 100 400 1200]); % figure('Position',[left bottom width height])
    patch(xRep,yRep,[0.7,0.7,0.7]) % patch sites in the range in grey
    axis tight
end

hold on;

xSitesRep = repmat(geometry(sites,1)',4,1) + repmat(xBounds,1,length(sites)); % x coordinates for sites
ySitesRep = repmat(geometry(sites,2)',4,1) + repmat(yBounds,1,length(sites)); % y coordinates for sites
patch(xSitesRep,ySitesRep,Z) % patch all sites in grey
axis tight
hold off;
colormap jet

end

