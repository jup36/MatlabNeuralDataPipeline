
function [dataOut] = plotChanFilesJP( pType, estZ, fname )

% Plot event rate vs position on probe, either for a single *_chan.txt file
% or multiple files. When calling

% read probe type and a single filename from the input arguments

% initialize variables
i = 0;

shank = [];
x = [];
y = [];
chan = [];
rms = [];
evtRate = [];

for iFile = 1:numel(fname)
    
    fid = fopen( fname{iFile}, 'r');
    %read header line
    tline = fgetl(fid);

    while ~feof(fid)
        tline = fgetl(fid);
        i = i + 1;
        cArray = textscan(tline, '%d\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f');
        shank(i) = cArray{1};
        x(i) = cArray{2};
        y(i) = cArray{3};
        chan(i)= cArray{4};
        rms(i) = cArray{5};
        evtRate(i) = cArray{6};
    end
    fclose(fid);

end

if pType == 1
    dataOut = XYPlot10( x, estZ(chan+1), evtRate, estZ, fname);  
elseif (pType == 21) || (pType == 24)
    XYPlot20(pType, shank, x, y, evtRate);
else
    fprintf( "unknown probe type\n");
    return
end
        
end


% =========================================================
% plot data on top of all electrode positions for 2.0
%
% 
function XYPlot20(pType, shankind, xCoord, yCoord, data)   

   
    nElec = 1280;   %per shank; pattern repeats for the four shanks
    vSep = 15;   % in um
    hSep = 32;

    elecPos = zeros(nElec, 2);   

    elecPos(1:2:end,1) = 0;         %sites 0,2,4...
    elecPos(2:2:end,1) = hSep;      %sites 1,3,5...

    % fill in y values        
    viHalf = (0:(nElec/2-1))';                %row numbers
    elecPos(1:2:end,2) = viHalf * vSep;       %sites 0,2,4...
    elecPos(2:2:end,2) = elecPos(1:2:end,2);  %sites 1,3,5...
    
    cMap = data;  

    f1 = figure('Units','Centimeters','Position',[15,4,10,18]);
    if pType == 21
        % single shank probe. Plot only lowest selected electrode        
        % plot all positions
        minX = min(elecPos(:,1));
        maxX = max(elecPos(:,1));
        xMargin = 8;
        scatter( elecPos(:,1), elecPos(:,2), 150, 'k', 'square' ); hold on;
        scatter( xCoord, yCoord, 100, cMap, 'square', 'filled' );hold on;   
        xlim([minX - xMargin, maxX+xMargin]);
        ylim([-10,10000]);
        colormap cool;
        colorbar;
        xlabel('x position (um)');
        ylabel('z position (um)');
        title('NP 2.0 single shank view');
        hold off;
    else
        % four shank probe, no multiple connections   
        % replace xCoordinates with shank coordinates (more useful plots)
        elecPos(:,1) = elecPos(:,1)/hSep; % normalize
        offset = 0.1;
        elecPos(:,1) = elecPos(:,1)*2*offset - offset;         
        shankSep = 250;
        newX = xCoord - shankind*shankSep;
        newX = newX/hSep;
        newX = newX*2*offset - offset;
        
        for sI = 0:3
            cc = find(shankind == sI);
            scatter( sI + elecPos(:,1), elecPos(:,2), 30, 'k', 'square' ); hold on;
            scatter( sI + newX(cc), yCoord(cc), 20, cMap(cc), 'square', 'filled' ); hold on; 
        end
        xlim([-2*offset, 3+2*offset]);
        ylim([-10,10000]);
        colormap cool;
        colorbar;
        title('NP2.0 MS shank view');
        xlabel('shank index');
        ylabel('z position (um)');
        xticks([0 1 2 3])
        hold off;
    end

      
end % XY20Coord

% =========================================================
% plot data on top of all electrode positions for 2.0
%
%
function [dOut] = XYPlot10( xCoord, yCoord, data, estZ, fname )

    nElec = 384; %960;   %per shank; pattern repeats for the four shanks
    vSep = 20;   % in um
    hSep = 32;

    elecPos = zeros(nElec, 2);
    
    elecPos(1:4:end,1) = hSep/2;            %sites 0,4,8...
    elecPos(2:4:end,1) =  (3/2)*hSep;       %sites 1,5,9...
    elecPos(3:4:end,1) = 0;                 %sites 2,6,10...
    elecPos(4:4:end,1) =  hSep;             %sites 3,7,11...
    elecPos(:,1) = elecPos(:,1) + 11;       %x offset on the shank
    
    % fill in y values        
    viHalf = (0:(nElec/2-1))';                %row numbers
    elecPos(1:2:end,2) = viHalf * vSep;       %sites 0,2,4...
    elecPos(2:2:end,2) = elecPos(1:2:end,2);  %sites 1,3,5...
    
    % convert y values 
    convY = -estZ; 
    elecPos(:,2) = convY; 
    cMap = data;  

    f1 = figure('Units','Centimeters','Position',[15,4,10,18]);    

    % plot all positions
    scatter( elecPos(:,1), elecPos(:,2), 150, 'k', 'square' ); hold on;
    scatter( xCoord, -yCoord, 100, cMap, 'square', 'filled' );hold on;   
    xlim([0,70]);
    ylim([min(convY)-200,max(convY)+200]);
    %axis tight
    colormap cool;
    colorbar;
    xlabel('x position (um)');
    ylabel('z position (um)');
    title('3A or NP 1.0 shank view');
    hold off;
    caxis([0 15])
    
    [fDir,fName,~] = fileparts(fname{1}); 
    fileNameMarker = strfind(fName,'_'); 
    print(fullfile(fDir,strcat('eRateProbeMap_',fName(1:fileNameMarker(2)-1))),'-dpdf','-painters')
    
    dOut.elecPos = elecPos; 
    dOut.xCood = xCoord; 
    dOut.yCood = yCoord; 
    dOut.cMap = cMap; 
    
end % XY10Coord
