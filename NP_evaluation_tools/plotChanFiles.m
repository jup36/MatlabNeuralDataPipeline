
function plotChanFiles( varargin )

% Plot event rate vs position on probe, either for a single *_chan.txt file
% or multiple files. When calling

if (numel(varargin) == 0)
    % set probe type and list of chan.txt (generated by detect_merge_peaks and
    % allSpike) files manually.
    % you can load as many *_chan.txt files as you like
    % If sites overlap, value at that site will be the last file loaded.
    pType = 1;
    fname{1} = "Z:\JP_NP24_test\catgt_WR39_100419_g0\WR39_100419_g0_tcat.imec.ap_chan.txt";
%     fname{1} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_shank0_bottomrow0_g0\anm_110_20191010_shank0_bottomrow0_g0_tcat.imec0.ap_chan.txt";
%     fname{2} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_shank0_bottomrow192_g0\anm_110_20191010_shank0_bottomrow192_g0_tcat.imec0.ap_chan.txt";
%     fname{3} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_shank1_bottomrow0_g0\anm_110_20191010_shank1_bottomrow0_g0_tcat.imec0.ap_chan.txt";
%     fname{4} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_shank1_bottomrow129_g0\anm_110_20191010_shank1_bottomrow129_g0_tcat.imec0.ap_chan.txt";
%     fname{5} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_shank2_bottomrow0_g0\anm_110_20191010_shank2_bottomrow0_g0_tcat.imec0.ap_chan.txt";
%     fname{6} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_shank2_bottomrow192_g0\anm_110_20191010_shank2_bottomrow192_g0_tcat.imec0.ap_chan.txt";
%     fname{7} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_shank3_bottomrow0_g0\anm_110_20191010_shank3_bottomrow0_g0_tcat.imec0.ap_chan.txt";
%     fname{8} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_shank3_bottomrow192_g0\anm_110_20191010_shank3_bottomrow192_g0_tcat.imec0.ap_chan.txt";
%     fname{9} = "Z:\CB_NP24\anm110\101019\catgt_anm_110_20191010_stripe_bottomrow384_g0\anm_110_20191010_stripe_bottomrow384_g0_tcat.imec0.ap_chan.txt";
%     fname{1} = "Z:\CB_NP24\anm112\101019\catgt_anm_112_20191010_shank0_bottomrow_0_g0\anm_112_20191010_shank0_bottomrow_0_g0_tcat.imec0.ap_chan.txt";
%     fname{2} = "Z:\CB_NP24\anm112\101019\catgt_anm_112_20191010_shank0_bottomrow_192b_g0\anm_112_20191010_shank0_bottomrow_192b_g0_tcat.imec0.ap_chan.txt";
%     fname{3} = "Z:\CB_NP24\anm112\101019\catgt_anm_112_20191010_shank1_bottomrow_0_g0\anm_112_20191010_shank1_bottomrow_0_g0_tcat.imec0.ap_chan.txt";
%     fname{4} = "Z:\CB_NP24\anm112\101019\catgt_anm_112_20191010_shank1_bottomrow_192_g0\anm_112_20191010_shank1_bottomrow_192_g0_tcat.imec0.ap_chan.txt";
%     fname{5} = "Z:\CB_NP24\anm112\101019\catgt_anm_112_20191010_shank2_bottomrow_0_g0\anm_112_20191010_shank2_bottomrow_0_g0_tcat.imec0.ap_chan.txt";
%     fname{6} = "Z:\CB_NP24\anm112\101019\catgt_anm_112_20191010_shank2_bottomrow_192_g0\anm_112_20191010_shank2_bottomrow_192_g0_tcat.imec0.ap_chan.txt";
%     fname{7} = "Z:\CB_NP24\anm112\101019\catgt_anm_112_20191010_shank3_bottomrow_0_g0\anm_112_20191010_shank3_bottomrow_0_g0_tcat.imec0.ap_chan.txt";
%     fname{8} = "Z:\CB_NP24\anm112\101019\catgt_anm_112_20191010_shank3_bottomrow_192_g0\anm_112_20191010_shank3_bottomrow_192_g0_tcat.imec0.ap_chan.txt";
else
    % read probe type and a single filename from the input arguments
    inputCell = varargin(1);
    pType = inputCell{1};
    for j = 2:numel(varargin)
        inputCell = varargin(2);
        fname{j-1} = inputCell{1};
    end
    
end

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
    XYPlot10( x, y, evtRate);  
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
        xlim([minX - xMargin, maxX + xMargin]);
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
function XYPlot10( xCoord, yCoord, data)

    nElec = 960;   %per shank; pattern repeats for the four shanks
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
    
    cMap = data;  

    f1 = figure('Units','Centimeters','Position',[15,4,10,18]);    

    % plot all positions
    scatter( elecPos(:,1), elecPos(:,2), 150, 'k', 'square' ); hold on;
    scatter( xCoord, yCoord, 100, cMap, 'square', 'filled' );hold on;   
    xlim([0,70]);
    ylim([-10,8000]);
    colormap cool;
    colorbar;
    xlabel('x position (um)');
    ylabel('z position (um)');
    title('3A or NP 1.0 shank view');
    hold off;
    
end % XY10Coord