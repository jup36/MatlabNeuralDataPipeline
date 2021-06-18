function [mapName] = TNC_CreateRBColormapJP(numValues,type)
% FUNCTION DETAILS: Simple utility to create a RWB style color map
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
%
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% modified by Junchol Park to add hot colormap

switch lower(type)
    
    case 'rb'
        mapLength = numValues;
        halfMap = floor(numValues./2);
        incr = 1./halfMap;
        increaser = (0:incr:1-incr)';
        tmp  = [ increaser, increaser, ones(halfMap,1)];
        tmp2 = [ ones(halfMap,1) , 1-increaser , 1-increaser ];
        mapName = [tmp;tmp2];
        
    case 'cb'
        mapLength = numValues;
        halfMap = floor(numValues./2);
        incr = 1./halfMap;
        incr2 = 0.33./halfMap;
        increaser = (0:incr:1-incr)';
        increaser2 = (0:incr2:0.33-incr2)';
        tmp  = [ increaser, 0.67+increaser2, ones(halfMap,1)];
        tmp2 = [ ones(halfMap,1) , 1-increaser , 1-increaser ];
        mapName = [tmp;tmp2];
        
    case 'wblack'
        incr = 1./numValues;
        increaser = (1-incr:-incr:0)'; %(0:incr:1-incr)';
        tmp  = [ increaser, increaser, increaser ];
        mapName = tmp;
        
    case 'wred'
        incr = 1./numValues;
        increaser = (1-incr:-incr:0)';
        %increaser = (0:incr:1-incr)';
        tmp  = [ ones(numValues,1), increaser, increaser ];
        mapName = tmp;
        
    case 'wblue'
        incr = 1./numValues;
        increaser = (1-incr:-incr:0)';
        tmp  = [ increaser, increaser, ones(numValues,1) ];
        mapName = tmp;
        
    case 'bo'
        
        if numel(numValues)==2
            incrL = 0.5./abs(numValues(1));
            increaserL = (0:incrL:0.5)';
            tmpL  = [ increaserL , ones(numel(increaserL),1).*0.25, 1-increaserL ];
            
            incrR = 0.5./abs(numValues(2));
            increaserR = (0.5:incrR:1)';
            tmpR  = [ increaserR , ones(numel(increaserR),1).*0.25, 1-increaserR ];
            
            mapName = [tmpL ; tmpR];
            
        else
            
            incr = 1./numValues;
            increaser = (0:incr:1-incr)';
            tmp  = [ increaser , ones(numValues,1).*0.5, 1-increaser ];
            mapName = tmp;
            
        end
        
    case 'cpb'
        
        if numel(numValues)==2
            incrL = 0.5./abs(numValues(1));
            increaserL = (0:incrL:0.5)';
            tmpL  = [ increaserL , ones(numel(increaserL),1).*0.67, 1-increaserL ];
            
            incrR = 0.5./abs(numValues(2));
            increaserR = (0.5:incrR:1)';
            tmpR  = [ increaserR , ones(numel(increaserR),1).*0.67, 1-increaserR ];
            
            mapName = [tmpL ; tmpR];
            
        else
            
            incr = 1./numValues;
            increaser = (0:incr:1-incr)';
            tmp  = [ increaser , (1-increaser).*0.67, 1-increaser ];
            mapName = tmp;
            
        end
        
        
    case 'mbr'
        mapName = [ 103,  0, 31;
            178, 24, 43;
            214, 96, 77;
            244,165,130;
            253,219,199;
            247,247,247;
            209,229,240;
            146,197,222;
            67,147,195;
            33,102,172;
            5, 48, 97 ] ./ 256;
        mapName = flipud(mapName);
        
    case 'gp'
        mapName = [ 64,  0, 75;
            118, 42,131;
            153,112,171;
            194,165,207;
            231,212,232;
            247,247,247;
            217,240,211;
            166,219,160;
            90,174, 97;
            27,120, 55;
            0, 68, 27] ./ 256;
        
        
    case 'map'
        mapName = [ 166,206,227;
            31,120,180;
            178,223,138;
            51,160, 44;
            251,154,153;
            227, 26, 28;
            253,191,111;
            255,127,  0;
            202,178,214;
            106, 61,154;
            255,255,153;
            177, 89, 40] ./ 256;
        
        
    case 'mapb'
        
        mapName = [ 31,120,180;
            
        51,160, 44;
        
        106, 61,154;
        
        227, 26, 28;
        
        255,127,  0] ./ 256;
    
    case 'hot'
        
        mapName = colormap(hot(numValues));
        
    case 'cool'
        
        mapName = colormap(cool(numValues));
        
    case 'parula'
        
        mapName = colormap(parula(numValues));
        
    case 'jet'
        
        mapName = colormap(jet(numValues));
        
    case 'summer'
        
        mapName = colormap(summer(numValues));
        
        
    case 'spring'
        
        mapName = colormap(spring(numValues));
        
    case 'winter'
        
        mapName = colormap(winter(numValues));
        
    case 'autumn'
        
        mapName = colormap(autumn(numValues));
        
    case 'gray'
        
        mapName = colormap(gray(numValues));
        
        
end
