function [mov] = TNC_WatchAllChanEvolve(data1,data2,data3,matDispArray);
% FUNCTION DETAILS: animates a cell array of matrices
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

dataDims = size(data1);
dispDims = size(matDispArray);

figure(2);
clf;

% Create a blue->white->red colormap
colormap = [0 0 1; 1 1 1; 1 0 0];

if length(matDispArray) > 1
    
    for index = 1:dispDims(1,2)
        subplot(3,dispDims(1,2),index);
        surf(data1{1,matDispArray(index)}(:,:));        
        view([0 90]);
        shading flat;
        axis off;
        zlim([-10 10]);

        subplot(3,dispDims(1,2),index+dispDims(1,2));
        surf(data2{1,matDispArray(index)}(:,:));        
        view([0 90]);
        shading flat;
        axis off;
        zlim([-10 10]);

        subplot(3,dispDims(1,2),index+dispDims(1,2)+dispDims(1,2));
        surf(data3{1,matDispArray(index)}(:,:));        
        view([0 90]);
        shading flat;
        axis off;
        zlim([-10 10]);

    end

    if nargout>0
        mov(1) = getframe;
    end
    
else
    
    for index=1:dataDims(1,2)
        surf(data1{1,index}(:,:));
        shading interp;
        title(num2str(index));
        zlim([-10 10]);
        drawnow;
        if nargout>0
            mov(index) = getframe;
        end
    end
    
%     if nargout>0
%         movie(mov,2);
%     end
    
end

