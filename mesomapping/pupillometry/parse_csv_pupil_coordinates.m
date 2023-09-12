    function [ellipse_area, faceCamPt] = parse_csv_pupil_coordinates(csvTable, faceCamPt)
    % match frame numbers
    if length(faceCamPt) > size(csvTable, 1)
        faceCamPt = faceCamPt(1:size(csvTable, 1)); 
    end
    
    ellipse_area = nan(size(csvTable, 1), 1); 

    % fit ellipse frame by frame
    for fr = 1:size(csvTable)
        tabFr = csvTable(fr, :); 
        X = []; 
        Y = []; 
        % L coordinates (1)
        if tabFr.L_like >= 0.5
            X(end+1) = tabFr.L_x; 
            Y(end+1) = tabFr.L_y; 
        end
        % LD coordinates (2)
        if tabFr.LD_like >= 0.5
            X(end+1) = tabFr.LD_x; 
            Y(end+1) = tabFr.LD_y; 
        end
        % D coordinates (3)
        if tabFr.D_like >= 0.5
            X(end+1) = tabFr.D_x; 
            Y(end+1) = tabFr.D_y; 
        end
        % DR coordinates (4)
        if tabFr.DR_like >= 0.5
            X(end+1) = tabFr.DR_x; 
            Y(end+1) = tabFr.DR_y; 
        end
        % R coordinates (5)
        if tabFr.R_like >= 0.5
            X(end+1) = tabFr.R_x; 
            Y(end+1) = tabFr.R_y; 
        end
        % RV coordinates (6)
        if tabFr.RV_like >= 0.5
            X(end+1) = tabFr.RV_x; 
            Y(end+1) = tabFr.RV_y; 
        end
        % V coordinates (7)
        if tabFr.V_like >= 0.5
            X(end+1) = tabFr.V_x; 
            Y(end+1) = tabFr.V_y; 
        end
        % VL coordinates (8)
        if tabFr.VL_like >= 0.5
            X(end+1) = tabFr.VL_x; 
            Y(end+1) = tabFr.VL_y; 
        end

        ellipse_fit = fit_ellipse(X, Y); 
        ellipse_area(fr) = ellipse_fit.long_axis * ellipse_fit.short_axis * pi; 

    end