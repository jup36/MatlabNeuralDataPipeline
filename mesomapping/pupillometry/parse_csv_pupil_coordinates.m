function [ellipse_area, faceCamPt, ellipse_fit, eyeOpenPer] = parse_csv_pupil_coordinates(csvTable, faceCamPt, height)
    % match frame numbers
    if length(faceCamPt) > size(csvTable, 1)
        faceCamPt = faceCamPt(1:size(csvTable, 1)); 
    end
    
    ellipse_area = nan(size(csvTable, 1), 1); 
    eyeOpenPer = nan(size(csvTable, 1), 1); 
    
    % fit ellipse frame by frame
    for fr = 1:size(csvTable, 1)
        tabFr = csvTable(fr, :); 
        X = []; 
        Y = []; 
        % L coordinates (1)
        if tabFr.L_like >= 0.5
            X(end+1) = tabFr.L_x; 
            Y(end+1) = height-tabFr.L_y; 
        end
        % LD coordinates (2)
        if tabFr.LD_like >= 0.5
            X(end+1) = tabFr.LD_x; 
            Y(end+1) = height-tabFr.LD_y; 
        end
        % D coordinates (3)
        if tabFr.D_like >= 0.5
            X(end+1) = tabFr.D_x; 
            Y(end+1) = height-tabFr.D_y; 
        end
        % DR coordinates (4)
        if tabFr.DR_like >= 0.5
            X(end+1) = tabFr.DR_x; 
            Y(end+1) = height-tabFr.DR_y; 
        end
        % R coordinates (5)
        if tabFr.R_like >= 0.5
            X(end+1) = tabFr.R_x; 
            Y(end+1) = height-tabFr.R_y; 
        end
        % RV coordinates (6)
        if tabFr.RV_like >= 0.5
            X(end+1) = tabFr.RV_x; 
            Y(end+1) = height-tabFr.RV_y; 
        end
        % V coordinates (7)
        if tabFr.V_like >= 0.5
            X(end+1) = tabFr.V_x; 
            Y(end+1) = height-tabFr.V_y; 
        end
        % VL coordinates (8)
        if tabFr.VL_like >= 0.5
            X(end+1) = tabFr.VL_x; 
            Y(end+1) = height-tabFr.VL_y; 
        end
        
        eyeOpenPer(fr) = length(X)/8; 

        if length(X) > 5 && length(Y) > 5
            ellipse_fit = fit_ellipse(X, Y); 
            %ellipse_fit_status{fr} = ellipse_fit.status; 
            ellipse_area(fr) = ellipse_fit.long_axis * ellipse_fit.short_axis * pi; 
        end
    end
    
    end