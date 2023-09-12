function organize_tByt_videos(filePath)
%% get tbytDat
tbytDatPath = GrabFiles_sort_trials('tbytDat', 0, {fullfile(filePath, 'Matfiles')}); 
load(fullfile(tbytDatPath{1}), 'tbytDat')

%% get raw video files and csv bodypart coordinate files
% raw video files
filePath_raw_vid = findFolderWithString(filePath, '_vid');
if isempty(filePath_raw_vid)
    filePath_raw_vid = uigetdir(filePath); 
end
vidFiles = GrabFiles_sort_trials('behvid', 0, {filePath_raw_vid}); 

% DLC bodypart coordinate csv files
filePath_cropped_vid = findFolderWithString(filePath, '_vid_cropped');
csvFiles = dir(fullfile(filePath_cropped_vid, 'behvid*.csv'));

if numel(tbytDat) == numel(csvFiles)
    csvNames = {csvFiles.name}; 
    csvTimeC = cellfun(@(a) hmsToDateTime(a(8:15)), csvNames, 'un', 0); 
    [~, csvTimeI] = sort([csvTimeC{:}]); % to ensure the csv files are in correct timely order
    csvNamesSorted = csvNames(csvTimeI); 
    [tbytDat(:).csvFile] = deal(csvNamesSorted{:}); 
else
    warning('The number of csv files does not match with the number of trials!') % If some videos are missing each videoes need to be assigned using the trial indices in their names (TO DO)
end

%% parse 
ellipse_areaC = cell(numel(tbytDat), 1); 
for t = 1:numel(tbytDat)

    csvTab = readDlcCsv(fullfile(filePath_cropped_vid, tbytDat(t).csvFile)); 
    csvTab = removevars(csvTab, 'bodyparts_coords'); % drop redundant info
    newVarNames = {'L_x', 'L_y', 'L_like', ...
               'LD_x', 'LD_y', 'LD_like', ... 
               'D_x', 'D_y', 'D_like', ... 
               'DR_x', 'DR_y', 'DR_like', ... 
               'R_x', 'R_y', 'R_like', ... 
               'RV_x', 'RV_y', 'RV_like', ... 
               'V_x', 'V_y', 'V_like', ... 
               'VL_x', 'VL_y', 'VL_like'};
    csvTab.Properties.VariableNames = newVarNames;



    ellipse_areaC{t} = parse_csv_pupil_coordinates(csvTable, tbytDat(t).faceCam); 
    fprintf('Frames in trial #%d is parsed for pupil size.\n', t);
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        ellipse_area(fr) = ellipse_fr.long_axis * ellipse_fr.short_axis * pi; 

    end
    
    


    end







end