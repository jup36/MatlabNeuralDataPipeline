function getBehVariablesJs(startDir, checkSubDir)
    % first search starting directory
    fileList = {};
    if checkSubDir
        subpath = strsplit(genpath(startDir), ';');
        isOutDir = cellfun(@(x) contains(x, {'.'}), subpath); % I remove folder with dots.
        subpath(isOutDir) = [];
        nSub = length(subpath) - 1; % the last data is empty
        files = [];
        for iS = 1:nSub
            files = [files; dir(fullfile(subpath{iS}, fileType))];
        end
    else                         
        files = dir(fullfile(startingDirectory, fileType));
    end
        


end