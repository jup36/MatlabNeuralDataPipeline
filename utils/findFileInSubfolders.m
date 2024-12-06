function filePath = findFileInSubfolders(rootFolder, fileName)
% findFileInSubfolders searches for a file with a specific name
% within a folder and all its valid subfolders.
%
% INPUT:
%   rootFolder - The root directory to start searching from (string).
%   fileName   - The name of the file to search for (string).
%
% OUTPUT:
%   filePath   - The full path to the file if found, or empty if not found.

    % Initialize the output
    filePath = '';
    
    % Check if the root folder exists
    if ~isfolder(rootFolder)
        error('The specified root folder does not exist.');
    end
    
    % Get all items in the root folder
    items = dir(rootFolder);
    
    % Loop through each item in the folder
    for i = 1:length(items)
        % Skip dummy folders '.' and '..'
        if items(i).isdir && (strcmp(items(i).name, '.') || strcmp(items(i).name, '..'))
            continue;
        end
        
        % Construct the full path of the current item
        currentPath = fullfile(rootFolder, items(i).name);
        
        % If the item is a folder, recursively search it
        if items(i).isdir
            filePath = findFileInSubfolders(currentPath, fileName);
            if ~isempty(filePath)
                return; % Stop searching once the file is found
            end
        % If the item is a file, check if the name matches
        elseif strcmp(items(i).name, fileName)
            filePath = currentPath;
            return; % Stop searching once the file is found
        end
    end
end
