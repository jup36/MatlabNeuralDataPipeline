function folderPath = findSubFolders(rootFolder, folderName)
% findSubFolders searches for a folder with a specific name
% within a root folder and all its valid subfolders.
%
% INPUT:
%   rootFolder - The root directory to start searching from (string).
%   folderName - The name of the folder to search for (string).
%
% OUTPUT:
%   folderPath - The full path to the folder if found, or empty if not found.

    % Initialize the output
    folderPath = '';
    
    % Check if the root folder exists
    if ~isfolder(rootFolder)
        error('The specified root folder does not exist.');
    end
    
    % Check if the root folder is the target folder
    [~, currentFolderName] = fileparts(rootFolder);
    if strcmp(currentFolderName, folderName)
        folderPath = rootFolder;
        return;
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
            if strcmp(items(i).name, folderName)
                folderPath = currentPath;
                return; % Stop searching once the folder is found
            else
                folderPath = findSubFolders(currentPath, folderName);
                if ~isempty(folderPath)
                    return; % Stop searching once the folder is found
                end
            end
        end
    end
end
