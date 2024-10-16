function keyword_file_paths = find_keyword_file(start_path, keyword, subFolderLogic)
    % Recursively search for all files that contain the keyword in the name
    % within the specified start_path.
    % Inputs:
    %   start_path - the path to start searching from
    %   keyword - the keyword to search for in file names
    %   subFolderLogic - boolean to indicate whether to search in subfolders
    % Output:
    %   keyword_file_paths - a cell array of full paths to the files containing the keyword

    assert(ischar(keyword), 'Keyword must be a character array');
    assert(islogical(subFolderLogic), 'subFolderLogic must be a boolean');
    
    % Initialize an empty cell array to store the paths of matching files
    keyword_file_paths = {};
    
    % Get the list of subdirectories and files in the current directory
    folder_info = dir(start_path);
    
    % Filter out the current and parent directory links ('.' and '..')
    folder_info = folder_info(~ismember({folder_info.name}, {'.', '..'}));
    
    % Iterate through the folder_info to find files
    for i = 1:length(folder_info)
        if ~folder_info(i).isdir
            % Check if the file name contains the keyword
            if contains(folder_info(i).name, keyword)
                % If keyword is found, store the full path in the cell array
                keyword_file_paths{end+1} = fullfile(start_path, folder_info(i).name);
            end
        elseif subFolderLogic
            % If subFolderLogic is true, recursively search subdirectories
            subfolder_path = fullfile(start_path, folder_info(i).name);
            subfolder_keyword_paths = find_keyword_file(subfolder_path, keyword, subFolderLogic);
            
            % Append any found paths from the subfolder search
            keyword_file_paths = [keyword_file_paths, subfolder_keyword_paths]; %#ok<AGROW>
        end
    end
end
