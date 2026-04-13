function keyword_folder_paths = find_keyword_containing_folder(start_path, keyword, varargin)
    % Search for all folders that contain the keyword in the specified start_path
    % Input:
    %   start_path - the path to start searching from
    %   keyword - the substring to search for in folder names
    %   varargin - name-value pair arguments (e.g., 'recursive', true/false)
    % Output:
    %   keyword_folder_paths - a cell array of full paths to the folders containing 'keyword'

    assert(ischar(keyword), 'Keyword must be a character array');
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'recursive', true, @islogical);
    parse(p, varargin{:});
    recursive = p.Results.recursive;
    
    % Initialize an empty cell array to store the paths of matching folders
    keyword_folder_paths = {};
    
    % Get the list of subdirectories and files in the current directory
    folder_info = dir(start_path);
    
    % Filter out the current and parent directory links ('.' and '..')
    folder_info = folder_info(~ismember({folder_info.name}, {'.', '..'}));
    
    % Iterate through the folder_info to find directories
    for i = 1:length(folder_info)
        if folder_info(i).isdir
            % Check if the folder name contains the keyword
            if contains(folder_info(i).name, keyword)
                % If keyword is found, store the full path in the cell array
                keyword_folder_paths{end+1} = fullfile(start_path, folder_info(i).name); %#ok<AGROW>
            end

            % If recursive search is enabled, search subdirectories
            if recursive
                subfolder_path = fullfile(start_path, folder_info(i).name);
                subfolder_keyword_paths = find_keyword_containing_folder(subfolder_path, keyword, 'recursive', true);

                % Append any found paths from the subfolder search
                keyword_folder_paths = [keyword_folder_paths, subfolder_keyword_paths]; %#ok<AGROW>
            end
        end
    end
end
