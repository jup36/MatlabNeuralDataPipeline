function keyword_file_paths = find_keyword_containing_files(start_path, keyword, varargin)
    % Search for all files that contain the keyword in the specified start_path
    % Input:
    %   start_path - the path to start searching from
    %   keyword - the substring to search for in file names
    %   varargin - name-value pair arguments (e.g., 'recursive', true/false)
    % Output:
    %   keyword_file_paths - a cell array of full paths to the files containing 'keyword'

    assert(ischar(keyword), 'Keyword must be a character array');
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'recursive', true, @islogical);
    parse(p, varargin{:});
    recursive = p.Results.recursive;
    
    % Initialize an empty cell array to store the paths of matching files
    keyword_file_paths = {};
    
    % Get the list of files and folders in the current directory
    folder_info = dir(start_path);
    
    % Filter out the current and parent directory links ('.' and '..')
    folder_info = folder_info(~ismember({folder_info.name}, {'.', '..'}));
    
    % Iterate through the folder_info to find files
    for i = 1:length(folder_info)
        full_path = fullfile(start_path, folder_info(i).name);

        if ~folder_info(i).isdir
            % Check if the file name contains the keyword
            if contains(folder_info(i).name, keyword)
                % If keyword is found, store the full path in the cell array
                keyword_file_paths{end+1} = full_path; %#ok<AGROW>
            end
        elseif recursive
            % If recursive search is enabled, search subdirectories
            subfolder_file_paths = find_keyword_containing_files(full_path, keyword, 'recursive', true);

            % Append any found paths from the subfolder search
            keyword_file_paths = [keyword_file_paths, subfolder_file_paths]; %#ok<AGROW>
        end
    end
end
