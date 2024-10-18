function keyword_folder_paths = find_keyword_folder_sortByDate(start_path, keyword)
    % Recursively search for all folders named keyword in the specified start_path
    % Input:
    %   start_path - the path to start searching from
    %   keyword - the folder name to search for
    % Output:
    %   keyword_folder_paths - a cell array of full paths to the folders named 'keyword'

    assert(ischar(keyword), 'Keyword must be a character array');
    
    % Initialize an empty cell array to store the paths of matching folders
    keyword_folder_paths = {};
    
    % Get the list of subdirectories and files in the current directory
    folder_info = dir(start_path);
    
    % Filter out the current and parent directory links ('.' and '..')
    folder_info = folder_info(~ismember({folder_info.name}, {'.', '..'}));
    
    % Iterate through the folder_info to find directories
    for i = 1:length(folder_info)
        if folder_info(i).isdir
            % Check if the folder name matches the keyword
            if strcmp(folder_info(i).name, keyword)
                % If keyword is found, store the full path in the cell array
                keyword_folder_paths{end+1} = fullfile(start_path, keyword);
            end
            
            % Recursively search this subdirectory
            subfolder_path = fullfile(start_path, folder_info(i).name);
            subfolder_keyword_paths = find_keyword_folder(subfolder_path, keyword);
            
            % Append any found paths from the subfolder search
            keyword_folder_paths = [keyword_folder_paths, subfolder_keyword_paths]; %#ok<AGROW>
        end
    end
    
    % After all folders are found, sort them based on the extracted date
    if ~isempty(keyword_folder_paths)
        headers = cellfun(@extract_date_animalID_header, keyword_folder_paths, 'UniformOutput', false);
        validHeaders = ~cellfun(@isempty, headers);  % Filter out invalid paths
        keyword_folder_paths = keyword_folder_paths(validHeaders);  % Keep only valid paths
        headers = headers(validHeaders);  % Filter headers

        % Use regexp to extract the 6-digit date part
        dateStrings = cellfun(@(x) regexp(x, '[0-9]{6}', 'match', 'once'), headers, 'UniformOutput', false);
        
        % Convert the extracted date strings to datetime objects
        dates = datetime(dateStrings, 'InputFormat', 'MMddyy');  % Convert to datetime
        
        % Sort by date
        [~, sortIdx] = sort(dates);
        keyword_folder_paths = keyword_folder_paths(sortIdx);  % Sort the folder paths by date
    end
end
