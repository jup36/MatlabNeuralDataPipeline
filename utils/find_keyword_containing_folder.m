function keyword_folder_paths = find_keyword_containing_folder(start_path, keyword, varargin)
    % Search for folders that contain or end with the keyword in start_path
    %
    % Input:
    %   start_path - the path to start searching from
    %   keyword - the substring to search for in folder names
    %
    % Name-value pairs:
    %   'recursive' - true/false, whether to search subdirectories
    %   'endsWithExactKeyword' - true/false
    %       false: return folders that contain keyword anywhere
    %       true:  return folders whose name ends exactly with keyword
    %
    % Output:
    %   keyword_folder_paths - cell array of full paths to matching folders

    assert(ischar(keyword) || isstring(keyword), 'Keyword must be a character array or string');

    % Convert keyword to char for compatibility with older MATLAB versions
    keyword = char(keyword);
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'recursive', true, @islogical);
    addParameter(p, 'endsWithExactKeyword', false, @islogical);
    parse(p, varargin{:});

    recursive = p.Results.recursive;
    endsWithExactKeyword = p.Results.endsWithExactKeyword;
    
    % Initialize output
    keyword_folder_paths = {};
    
    % Get the list of subdirectories and files
    folder_info = dir(start_path);
    
    % Filter out '.' and '..'
    folder_info = folder_info(~ismember({folder_info.name}, {'.', '..'}));
    
    % Iterate through folders
    for i = 1:length(folder_info)
        if folder_info(i).isdir

            folder_name = folder_info(i).name;

            % Decide matching rule
            if endsWithExactKeyword
                is_match = endsWith(folder_name, keyword);
            else
                is_match = contains(folder_name, keyword);
            end

            % Store matching folder path
            if is_match
                keyword_folder_paths{end+1} = fullfile(start_path, folder_name); %#ok<AGROW>
            end

            % Recursive search
            if recursive
                subfolder_path = fullfile(start_path, folder_name);
                subfolder_keyword_paths = find_keyword_containing_folder( ...
                    subfolder_path, keyword, ...
                    'recursive', true, ...
                    'endsWithExactKeyword', endsWithExactKeyword);

                keyword_folder_paths = [keyword_folder_paths, subfolder_keyword_paths]; %#ok<AGROW>
            end
        end
    end
end