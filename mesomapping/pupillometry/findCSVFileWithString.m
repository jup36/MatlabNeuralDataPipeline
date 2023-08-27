function filePath = findCSVFileWithString(directory, searchString)
    % List all items in the directory
    dirContents = dir(directory);
    
    % Initialize an empty array to store file paths
    filePaths = {};
    
    % Loop through each item in the directory
    for i = 1:length(dirContents)
        itemName = dirContents(i).name;
        
        % Check if the item is a file and its name contains the search string
        if ~dirContents(i).isdir && contains(itemName, searchString) && endsWith(itemName, '.csv')
            filePaths{end+1} = fullfile(directory, itemName);
        end
    end
    
    % Check if any matching CSV files were found
    if ~isempty(filePaths)
        % Display the found file paths
        disp('Matching CSV files:');
        disp(filePaths);
        
        % Return the first matching file path (you can modify this if needed)
        filePath = filePaths{1};
    else
        disp('No matching CSV files found.');
        filePath = '';
    end
end