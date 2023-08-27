function folderPath = findFolderWithString(directory, searchString)
    % List all items in the directory
    dirContents = dir(directory);
    
    % Initialize an empty array to store folder names
    folderNames = {};
    
    % Loop through each item in the directory
    for i = 1:length(dirContents)
        itemName = dirContents(i).name;
        
        % Check if the item is a directory (folder) and its name contains the search string
        if dirContents(i).isdir && contains(itemName, searchString)
            folderNames{end+1} = fullfile(directory, itemName);
        end
    end
    
    % Check if any matching folders were found
    if ~isempty(folderNames)
        % Display the found folders
        %disp('Matching folders:');
        %disp(folderNames);
        
        % Return the first matching folder path (you can modify this if needed)
        folderPath = folderNames{1};
    else
        %disp('No matching folders found.');
        folderPath = '';
    end
end
