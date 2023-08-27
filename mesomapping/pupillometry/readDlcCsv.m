function dataTable = readDlcCsv(filePath)

%% Read the header part
% Open the file for reading
fileID = fopen(filePath, 'r');

% Read the first three lines
lines = cell(3, 1);
for i = 1:3
    lines{i} = fgets(fileID);
    if lines{i} == -1
        break;  % Stop if end of file is reached
    end
end

% Close the file
fclose(fileID);

% Split the contents of each line by commas
splitContents = [];
for i = 1:length(lines)
    if ~isempty(lines{i})
        splitContents = [splitContents; strsplit(lines{i}, ',')];
    end
end

combinedHeader = cellfun(@(a, b) strcat(a, '_', b), splitContents(2, :), splitContents(3, :), 'un', 0);

%% Read the data part
data = readmatrix(filePath);

% Create a table using the combined header and the read data
dataTable = array2table(data, 'VariableNames', combinedHeader);

end


