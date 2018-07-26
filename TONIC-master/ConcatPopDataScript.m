%% Combine PopData structures
% 
% Population data has the structure:
% PopData.session(i), where i is an indexed list of sessions
% 

% Run through file list and load data
% create a list of all the files
allFiles = dir('D*.mat');
disp('---------------------------------------------');
disp('Compiling into a complete dataset');
disp('---------------------------------------------');
disp(['Files to load: ' num2str(size(allFiles,1))]);
disp('---------------------------------------------');
    
for i = 1:size(allFiles,1)
    
    % load the file
    disp(allFiles(i).name);
    eval(['load ' allFiles(i).name]);
    
    % Create target structure:
    if exist('CompleteDataSet')
        sizeExist = size(CompleteDataSet.session,2);
    else
        sizeExist = 0;
    end

    sizeToAdd = size(PopData.session,2);

    for j=1:sizeToAdd

        CompleteDataSet.session(j+sizeExist) = PopData.session(j);

    end

    disp(['Current dataset size: ' num2str(size(CompleteDataSet.session,2))]);
    disp('  ');
    
% clear the last file and prepare for the next one
    clear PopData
    
end

disp('---------------------------------------------');
disp('done.');
disp('---------------------------------------------');
