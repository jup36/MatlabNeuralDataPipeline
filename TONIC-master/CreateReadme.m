% SCRIPT DETAILS: A script to collect help descriptions from a directory (with subdirectories) of m-files.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% 

fid = fopen('README.txt','w');
basePath = pwd;

% create a list of all the folders in the directory
allHighLevel        = dir('*.m'); % creates a structure, use the "name" field
allHighLevelCount   = size(allHighLevel,1);

% First print the help of all the high level functions
for j=1:allHighLevelCount

    % write the function name, path
    if j==1 % write the directory name
        fprintf(fid, '\r_______________________TOP LEVEL______________________________________\r%s\r\r', basePath);
    end

    if length(strfind(allHighLevel(j).name,'.m'))>0
        if length(strfind(allHighLevel(j).name,'.m~'))==0
            helpStr = help(allHighLevel(j).name);
            helpStrTrunc = helpStr(1:strfind(helpStr,'_')-1);
            fprintf(fid, '%s\r', allHighLevel(j).name);
            fprintf(fid, '%s\r', helpStrTrunc);
        end
    end
    
end

% create a list of all the folders in the directory
allDirectories  = dir; % creates a structure, use the "name" field
dirs            = size(allDirectories,1);

% loop through all the directories
for i = 1:dirs
    
    % test for a valid directory
    if allDirectories(i).isdir == 1
        
        dirPath = allDirectories(i).name;
        if strncmp('.',dirPath,1)==0
            
            % create a list of all the files
            allFiles = dir(sprintf('%s/*.m',dirPath));
            
            for k = 1:size(allFiles,1)
                
                % for each directory loop through all the files
                currFuncName = allFiles(k).name;        

                % write the function name, path
                if k==1 % write the directory name
                    fprintf(fid, '\r______________________________________________________________________\r%s/%s\r\r', basePath, dirPath);
                end

                if length(strfind(allFiles(k).name,'.m'))>0
                    if length(strfind(allFiles(k).name,'.m~'))==0
                        helpStr = help(allFiles(k).name);
                        helpStrTrunc = helpStr(1:strfind(helpStr,'_')-1);
                        fprintf(fid, '%s\r', allFiles(k).name);
                        fprintf(fid, '%s\r', helpStrTrunc);
                    end
                end

            end
    
        fprintf(fid, '\r');
        end
        
    else % this is a file not a directory, now check if it is an m file and if so return the help
                
    end
end

% close the file
fclose(fid);
