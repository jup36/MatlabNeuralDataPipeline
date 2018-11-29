function  

vFiles = dir('**/*.avi'); % list all the video files
vFronFiles = vFiles(cellfun(@(c)contains(c,'cam0'), {vFiles(:).name})); % front cam files

for i = 1:length(vFronFiles)
    % videoInfo(i).info = videoReaderWrapper(fullfile(vFronFiles(i).folder, vFronFiles(i).name)); % this doesn't work!
    
    
    fprintf('completed video #%d\n', i);
end

mmread(fullfile(vFronFiles(1).folder, vFronFiles(1).name))

end