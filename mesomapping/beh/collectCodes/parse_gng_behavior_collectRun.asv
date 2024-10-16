
filePaths = {'Z:\Rodent Data\dualImaging_parkj\m1237_GCAMP', ...
    'Z:\Rodent Data\dualImaging_parkj\m1092_jRGECO', ...
    'Z:\Rodent Data\dualImaging_parkj\m1094_jRGECO'};

redetectLogic = false;

for f = 1:length(filePaths)
    filePath = filePaths{f};
    subFolders = find_keyword_folder(filePath, 'task');
    for ff = 1:length(subFolders)
        tbytDatFileC = find_keyword_file(subFolders{ff}, '_tbytDat', true);
        if isempty(tbytDatFileC) || redetectLogic
            parse_tbyt_auditory_gng_behavior_auto(subFolders{ff}, ...
                'preToneWin', 1, 'postToneWin', 4);
        end
    end
end
