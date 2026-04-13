
filePaths = {'Z:\Rodent Data\dualImaging_parkj\m1613_jRGECO_GRABda', ...
             'Z:\Rodent Data\dualImaging_parkj\m1859_jRGECO_GRABda', ...   
             'Z:\Rodent Data\dualImaging_parkj\m1873_jRGECO_GRABda'};
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

% parseAuditoryGngTrials
for f = 1:length(filePaths)
    filePath = filePaths{f};
    subFolders = find_keyword_folder(filePath, 'task');
    for ff = 1:length(subFolders)
        tbytDatFileC = find_keyword_file(subFolders{ff}, '_tbytDat.mat', true);
        if ~isempty(tbytDatFileC)
           load(tbytDatFileC{1}, 'tbytDat')
           tbytDat = parseAuditoryGngTrials(tbytDat); 
           header = extract_date_animalID_header(tbytDatFileC{1}); 
           fprintf(sprintf("Saving tbytDat parseGng session #%d of file #%d!\n", ff, f))
           save(fullfile(fileparts(tbytDatFileC{1}), strcat(header, '_tbytDat_parseGng')), 'tbytDat')
        end
    end
end
