%% WALK THROUGH ALL FILES IN A DIRECTORY AND LOAD DATA INTO MEMORY

clear

preCheck        = 0;
error           = 0;
stableTrials    = 10;
winPsth         = 10;

RunTonicInfo.activeSet = 1;
RunTonicInfo.loadedFiles = 0;
totalFiles = 0;
extLogic =0;

% create a list of all the folders in the directory
allDirectories  = dir; % creates a structure, use the "name" field
dirs            = size(allDirectories,1);
currSes         = 0;
unitCount       = 0;

disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(['________________LOADING ALL NEX FILES IN DIRECTORY_________________________']);
disp(['Directory: ' pwd]);
disp(['Timestamp: ' datestr(now)]);

fid = fopen([date '_' num2str(round(rand(1).*100)) '.log'],'w');

fprintf(fid,['\n________________Running WalkAllDir... Script_________________________\n']);
fprintf(fid,['Directory: ' pwd '\n']);
fprintf(fid,['Timestamp: ' datestr(now) '\n']);

disp(' ');

% loop through all the directories
for i = 1:dirs
    
    % test for a valid directory
    if allDirectories(i).isdir == 1
        
        dirPath = allDirectories(i).name;
        if strncmp('.',dirPath,1)==0
            disp(['________________' allDirectories(i).name '_________________________']);
            fprintf(fid,['\n________________' allDirectories(i).name '_________________________\n']);
            
            % create a list of all the files
            allFiles = dir(sprintf('%s/*.nex',dirPath));
                        
            for j = 1:size(allFiles,1)
                disp(' ');
                disp(allFiles(j).name);   

                currSes = currSes+1;
                disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
                disp(['Current file index: ' num2str(currSes)]);
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____________APPLY ANALYSIS TO INDIVIDUAL FILES HERE_____________         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    clear propUnits; clear propEvents; p = 1; propUnits.unitInds = []; 

                    RunTonicInfo.activeSet = RunTonicInfo.loadedFiles;
                    RunTonicInfo.loadedFiles = RunTonicInfo.loadedFiles+totalFiles;

                    fileNameStr = [dirPath '/' allFiles(j).name];
                    tmpName = allFiles(j).name(1:strfind(allFiles(j).name,'.')-1);

                    % remove '-'s from the name
                    nonDashInds = strfind(tmpName,'-');
                    tmpName(nonDashInds) = '_';
                    % remove ' ' from the name
                    spaceInds = strfind(tmpName,' ');
                    tmpName(spaceInds) = '_';
                    nameInMem = tmpName;

                    RunTonicInfo.activeSet = RunTonicInfo.activeSet+1;
                    RunTonicInfo.fileNames(RunTonicInfo.activeSet).name     = nameInMem;
                    RunTonicInfo.fileNames(RunTonicInfo.activeSet).loadTime = datestr(now);

                    RD_Active = TNC_LoadData(0, 0, fileNameStr);

                    for m = 1:size(RD_Active.neurons,1)    
                        if findstr(RD_Active.neurons{m}.name,'elec') == 1% electrode

                            testi = findstr(RD_Active.neurons{m}.name,'U');
                            if size(testi,1) == 0
                                
                                propUnits.unitInds = [propUnits.unitInds,m];
                                propUnits.names(p).str = RD_Active.neurons{m}.name;
                                p = p+1;
                                
                                fprintf(fid,[RD_Active.neurons{m}.name '\n']);                                

                            end
                                                       
                        end
                        
                        if findstr(RD_Active.neurons{m}.name,'ainp9a')
                            disp(RD_Active.neurons{m}.name);
                        end
                    end
                    
                    % ALWAYS the tone is ainp9a and light is ainp13U
                    for n = 1:size(RD_Active.neurons,1) 
                        if findstr(RD_Active.neurons{n}.name,'ainp9a')==1
                            propEvents.toneInd = n;
                        end
                        
                        if findstr(RD_Active.neurons{n}.name,'ainp13U')==1
                            propEvents.lightInd = n;
                        end
                        
                    end
                    
                    % create a unique id for the session
                    PopData.session(currSes).sessId = [nameInMem '>>' num2str(currSes)];

                    disp(       ['Session id: ' PopData.session(currSes).sessId ' | units: ' num2str(size(propUnits.unitInds,2))]);
                    fprintf(fid,['Session id: ' PopData.session(currSes).sessId ' | units: ' num2str(size(propUnits.unitInds,2)) '\n']);
                    
                    % FLAG TO PREVENT COMPLETE DATA EXTRACTION

                    if preCheck == 0
                        
                        % for reference: [tonicDataStructure] = TNC_PhotoStimToStruct(dataStructure,unitArrayToLoad,toneID,liteID)
                        [TD_Active] = TNC_PhotoStimToStruct(RD_Active,propUnits.unitInds,propEvents.toneInd,propEvents.lightInd)

                        % Store the general session event data:
                        PopData.session(currSes).events         = TD_Active.events;
                        
                        kMax = size(TD_Active.unit,2);
                        disp(['Total units: ' num2str(kMax)]);
                        
                        for k=1:kMax

                            % create a unique id for the cell/session pair
                            PopData.session(currSes).unit(k).uID  = [num2str(currSes) '_' num2str(k)];    
                            PopData.session(currSes).unit(k).name = TD_Active.unit(k).name;    
                            PopData.session(currSes).unit(k).ts   = TD_Active.unit(k).ts;   

                            % extract and store waveforms
                            PopData.session(currSes).unit(k).WfMean = mean(TD_Active.unit(k).wf,2);
                            PopData.session(currSes).unit(k).WfStd  = std(TD_Active.unit(k).wf,0,2);

                            % calculate and store isis
                            [isi] = TNC_QuantISI(TD_Active.unit(k).ts);
                            PopData.session(currSes).unit(k).isi = isi;    

                        end

                    end

                    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
                    disp('  ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____________APPLY ANALYSIS TO INDIVIDUAL FILES HERE_____________         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end

            disp('  ');
        end
        
    end
    
end

disp(['Unit count: ' num2str(unitCount)]);
clear DA_Ma*

fprintf(fid,'_____________________________________________________________________________________\n');
fprintf(fid,['Unit count: ' num2str(unitCount) '\n']);
% close the file
fclose(fid);

disp(['Saving data as ' nameInMem(1:8) 'PD.mat']);
eval(['save ' nameInMem(1:8) 'PD PopData']);
