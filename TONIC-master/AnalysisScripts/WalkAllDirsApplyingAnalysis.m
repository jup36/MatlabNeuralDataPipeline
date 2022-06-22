% COLLECT PROPERTIES ACROSS A POPULATION OF CELLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From a given list of directories... 
%     enter the directory and store the file names,
% 
%     load the data for learning and extinction,
% 
%     plugin type calculation:
%     
%         module 1 for extinction cell paper:
%             psth aligned to cs (raw spike count)
%             isi histograms (log spacing)
%             spike waveforms
%             smoothed, z-scored psth correlations
%             create a unique cell id
%             aligned licking data (raw lick events)
%             electrode id [distance estimate]
%             mouse id
%             
%         module 2 for other signaling properties
%             psth aligned to us (raw spike count) [learning only]
%             correlated psthZ for simultaneous units
%             isi properties (+classifier)
%             electrode id [some estimate of distance]
%             mouse id
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
disp(['________________Running WalkAllDir... Script_________________________']);
disp(['Directory: ' pwd]);
disp(['Timestamp: ' datestr(now)]);

fid = fopen('WalkDirs.log','w');

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
                disp(allFiles(j).name);   

                currSes = currSes+1;
                disp(['Current session: ' num2str(currSes)]);
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____________APPLY ANALYSIS TO INDIVIDUAL FILES HERE_____________         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    propUA = []; propCS = 0; propEL = 0; clear propNames; p = 1;

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
                                propUA = [propUA,m];
                                propNames.names(p).str = RD_Active.neurons{m}.name;
                                p = p+1;
                                fprintf(fid,[RD_Active.neurons{m}.name '\n']);
                            end
                        end
                    end
                        
%%%%%%%%%%%% THIS IS GENERALLY UNIQUE TO A GIVEN MOUSE OR COLLECTION OF     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% DATA AND HAS CHANGED OVER THE COURSE OF OUR EXPERIMENTS:       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% digFlag is a 3 place vector that gives the datatype for CS, US, and EL (1==digital)
% % % % % % 
% % % % % % % % % %%%%%%%%%
% % % % % % % % % %%%%%%%%% DA-Mb-05, DA-Ma-10 Settings...
% % % % % % % % % '-0008', '-0004', '-0006', '-0007', ainp9a, ainp10Uall exist
% % % % % % % % % in extinction only 7 and 8 are still there
% % % % % % for q = 1:size(RD_Active.neurons,1) 
% % % % % %     if findstr(RD_Active.neurons{q}.name,'ainp10U')==1
% % % % % %         propEL = q;
% % % % % %     end
% % % % % % end
% % % % % % digUS = '-0006';
% % % % % % 
% % % % % % finalUS = digUS;
% % % % % % finalEL = propEL;
% % % % % % digFlag = [0,1,0];
% % % % % % latencyShift    = 0;


% % % % % % %%%%%%%%%%%%
% % % % % % %%%%%%%%%%%% DA-Ma-05 Settings...
% % % % % % % CS: ainp9a, US:?, EL: '-0001'?, DS: ainp11?, FL: ainp10a
% % % % % % % firsLick = 'ainp10a'
% % % % % % % SETTINGS FOR MOST OF THE DA-MA-05 DATA
% % % % % % for q = 1:size(RD_Active.neurons,1) 
% % % % % %     if findstr(RD_Active.neurons{q}.name,'ainp11a')==1
% % % % % %         propFL = q;
% % % % % %     end
% % % % % % end
% % % % % % digEL = '-0002';
% % % % % % 
% % % % % % finalUS = -2;
% % % % % % finalEL = digEL;
% % % % % % digFlag = [0,0,1];
% % % % % % latencyShift    = 1;


% Needed for one of the Da-Ma-05 with different digital signals.
for q = 1:size(RD_Active.neurons,1) 
    if findstr(RD_Active.neurons{q}.name,'ainp11a')==1
        propFL = q;
    end
end
digEL = '-0011';

finalUS = -2;
finalEL = digEL;
digFlag = [0,0,1];
latencyShift    = 0;

% % % % % % 
% % % % % % %%%%%%%%%%%%
% % % % % % %%%%%%%%%%%% DA-05 Settings...
% % % % % % % CS: ainp9a, US:?, EL: '-0001', FL: ainp11a
% % % % % % propUS = -1;
% % % % % % for q = 1:size(RD_Active.neurons,1) 
% % % % % %     if findstr(RD_Active.neurons{q}.name,'ainp11a')==1
% % % % % %         propUS = q;
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % propEL = 0;
% % % % % % for r = 1:size(RD_Active.neurons,1) 
% % % % % %     if findstr(RD_Active.neurons{r}.name,'ainp12a')==1
% % % % % %         propEL = r;
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % finalUS = propUS;
% % % % % % finalEL = '-0001';
% % % % % % digFlag = [0,0,1];
% % % % % % latencyShift    = 1;


% % % % % % %%%%%%%%%%%%
% % % % % % %%%%%%%%%%%% DA-02 Settings...
% % % % % % % CS: ainp9a, US:ainp11a, EL: ?, FL: none
% % % % % % propUS = -2;
% % % % % % for q = 1:size(RD_Active.neurons,1) 
% % % % % %     if findstr(RD_Active.neurons{q}.name,'ainp11a')==1
% % % % % %         propUS = q;
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % finalUS = propUS;
% % % % % % finalEL = 0; %doesn't exist
% % % % % % digFlag = [0,0,0];
% % % % % % latencyShift    = 1;
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % ALWAYS the CS is on ainp9a
                    for n = 1:size(RD_Active.neurons,1) 
                        if findstr(RD_Active.neurons{n}.name,'ainp9a')==1
                            propCS = n;
                        end
                    end

                    % Make sure the necessary data is found
                    if propCS==0
                        disp('ERROR: did not find the CS signal on ainp9a');
                    end
                    
                    
                    % classify the session
                    if size(findstr(allFiles(j).name,'-ex')) > 0
                        PopData.session(currSes).sessClass = 'extinction';
                        extLogic =1;
                    else
                        PopData.session(currSes).sessClass = 'learning';
                        extLogic =0;
                        % count units:
                        unitCount = unitCount + size(propUA,2);
                    end
                    
                    % create a unique id for the session
                    PopData.session(currSes).sessId = [nameInMem(1:8) '__' dirPath '__' num2str(extLogic)];

                    disp(       ['Session id: ' PopData.session(currSes).sessId ' |  class: ' PopData.session(currSes).sessClass ' | units: ' num2str(size(propUA,2))]);
                    fprintf(fid,['Session id: ' PopData.session(currSes).sessId ' |  class: ' PopData.session(currSes).sessClass ' | units: ' num2str(size(propUA,2)) '\n']);
                    
                    % FLAG TO PREVENT COMPLETE DATA EXTRACTION
                    if preCheck == 0
                        
                        [TD_Active] = TNC_ReadToStdRecStruct(RD_Active,propUA,propCS,finalUS,finalEL,digFlag,extLogic);

                        % Store the general session event data:
                        PopData.session(currSes).events         = TD_Active.events;
                        PopData.session(currSes).latencyShift   = latencyShift;
                        PopData.session(currSes).events.EL.ts   = TD_Active.events.EL.ts;
                        
                        if latencyShift
                            PopData.session(currSes).events.CS.ts = TD_Active.events.CS.ts + 100;
                        end

                        kMax = size(TD_Active.unit,2);
                        disp(['Total units: ' num2str(kMax)]);

                        if size(PopData.session(currSes).events.EL.ts,1) < 1
                            fprintf(fid,'No licks were found in this file.\n');
                        end
                        
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
                        
                    else
                        
                        % store the propUA vector to check
                        testingNames.sess(j).names = propNames.names;
                        
                        if j==2
                            
                            disp('   ');
                            numNames = size(testingNames.sess(1).names,2);
                            
                            if  size(testingNames.sess(1).names,2) ~= size(testingNames.sess(2).names,2)
                                disp('!!!!!!!!!!!!!! ERROR TYPE 2 !!!!!!!!!!!!!!');
                                disp('Mismatch in number of names found ... Please check file');
                                disp('!!!!!!!!!!!!!! ERROR TYPE 2 !!!!!!!!!!!!!!');
                                
                                fprintf(fid,'!!!!!!!!!!!!!! ERROR TYPE 2 !!!!!!!!!!!!!!\n');
                                fprintf(fid,'Mismatch in number of names found ... Please check file\n');
                                fprintf(fid,'!!!!!!!!!!!!!! ERROR TYPE 2 !!!!!!!!!!!!!!\n');
                                error = 2;
                            else

                                for y = 1:numNames  
                                    if ~strcmp(testingNames.sess(1).names(y).str,testingNames.sess(2).names(y).str)
                                        disp('!!!!!!!!!!!!!! ERROR TYPE 1 !!!!!!!!!!!!!!');
                                        disp('Mismatch in specific name found ... Please check file');
                                        disp('!!!!!!!!!!!!!! ERROR TYPE 1 !!!!!!!!!!!!!!');
                                        disp([testingNames.sess(1).names(y).str ' ... ' testingNames.sess(2).names(y).str]);

                                        fprintf(fid,'!!!!!!!!!!!!!! ERROR TYPE 1 !!!!!!!!!!!!!!\n');
                                        fprintf(fid,'Mismatch in specific name found ... Please check file\n');
                                        fprintf(fid,'!!!!!!!!!!!!!! ERROR TYPE 1 !!!!!!!!!!!!!!\n');
                                        fprintf(fid,[testingNames.sess(1).names(y).str ' ... ' testingNames.sess(2).names(y).str]);
                                        error = 1;
                                    end

                                end
                                
                            end
                            
                            if error==0;
                                disp('Directory is confirmed to have matched data between learning and extinction.');
                                fprintf(fid,'Directory is confirmed to have matched data between learning and extinction.\n');
                            else
                                error=0;
                            end
                            
                            disp('   ');
    
                        end
                        
                    end
                 
%                     PopData.session(currSes)
                    disp('  ');
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

system('open -e WalkDirs.log');