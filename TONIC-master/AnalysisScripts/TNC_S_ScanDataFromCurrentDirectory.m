

%% Screen through data for:
% For each prior used animal:
% 1 >>> Extract behavior from all datasets
% 2 >>> Select datasets with at least 50 reaches/rewards
% 3 >>> Scan online sorted units to look for largest number of online sorted stuff 
% 4 >>> Rank by the number of good online sorted units?
% 5 >>> Go through that rank order and look for the best stability over the duration of the recording session using Ns5Previewer
% 6 >>> Preprocess 1-5 files per mouse

%% 1 >>> ITERATE THROUGH THE LIST OF FILES EXTRACTING BEHAVIOR FROM THE NS4 DATA

clear baseName fileData

% Create a list of files which contain 'wide' or 'buz', but not 'bo'
allValid = [dir('./217547*.ns5') ; dir('./datafile*.ns5')];
% allValid = [dir('./*ch*.ns5')];
numFiles = numel(allValid);
disp(' ');
disp(' ');
disp(' ');
disp('________________________________________________________________________');
disp('________________________________________________________________________');
disp('________________________________________________________________________');
disp(['Begin processing ' num2str(numFiles) ' files...']);
disp(' ');
disp(' ');
disp(' ');

% % Old files
% chan.x = 137;
% chan.y = 138;
% chan.lick = 139;
% chan.rew = 142;
% chan.thr = 143;
    

%Katie's rig
chan.x = 129;
chan.y = 130;
chan.lick = 131;
chan.rew = 122;
chan.thr = 137;
    
for j=1:numFiles
    
    % Get the basename of the file
    numChar = numel(allValid(j).name);
    fileData(j).name = allValid(j).name(1:numChar-4);

    
    % Extract the behavior data
    [ContData] = TNC_MoverBehaviorExtract(allValid(j).name,fileData(j).name,'ns4',chan);
    
    % Store the number of reaches in the particular file
    fileData(j).numReaches = ContData.behavior.reach.numReaches;
    
    % Save a display of the basic behavior data
%     print(h0,'-depsc',[baseName '_bh.eps']);

end

%% 2 >>> SELECT DATASETS WITH AT LEAST 50 TRIALS OF DATA (AND RECORDING & VIDEO?)

clear validSession
validSession.count = 0;
reachThresh = 28;

disp(' ');
disp('________________________________________________________________________');

for  j=1:numFiles
    
   disp([num2str(j) ' >>> ' fileData(j).name ' contains ' num2str(fileData(j).numReaches) ' reaches.']);

   reachCount(j) = fileData(j).numReaches;
   
    if fileData(j).numReaches>reachThresh
        fileData(j).valid = 1;
        validSession.count = validSession.count+1;
        validSession.list(validSession.count) = j;
    else
        fileData(j).valid = 0;
    end
      
end
disp('________________________________________________________________________');
disp(['This directory contains ' num2str(validSession.count) ' valid sessions (>' num2str(reachThresh) ' reaches) for further analysis.']);

figure(100); hist(reachCount,50);

%% 3 >>> WALK THROUGH THE RANK ORDERED LIST OF BEHAVIOR DATA AND CREATE A PRINT OUT OF THE ONLINE SORTED UNITS

disp('________________________________________________________________________');
disp('Previewing data from valid files.');
disp(' ');

for k=1:numel(validSession.list)

    tmpNs5info  = openNSx('report','read',[fileData(validSession.list(k)).name '.ns5'],'t:1:3','sec');
    numChan     = size(tmpNs5info.MetaTags.ElecLabel,1);
    numShanks   = ceil(numChan ./ 8);

    if numel(strfind(fileData(validSession.list(k)).name,'-w-')) > 0
        arrayName = 'NN_w64';
        numShanks = 8;
    elseif numel(strfind(fileData(validSession.list(k)).name,'buz278')) > 0
        arrayName = 'NN_b64';
        numShanks = 8;
    elseif numel(strfind(fileData(validSession.list(k)).name,'128')) > 0
        arrayName = 'SingleSites';
        numShanks = 16;
    else
        arrayName = 'SingleSites';       
    end
    
    disp(['Treating as ' num2str(numShanks) ' shanks of 8 electrode sites.']);

    TNC_SSPL_Ns5Previewer([fileData(validSession.list(k)).name '.ns5'],15,numShanks,2,arrayName)
    
    disp(['File ' fileData(validSession.list(k)).name ' >>> array type is: ' arrayName]);
    fileData(validSession.list(k)).estimate = input('Quality estimate (1-10):');
    
end

%%
disp(' ');
disp(' ');

disp('________________________________________________________________________');
disp('Summary report on candidate files...');
disp(' ');
disp(' ');
disp(' ');for k=1:validSession.count
    disp([num2str(k) ' >>> ' fileData(validSession.list(k)).name ' Q> ' num2str(fileData(validSession.list(k)).estimate) ' R> ' num2str(fileData(validSession.list(k)).numReaches) ]);    
end
disp(' ');
disp('________________________________________________________________________');

save ValidSessionInfo validSession fileData

%% 5 >>> CHECK STABILITY BY NS5PREVIEW OF THE RANK ORDERED LIST OF BEHAVIOR AND NUMBER OF SORTED UNITS

%% 6 >>> PREPROCESS UP TO 5 OF THE BEST DATA SESSIONS PER MOUSE

