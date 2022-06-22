function [data] = TNC_LoadOpenFieldCheapTracker(filename)


%% PARAMETERS OF ANALYSIS
% Create the filename list from the current directory
% fileDirectory = 'C:\Users\babita\Documents\HHMI\Mitopark\Open Field\MP1\'
% fileDirectory = '/Users/dudmanj/Documents/Work/Janelia/_PROJECTS/MOVER_CODE/OPEN_FIELD_ANALYSIS/'
% fileDirectory = '/Users/dudmanj/Documents/Work/Janelia/_PROJECTS/APPROACHINGNATURAL/'
% fileDirectory = './'
% fileDirectory = '~/subversion_working_copies/TONIC_v2/TrajectoryAnalysis/'
% Number of samples in running median
runMedWin   = 5;
offset = 0;
threshold = 5000; % choose an aggressive threshold because will also look for only single blips

%% DATA EXTRACTION

% % create a list of all the files
% allFiles = dir(sprintf('%s*.csv',fileDirectory));
% numFiles = size(allFiles,1);
 
    dateData = dlmread(filename,',',[0,0,0,4]);
    filenameBase = filename(1,1:numel(filename)-4);
    
    
    data.year   = dateData(1);
    data.month  = dateData(2);
    data.day    = dateData(3);
    data.hour   = dateData(4);
    data.min    = dateData(5);

    % Script to load csv behavior data
    clear disp* trial* 
 
    skipInitTrials = 0;

    disp(['________________Loading file: ' filename ' __________________________________'])

    tmpData = dlmread(filename,',',2,0);
    numSamps = size(tmpData,1);
    
    %___________________________________________________________
    % TIMESTAMPS (MS)
    data.tStamps = cumsum(tmpData(:,1)) + offset;
 
    %___________________________________________________________
    % DIGITAL SIGNALS FROM BEHAVIOR SYSTEM
    data.dig    = tmpData(:,2);

    %___________________________________________________________
    % SMOOTH AND NICE UP THE POSITION DATA
    data.posRaw.x = tmpData(:,4);
    data.posRaw.y = tmpData(:,5);
    data.posRaw.w = tmpData(:,6);
    data.posRaw.h = tmpData(:,7);

    %___________________________________________________________
    % EUCLIDEAN DISTANCE WITH SMOOTHED POSITION DATA
    i=2:numel(data.posRaw.x);
    data.velocity.velRaw(i) = sqrt( (data.posRaw.x(i)-data.posRaw.x(i-1)).^2 + (data.posRaw.y(i)-data.posRaw.y(i-1)).^2 );

    %___________________________________________________________
    % CLEAN UP THE VIDEO BLIPS POSITION DATA    
    allBlips = find(data.posRaw.w.*data.posRaw.h>threshold);
    figure(1); clf;
    plot(data.tStamps,data.posRaw.x,'k',data.tStamps(allBlips),data.posRaw.x(allBlips),'ro'); drawnow; hold on; pause();
    for p=1:numel(allBlips)
        data.posRaw.x(allBlips(p)) = mean([data.posRaw.x(allBlips(p)-1),data.posRaw.x(allBlips(p)+1)]);
        data.posRaw.y(allBlips(p)) = mean([data.posRaw.y(allBlips(p)-1),data.posRaw.y(allBlips(p)+1)]);
        data.posRaw.w(allBlips(p)) = mean([data.posRaw.w(allBlips(p)-1),data.posRaw.w(allBlips(p)+1)]);
        data.posRaw.h(allBlips(p)) = mean([data.posRaw.h(allBlips(p)-1),data.posRaw.h(allBlips(p)+1)]);
    end
    plot(data.tStamps(allBlips),data.posRaw.x(allBlips),'bo'); drawnow;

    %___________________________________________________________
    % POSITION DATA APPEARS TO BE OVERSAMPLED AND SO I NEED TO REMOVE ZEROS
    % AND THEN CREATE THE CONTINUOUS FUNCTION AGAIN WITH SMOOTHING
    data.posSmth.x = sgolayfilt(data.posRaw.x,11,31);
    data.posSmth.y = sgolayfilt(data.posRaw.y,11,31);    
        
    %___________________________________________________________
    % FURTHER SMOOTH POSITION DATA BY CREATING RUNNING MEDIANS    
    data.posSmth.Xrmed=medfilt1(data.posSmth.x,runMedWin);
    data.posSmth.Yrmed=medfilt1(data.posSmth.y,runMedWin);
             
    %___________________________________________________________
    % FURTHER SMOOTH VELOCITY DATA 
    data.velocity.vel(i)    = sqrt( (data.posSmth.Xrmed(i)-data.posSmth.Xrmed(i-1)).^2 + (data.posSmth.Yrmed(i)-data.posSmth.Yrmed(i-1)).^2 );
    data.velocity.velRM     = abs(sgolayfilt(data.velocity.vel,5,11));

    %___________________________________________________________
    % TAKE CUMULATIVE SUMS OF DISTANCE TO GET DISTANCE FROM ORIGIN
    data.velocity.cumdistance = (cumsum(data.velocity.vel))';
    
    % Write the data structure to a specific name
    disp(' ');
    disp(' ');
    disp([filenameBase '.data = data;']);
    eval([filenameBase '.data = data;']);
    eval(['save ' filenameBase ' ' filenameBase]);
    
    disp(['Output data in the structure: ' filename]);
    disp(' ');
   
    clear tmpData;

    figure(1); clf;
    plot(data.posSmth.Xrmed,data.posSmth.Yrmed,'k');
    

disp(['________________Completed loading file: ' filename ' __________________________________'])
disp(' ');
disp(' ');
disp(' ');
disp(' ');
