% Script to load the standard head-fixed behavior rig recording session

%% Basic data loading
chunkDuration = sprintf('t:%g:%g',0,150);

% Load the continuous data recording channels
filenamestr = '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.ns5'
dataSpks = openNSx('report','read',filenamestr, 't:0:300', 'sec');

% Load the continuous behavior monitor channels
filenamestr = '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.ns4'
dataAinp = openNSx('report','read',filenamestr,'e:137:141','t:0:300', 'sec');

% Load behavior event times and digital stamps
filenamestr = '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.nev'
dataEvents = openNEV(filenamestr,'report','read','nomat');

