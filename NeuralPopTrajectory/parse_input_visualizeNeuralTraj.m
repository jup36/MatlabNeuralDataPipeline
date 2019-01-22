function p = parse_input_visualizeNeuralTraj( filePath, fileNameNeuralTrj, saveNameTag, eventMarkersRelative, vargs )
%parse input, and extract name-value pairs for the main function 'corrNeuralTrajMovKinematics.m'

default_whatToPlot = 'Trials';
default_dimPlotBy = 1; % the default PC or GPFA dimension to sort the trials/folds by
default_numbTrials = 10; % # of trials to plot
default_numbFolds = 4; % # of folds to plot
default_selectTrials = 1:10; % default select trials to plot

default_lineWidth  = 2;  % default lineWidth to be used to plot the neural trajectories
default_markerSize = 10; % default markerSize to be used to plot the neural trajectories
default_saveNtrjFigs = false; % By default, there's no need to save the 3-d nTrj figures
default_alpha = 0.05;  % default alpha value to just plot the significant correlations between neural population trajectories and the behavioral kinematics
default_trajCmap = 'parula'; % the default trajectory colormap

p = inputParser; % create parser object

addRequired(p,'filePath'); % file directory
addRequired(p,'fileNameNeuralTrj'); % fileName for neural population trajectories
addRequired(p,'saveNameTag'); % saveName used to save the outcomes
addRequired(p,'eventMarkersRelative'); % timing for specific events to be marked on each neural trajectory

addParameter(p,'whatToPlot', default_whatToPlot) % three different plot modes: 'Trials', 'Folds', 'selectTrials'
addParameter(p,'dimPlotBy', default_dimPlotBy) % the default PC or GPFA dimension to sort the trials/folds by
addParameter(p,'numbTrials', default_numbTrials) % the # of trials to plot
addParameter(p,'numbFolds', default_numbFolds) % the # of folds to plot
addParameter(p,'selectTrials', default_selectTrials) % the select trials to plot

addParameter(p,'lineWidth', default_lineWidth)   % the lineWidth to be used to plot nTrjs
addParameter(p,'markerSize', default_markerSize) % the markerSzie to be used to plot nTrjs
addParameter(p,'saveNtrjFigs', default_saveNtrjFigs) % By default, there's no need to save the 3-d nTrj figures
addParameter(p,'alpha', default_alpha) % By default, the alpha level of 0.01 is used as the criterion for significant correlations
addParameter(p,'trajCmap', default_trajCmap) % the default trajectory colormap

parse(p, filePath, fileNameNeuralTrj, saveNameTag, eventMarkersRelative, vargs{:})

end