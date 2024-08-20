%% prepare color
pastel1 = slanCM('Pastel1', 10); 
plotColorListWithNumbers(pastel1); 
pastel1(end, :) = [0.85 0.85 0.85]; 

% Compiling the number of training sessions per animal
trainDayC{1, 1} = 'm25'; 
trainDayC{1, 2} = 40; 

trainDayC{2, 1} = 'm28'; 
trainDayC{2, 2} = 35; 

trainDayC{3, 1} = 'm37'; 
trainDayC{3, 2} = 28; 

trainDayC{4, 1} = 'm38'; 
trainDayC{4, 2} = 23; 

trainDayC{5, 1} = 'm39'; 
trainDayC{5, 2} = 24; 

trainDayC{6, 1} = 'm40'; 
trainDayC{6, 2} = 20; 

trainDayC{7, 1} = 'm44'; 
trainDayC{7, 2} = 10; 

trainDayC{8, 1} = 'm45'; 
trainDayC{8, 2} = 20; 

trainDayC{9, 1} = 'm46'; 
trainDayC{9, 2} = 16; 

% Descriptive stat
dat = cell2mat(trainDayC(:, 2)); 
[meanNT, stdNT, semNT] = meanstdsem(dat); 

idx = randperm(length(dat));
shufDat = dat(idx); 

% plot
figure;
hBar = bar(shufDat, 'BarWidth', 0.5, 'EdgeColor', 'none');
for k = 1:size(trainDayC)
    hBar.FaceColor = 'flat';  % Enable individual bar colors
    hBar.CData(k, :) = pastel1(k, :);  % Apply color to each bar
end

xlim([0 10])
pbaspect([1 1 1])
set(gca, 'TickDir', 'out')

print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/trainingDays'), '-dpdf', '-vector')

