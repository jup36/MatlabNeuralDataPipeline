function [] = TNC_PlotNCTraces(neuroCubeDataFile)   

sampling_rate = 24000;   % Sampling rate (# samples per sec)

%% Build up the data array from structure
ncdata = load(neuroCubeDataFile);
num_points = size(ncdata.data,2);
time = num_points/sampling_rate;

figure(1);

% Plot spiketimes
min_x = Inf;
max_x = -Inf;
for i=1:length(ncdata.Spiketimes)
    for j=1:length(ncdata.Spiketimes{i})
        x = [ncdata.Spiketimes{i}(j) ncdata.Spiketimes{i}(j)]/1000;
        y = [-4 2];
        if i>1 || j>1
            hold on; axis tight;
        end
        axis tight; 
        plot(x, y, 'color', 'red');
        hold on;  
        if min_x > x
           min_x = x;
        end
        if max_x < x
           max_x = x;
        end
    end
end
disp([' min_x=' num2str(min_x) ' max_x =' num2str(max_x )]);

% Plot traces
x = [1:num_points]*time/num_points;
disp(['min(x)=' num2str(min(x)) ' max(x)=' num2str(max(x))]);
for i=1:4
    axis tight;
    y = ncdata.data(i, :)/500;  
    plot(x, y); 
    hold on;  
end
xlabel('Time, sec');
ylabel('Potential');
title('NeuroCube traces');

drawnow;
zoom xon;

