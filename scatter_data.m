function scatter_data(dataset1, varargin)
numb_sets = length(varargin) + 1; 

figure; hold on; 

n = numel(dataset1); 

for k = 1:numb_sets
    if k == 1
        x_jitter = (rand(length(dataset1), 1)-0.5)/10; 
        x = k*ones(length(dataset1), 1) + x_jitter; 
        y = dataset1; 
        scatter(x, y, 'k'); 
        pre_data = {x, y}; 
    else
        assert(n==numel(varargin{k-1}))
        x_jitter = (rand(length(dataset1), 1)-0.5)/10; 
        x = k*ones(length(varargin{k-1}), 1) + x_jitter;
        y = varargin{k-1}; 
        scatter(x, y, 'k'); 
        % draw line plot
        for kk = 1:n
            plot([pre_data{1}(kk), x(kk)], [pre_data{2}(kk), y(kk)], 'k:')
        end
        pre_data = {x, y}; 
    end
end
hold off; 
xlim([0.5 numb_sets+0.5])

end 