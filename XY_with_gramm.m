function XY_with_gramm(varargin)
numb_sets = length(varargin); 

x = {};
y = [];

for jj = 1:numb_sets
    n = numel(varargin{jj});
    x = [x; repmat({sprintf('set_%d', jj)}, n, 1)];
    y = [y; varargin{jj}];
end

clear g

g(1,1)=gramm('x',x,'y',y,'color',x);
g(1,2)=copy(g(1));
g(1,3)=copy(g(1));
g(2,1)=copy(g(1));
g(2,2)=copy(g(1));

%Raw data as scatter plot
g(1,1).geom_point();
g(1,1).set_title('geom_point()');

%Jittered scatter plot
g(1,2).geom_jitter('width',0.4,'height',0);
g(1,2).set_title('geom_jitter()');

%Averages with confidence interval
g(1,3).stat_summary('geom',{'bar','black_errorbar'});
g(1,3).set_title('stat_summary()');

%Boxplots
g(2,1).stat_boxplot();
g(2,1).set_title('stat_boxplot()');

%Violin plots
g(2,2).stat_violin('fill','transparent');
g(2,2).set_title('stat_violin()');

%These functions can be called on arrays of gramm objects
g.set_names('x','data_sets','y','y');
g.set_title('');

gf = copy(g);

figure('Position',[100 100 800 550]);
g.draw();

gf.set_title('');
figure('Position',[100 100 800 550]);
gf.coord_flip();
gf.draw();

hold off; 

end 