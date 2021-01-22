function [bump, bin, cos_func] = log_cos(n_bump, range, dt, bias, normalize)

if nargin < 2
    range = [0, 100];
end
if length(range) == 1
    range = [0, range];
end
assert(range(1) >=0, 'range should be larger than 0');

if nargin < 3
    dt = 1;
end

bias_min = 0.1;
if nargin < 4 || bias < bias_min
    bias = bias_min;
end
bias = bias * diff(range) / 100; % a temporal bias, e.g. 240ms

if nargin < 5
    normalize = true;
end

log_func = @(x) log(x + bias);
exp_func = @(x) exp(x) - bias;

range_log = log_func(range);
gap_log = diff(range_log) / (n_bump + 1);

bin = (range(1):dt:range(2))';
peak_log = range_log(1):gap_log:range_log(2)-2*gap_log;

cos_func = @(x) raised_cosine(x, peak_log, gap_log, log_func);
bump = cos_func(bin);

if normalize
    bump = bump ./ sum(bump, 1);
end


function bump = raised_cosine(x, peak, gap_log, log_func)
x = x(:);
n_x = length(x);

peak = peak(:)';
n_peak = length(peak);

x_mat = repmat(x, 1, n_peak);
peak_mat = repmat(peak, n_x, 1);

cos_func = @(x, peak) (cos(max(-pi, min(pi, (log_func(x) - peak) * pi / (2 * gap_log)))) + 1) / 2; % / (2 * gap_log): this is a clever way to sharpen the cosine functions
%cos_func1 = @(x, peak) (cos(max(-pi, min(pi, (log_func(x) - peak) * pi)))+1) / 2; % Try without this term ('/(2 * gap_log)')

bump = cos_func(x_mat, peak_mat);
%bump1 = cos_func1(x_mat,peak_mat); 

% fPath ='/Users/parkj/Documents/MATLAB/MatlabNeuralDataPipeline/reachPull_glm/understand_log_consine';
% % figure; plot(x); axis tight
% figure; imagesc(log_func(x_mat)); title('log_func(X_mat)','Interpreter','none'); colorbar; print(fullfile(fPath,'x_mat'),'-dpdf','-painters')
% figure; imagesc(peak_mat); title('log(peak_mat)','Interpreter','none'); colorbar; print(fullfile(fPath,'peak_mat'),'-dpdf','-painters')
% figure; imagesc((log_func(x) - peak)); title('Subtraction of peak from X','Interpreter','none'); colorbar; print(fullfile(fPath,'subtraction_X_Peak_mat'),'-dpdf','-painters')
% figure; imagesc((log_func(x) - peak)*pi); title('Subraction_multiplied_by_pi','Interpreter','none'); colorbar; print(fullfile(fPath,'subraction_multiplied_by_pi'),'-dpdf','-painters')
% figure; imagesc((log_func(x) - peak)*pi/(2*gap_log)); title('Subraction_multiplied_by_pi_expand','Interpreter','none'); colorbar; print(fullfile(fPath,'subraction_multiplied_by_pi_expand'),'-dpdf','-painters')
% figure; imagesc(min(pi, (log_func(x) - peak) * pi / (2 * gap_log))); title('Apply the upper limit of pi','Interpreter','none'); colorbar; print(fullfile(fPath,'Apply the upper limit of pi'),'-dpdf','-painters')
% figure; imagesc(max(-pi, min(pi, (log_func(x) - peak) * pi / (2 * gap_log)))); title('Apply the lower limit of pi','Interpreter','none'); colorbar; print(fullfile(fPath,'Apply the lower limit of pi'),'-dpdf','-painters')
% figure; imagesc((cos(max(-pi, min(pi, (log_func(x) - peak) * pi / (2 * gap_log)))))); title('Take cosine','Interpreter','none'); colorbar; print(fullfile(fPath,'Take cosine'),'-dpdf','-painters')
