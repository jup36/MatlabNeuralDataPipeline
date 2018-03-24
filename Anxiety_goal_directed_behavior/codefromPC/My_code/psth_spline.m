function data = psth_spline(data, knots_flag)
%% prepare spline basis for computing psth

%knots_flag=0: select knots automatically;
%knots_flag=1: select knots manually
%knots_flag>1: select m=knots_flag equally spaced knots
if nargin<2
    knots_flag = 0;
end

%if you already compute psth basis, ask to continuing.
if isfield(data, 'psth')
    temp = input('you have computed psth basis already, do you want to do it again? (y/n):  ', 's');
    if ~strcmp(temp, 'y')
        return;
    end
end

%create binary matrix for representing spikes
t0t1 = data.t0t1;
t0 = t0t1(1);
if t0==0
    t0=0.001;
end
t1 = t0t1(2);
tsp = data.tsp;
trial_num = length(data.trial_id);

T = ceil(t1*1000)-ceil(t0*1000)+1;
% if t0t1==0
%     T = T-1; 
% end
%create bin_spk and one vector storing all spike times
tsp_total = zeros(T*trial_num/10,1); 
k = 0; 
if ~isfield(data, 'bin_spk')
    bin_spk = uint8(zeros(T,trial_num));
    for m=1:trial_num
        tmp_tsp = tsp{m};
        if isempty(tmp_tsp)
            continue;
        end
        if isnan(tmp_tsp(end))
            tmp_tsp(end) = [];
        else
            tsp{m}(end+1) = nan;
            data.tsp{m}(end+1) = nan;
        end
        tmp_tsp(isnan(tmp_tsp)) = [];
        tsp_total((k+1):(k+length(tmp_tsp))) = tmp_tsp; 
        k = k+length(tmp_tsp); 
        ind = ceil(tmp_tsp*1000)-ceil(t0*1000)+1;
        ind(ind>T) = [];
        ind(ind<1) = []; 
        bin_spk(ind,m) = 1;
    end
    data.bin_spk = bin_spk;
    data.tsp_total = tsp_total(1:k); 
end
bin_spk = double(data.bin_spk);

if knots_flag<0
    return; 
end
%% create basis
xx = 1:T;
nOrder = 4; %cubic spline

%choose knots
if knots_flag==0
    knots = linspace(xx(1), xx(end), min(ceil(T/100),30)); %1 knots in every 100 ms.
elseif length(knots_flag)>1
    knots = knots_flag; 
elseif knots_flag==1
    figure;
    yy = mean(bin_spk,2);
    plot(xx, smooth(yy,30));
    [tmp_x, ~] = ginput();
    close;
    knots = [xx(1), tmp_x', xx(end)] ;
elseif knots_flag>1
    knots = linspace(xx(1), xx(end), ceil(knots_flag)); %1 knots in every 100 ms.
end
knots(knots>T) = [];
knots(knots<1) = [];
data.psth.knots = knots;

%create basis
nBasis = length(knots)+nOrder-2;
basisObj = create_bspline_basis([xx(1),xx(end)], nBasis, nOrder, knots);
data.psth.basisObj = basisObj;

%expand data on basis space
trial_num = size(bin_spk,2);
tvec = repmat(xx', trial_num, 1);
X = full(eval_basis(tvec, basisObj));
basisT = full(eval_basis(xx, basisObj));
data.psth.X = X;
data.psth.basisT = basisT;






