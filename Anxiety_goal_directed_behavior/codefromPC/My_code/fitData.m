function data = fitData(data, penalty, offset, maxIter)
%% fit glm model of neuron
fitVar = data.fitVar;       %factors to be fitted
numVar = length(fitVar);    %number of factors to be fitted

%time interval
temp = data.t0t1;
t0 = temp(1);
t1 = temp(2);
Fs = data.Fs;           %sampling frequency

if ~exist('penalty', 'var') || isempty(penalty)
    penalty = 0;
end
%response variable
if ~isfield(data, 'binMat')
    data.binMat = tsp2bin(data.tsp, t0, t1, Fs);    %binary matrix representing spike trains
end
Y = reshape(full(data.binMat), [], 1);
Y = double(Y);

%offset
if ~exist('offset', 'var')
    offset = zeros(size(Y));
end
fitPsth = false;
for m=1:numVar
    %creat offset for each variable
    eval(sprintf('offset_%s=zeros(size(Y)); ', fitVar{m}));
    if strcmp(fitVar{m}, 'psth')
        fitPsth = true;
    end
end

%% fit
if ~exist('maxIter', 'var')
    maxIter = 1;
end
warning off all
delta = 0.0001;
extra_offset = 0;

for iter = 1:maxIter
    for m=1:numVar
        varName = fitVar{m};
        tmpData = eval(sprintf('data.%s', varName));
        offset = offset-eval(sprintf('offset_%s', varName));
        
        %         if ~isfield(tmpData, 'beta')
        %             if strcmp(varName, 'psth')
        %                 temp = smooth(double(full(mean(data.binMat, 2))), 50);
        %                 temp(temp<=0) = 10^(-7);
        %                 x = tmpData.basis.basisX;
        %                 x = [ones(size(x,1), 1), x];  %#ok<AGROW>
        %                 tmpData.beta = x\log(temp);
        %            %                 tmpData.beta = rand(size(tmpData.basis.basisX, 2)+1, 1);
        %
        %                 tic;
        %             else
        %             tmpData.beta = rand(size(tmpData.basis.basisX, 2)+1, 1);
        %             end
        %         end
        if ~isfield(tmpData, 'beta')
            tmpData.beta = rand(size(tmpData.basis.basisX, 2)+1, 1);
        end
        if tmpData.type==0
            %use irls.m to fit
            if isfield(tmpData, 'X')
                X = tmpData.X;
                if iter==maxIter
                    %after finishing fitting, delete this variable to
                    %decrease data size
                    tmpData = rmfield(tmpData, 'X');
                end
            else
                X = full(eval_basis(tmpData.xind, tmpData.basis.basisObj));
                X = X+(rand(size(X))-0.5)*delta;
                tmpData.X = X;
            end
            
            if isfield(tmpData, 'ind')
                b = irls(X(tmpData.ind, :), Y(tmpData.ind), tmpData.beta, offset(tmpData.ind), penalty);
            else
                b = irls(X, Y, tmpData.beta, offset, penalty);
            end
            
            tmpData.mod = exp(xbeta(tmpData.basis.basisX, b));
            tmpData.beta = b;
            tmp_offset = xbeta(X, b);
            offset = offset+tmp_offset;
            eval(sprintf('offset_%s=tmp_offset; ', varName));
        elseif tmpData.type==1
            x = tmpData.basis.basisX;
            x = x+(rand(size(x))-0.5)*delta;
            xind = tmpData.xind;
            if isfield(tmpData, 'ind')
                ind = tmpData.ind;
                b = irls_repeat(x, Y(ind), xind(ind), tmpData.beta, offset(ind), penalty);
            else
                b = irls_repeat(x, Y, xind, tmpData.beta, offset, penalty);
            end
            
            %             norm_b0 = log(mean(exp(xbeta(tmpData.basis.basisX, b))));
            %             norm_b0 = 0;
            %             b(1) = b(1)-norm_b0;
            if isfield(tmpData, 'ind')
               xind(~ind) = 1;
               tmp_offset = xbeta(x, b, xind);
               % norm_b0 = mean(tmp_offset(ind)); 
                norm_b0 = log(mean(exp(tmp_offset(ind))));
               %norm_b0 = log(mean(exp(xbeta(tmpData.basis.basisX, b))));
                b(1) = b(1)-norm_b0;
                tmp_offset = xbeta(x, b, xind);
                tmp_offset(~ind) = 0;
            else
             %   norm_b0 = log(mean(exp(xbeta(x, b, xind)))); 
                norm_b0 = log(mean(exp(xbeta(tmpData.basis.basisX, b))));
                %             norm_b0 = 0;
                b(1) = b(1)-norm_b0;
                tmp_offset = xbeta(x, b, xind);
            end
           % extra_offset = extra_offset+norm_b0;
            tmpData.mod = exp(xbeta(tmpData.basis.basisX, b));
            tmpData.beta = b;
            offset = offset+tmp_offset;
            eval(sprintf('offset_%s=tmp_offset; ', varName));
        end
        eval(sprintf('data.%s=tmpData; ', varName));
        eval(sprintf('data.%s.offset = tmp_offset; ', varName));
    end
end

data.offset = offset;
%% normalize fit
%data.extra_offset = extra_offset;

if fitPsth && (numVar>1)
    data_norm = data;
    data_norm.fitVar = {'psth'};
    data_norm = fitData(data_norm, penalty, offset-offset_psth);
    data.psth = data_norm.psth;
end
data.penalty = penalty;
