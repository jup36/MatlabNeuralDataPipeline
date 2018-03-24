function [data, y] = getMod(varName, data)
%% get modulation function of one variable 
if ~exist('data', 'var')
    data = evalin('base', 'data'); 
end
tmpVar = eval(sprintf('data.%s', varName)); 
X = tmpVar.basisX; 
b = tmpVar.beta; 

y = exp(xbeta(X, b)); 
tmpVar.mod = y; 
eval(sprintf('data.%s=tmpVar; ', varName)); 
