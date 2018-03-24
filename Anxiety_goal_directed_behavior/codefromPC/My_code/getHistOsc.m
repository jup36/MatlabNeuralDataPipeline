function histOsc = getHistOsc(oscPhase, binMat, oscBin, t0t1, Fs, nsmuth)
if ~exist('nsmuth', 'var')
    nsmuth = 1; 
end

if ~exist('Fs', 'var')
    Fs = 1000; 
end

if exist('t0t1', 'var')
   ind0 = floor(t0t1(1)*Fs)+1;  
   ind1 = ceil(diff(t0t1)*Fs); 
   oscPhase = oscPhase(ind0:ind1,:); 
   binMat = binMat(ind0:ind1, :); 
end
   
%get modulatory function with histogram method
if ~exist('oscBin', 'var')||isempty(oscBin)
    oscBin = linspace(-pi, pi, 100);
end
temp = hist(oscPhase(binMat), oscBin);
oscHist = temp/mean(temp); 
% temp(1) = temp(1)+temp(end);
% temp(end) = [];
% temp1 = smooth(repmat(temp, 1, 3), nsmuth);
% nn = length(temp);
% oscHist = temp1((nn+1):(2*nn));
% oscHist = oscHist/mean(oscHist);
% oscHist(end+1) = oscHist(1);

histOsc.x = oscBin; 
histOsc.y = oscHist; 