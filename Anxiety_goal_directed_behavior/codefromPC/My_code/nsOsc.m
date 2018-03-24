function tmpBasis = nsOsc(data, flag, Fs)
%% creating periodic spline function 

if ~exist('Fs', 'var')
    Fs = 1000; 
end
if ~exist('flag', 'var')
    flag = 0; 
end

%choose knots
if length(flag)>1
    %knots are pre-assigned
    knots = flag;
elseif flag==0
   knots = linspace(-pi, pi, 7); 
%      tmpN = 7; 
%      dphi = 2*pi/tmpN; 
%      knots = (dphi/2-pi):dphi:(pi-dphi/2); 
elseif flag==1
    oscPhase = data.oscPhase;
    tsp = data.tsp;
    temp = data.t0t1;
    t0 = temp(1);
    t1 = temp(2);
    binMat = full(tsp2bin(tsp, t0, t1, Fs));
    histOsc = getHistOsc(oscPhase, binMat);
    figure;
    plot(histOsc.x, histOsc.y);
    [knots, ~] = ginput();
    close;
else
    knots = linspace(-pi, pi, flag); 
%     tmpN = flag; 
%     dphi = 2*pi/tmpN; 
%     knots = (-pi+dphi/2):dphi:(pi-dphi/2); 
end

knots = mod(knots+pi, 2*pi)-pi;
phi = linspace(-pi, pi, 100);
basisX = eval_circular_basis(phi, knots);
tmpBasis.knots = knots;
tmpBasis.phi = phi;
% tmpBasis.basisX = basisX; 
basisX = bsxfun(@minus, basisX, mean(basisX, 1)); 
tmpBasis.basisX = bsxfun(@times, basisX, 1./std(basisX)); 
