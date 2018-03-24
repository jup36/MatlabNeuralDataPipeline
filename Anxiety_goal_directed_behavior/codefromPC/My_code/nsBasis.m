function basisStruct = nsBasis(var_name, var_range, flag, Fs, wd)
%% create natural spline basis for specific variable
% now only 'psth', 'srf', 'osc' can be analyzed 

if ~exist('flag', 'var')||isempty(flag)
    flag = 0;
end

if ~exist('Fs', 'var')||isempty(Fs)
    %sampling rate 
    Fs = 1000;  
end

if strcmp(var_name, 'psth')
    %psth basis 
    if ~exist('wd', 'var')||isempty(wd)
        %smoothing width 
        wd = 0.1;
    end
    basisStruct = nsPsth(var_range, flag, Fs, wd);
    
elseif strcmp(var_name, 'srf')
    %self-recovery function basis 
    basisStruct = nsSrf(var_range, flag, Fs);
    
elseif strcmp(var_name, 'osc')
    %oscillation basis
%     if ~exist('var_range', 'var')||isempty(var_range)
%         var_range = [-pi, pi];
%     end
    basisStruct = nsOsc(flag, Fs);
else
    disp('no related variable');
end









