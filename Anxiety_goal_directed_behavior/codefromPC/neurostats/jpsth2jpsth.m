function J = jpsth2jpsth(J,instr,str,T,xpsth,ypsth,xpsthSD,ypsthSD)
% function J = jpsth2jpsth(J,instr,str,T,xpsth,ypsth,xpsthSD,ypsthSD)
%
% converts from one type of jpsth to another
%
% J is the jpsth matrix
%
% instr is a string describing the type of input jpsth
% the options are
%
% [] : uses default
% raw : (default) input J is a raw jpsth
% cov : input J is a shuffle corrected jpsth
% covub : input J is an unbiased shuffle corrected jpsth
% cor : input J is a shuffle corrected and normalized jpsth 
% corub: input J is an unbiased shuffle corrected and normalized jpsth
%
% str is a string describing the type of output jpsth
% the options are
%
% 'cov' : the shuffle corrected jpsth obtained by subtracting the
%         product of the psths; needs xpsth and ypsth; may need T
% 'covub' : T/(T-1) times the result of 'cov'; needs T, xpsth and ypsth
% 'cor': the result of 'cov' divided by the product of the psth standard
%        deviations; needs xpsth, ypsth, xpsthSD, ypsthSD; may need T
% 'corub': the result of 'covub' divided by the product of the psth standard
%        deviations (using the ub option); needs T, xpsth, ypsth, xpsthSD, ypsthSD
%        this should be same result (up to rounding errors) as cor 
%        (but you must provide different psthSD's for cor and corub) 
% 'covNS' : amplitude non-stationary correction to 'cov'
% 'covubNS' : amplitude non-stationary correction to 'covub'
% 'corNS' : amplitude non-stationary correction to 'cor'
% 'corubNS' : amplitude non-stationary correction to 'covub'
%
% the algorithm only goes from less complicated to more complicated, so,
% for example, instr = 'cor' and outstr = 'cov' are incompatible
%
% the folloing inputs are needed by certain of the conversions
% T = number of trials
% xpsth = bins2psth(x); ypsth = bins2psth(y);
% xpsthSD = bins2psthSD(x); ypsthSD = bins2psthSD(y); % cor option
% xpsthSD = bins2psthSD(x,'ub'); ypsthSD = bins2psthSD(y,'ub'); % corub option

if isempty(instr)
    instr = 'raw';
end

instr = lower(instr);

inrawF = strcmp(instr,'raw');
incovF = strcmp(instr,'cov');
incovubF = strcmp(instr,'covub');
incorF = strcmp(instr,'cor');
incorubF = strcmp(instr,'corub');

if ~(inrawF || incovF || incovubF || incorF || incorubF)
    error('unrecognized instr option')
end

str = lower(str);

covF = strcmp(str,'cov');
covubF = strcmp(str,'covub');
corF = strcmp(str,'cor');
corubF = strcmp(str,'corub');
covNSF = strcmp(str,'covns');
covubNSF = strcmp(str,'covubns');
corNSF = strcmp(str,'corns');
corubNSF = strcmp(str,'corubns');

if covNSF || covubNSF || corNSF || corubNSF
    error('option unimplemented')
end
if ~(covF || covubF || corF || corubF || covNSF || covubNSF || corNSF || corubNSF)
    error('unrecognized option')
end

% identical in and out
if incovF && covF, return, end
if incovubF && covubF, return, end
if (incorF || incorubF) && (corF || corubF), return, end

% unbiased switching only
if incovF && covubF, J = (T./(T-1)).*J; return, end
if incovubF && covF, J = ((T-1)./T).*J; return, end

% unbiased switching with conversion to cor (saves a multiply)
if (inrawF && corubF) || (incovF && corubF)
    xpsthSD = ((T-1)./T).*xpsthSD;
end
if incovubF && corF
    xpsthSD = (T./(T-1)).*xpsthSD;
end

% divide by zero prevention
if corF || corubF 
    xpsthSD = xpsthSD + sqrt(eps(xpsthSD));
    ypsthSD = ypsthSD + sqrt(eps(ypsthSD));
end

% raw to cov
if inrawF && covF
    J = J - double(xpsth)*(ypsth.'); 
    return
end
if inrawF && covubF
    J = (T./(T-1)).*(J - double(xpsth)*(ypsth.')); 
    return
end

% raw to cor
if inrawF && (corF || corubF)
    J = (J - double(xpsth)*(ypsth.'))./(xpsthSD*(ypsthSD.'));
    return
end

% cov to cor
if (incovF || incovubF) && (corF || corubF)
    J = J ./ (xpsthSD*(ypsthSD.'));
    return
end
