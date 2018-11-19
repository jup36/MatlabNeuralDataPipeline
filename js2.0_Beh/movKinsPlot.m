function movKinsPlot(movKins)
%This function takes the structure movKins and plots Js movement
% trajectories (pos, vel, accel with reach start/stop/maxVel points marked on the trajectories). 

switch movKins.trialType
    case 'sp' % for a successful pull
        subplot(2,2,1); plot(movKins.smJsTraj); hold on; plot(movKins.pullStart,movKins.smJsTraj(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsTraj(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsTraj(movKins.pullMaxVelI),'og'); hold off;
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pullStart,movKins.smJsVel(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsVel(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsVel(movKins.pullMaxVelI),'og'); hold off; 
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pullStart,movKins.smJsAcl(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsAcl(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsAcl(movKins.pullMaxVelI),'og'); hold off; 
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pullStart,movKins.periodicAbsVelSum(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.periodicAbsVelSum(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.periodicAbsVelSum(movKins.pullMaxVelI),'og'); hold off; 
        title('absVelSum10msCum')
    case 'ps' % for a push
        subplot(2,2,1); plot(movKins.smJsTraj); hold on; plot(movKins.pushStart,movKins.smJsTraj(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.smJsTraj(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.smJsTraj(movKins.pushMaxVelI),'og'); hold off;
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pushStart,movKins.smJsVel(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.smJsVel(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.smJsVel(movKins.pushMaxVelI),'og'); hold off; 
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pushStart,movKins.smJsAcl(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.smJsAcl(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.smJsAcl(movKins.pushMaxVelI),'og'); hold off; 
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pushStart,movKins.periodicAbsVelSum(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.periodicAbsVelSum(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.periodicAbsVelSum(movKins.pushMaxVelI),'og'); hold off; 
        title('absVelSum10msCum')        
    case 'pm' % for a premature pull
        subplot(2,2,1); plot(movKins.smJsTraj); hold on; plot(movKins.pullStart,movKins.smJsTraj(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsTraj(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsTraj(movKins.pullMaxVelI),'og'); hold off;
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pullStart,movKins.smJsVel(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsVel(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsVel(movKins.pullMaxVelI),'og'); hold off; 
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pullStart,movKins.smJsAcl(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsAcl(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsAcl(movKins.pullMaxVelI),'og'); hold off; 
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pullStart,movKins.periodicAbsVelSum(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.periodicAbsVelSum(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.periodicAbsVelSum(movKins.pullMaxVelI),'og'); hold off; 
        title('absVelSum10msCum')

    case 'pmpp' % for a premature pull then push
        if ~isnan(movKins.pull.startI)&&~isnan(movKins.push.startI)
            subplot(2,2,1); plot(movKins.smJsTraj); hold on; plot(movKins.pull.startI,movKins.smJsTraj(movKins.pull.startI),'*b'); plot(movKins.pull.stopI,movKins.smJsTraj(movKins.pull.stopI),'*r');  plot(movKins.pull.maxVelI, movKins.smJsTraj(movKins.pull.maxVelI),'og');
            plot(movKins.push.startI,movKins.smJsTraj(movKins.push.startI),'*b'); plot(movKins.push.stopI,movKins.smJsTraj(movKins.push.stopI),'*r'); plot(movKins.push.maxVelI, movKins.smJsTraj(movKins.push.maxVelI),'og'); hold off;
            title('jsTrajmm (mm)')
            subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pull.startI,movKins.smJsVel(movKins.pull.startI),'*b'); plot(movKins.pull.stopI,movKins.smJsVel(movKins.pull.stopI),'*r'); plot(movKins.pull.maxVelI, movKins.smJsVel(movKins.pull.maxVelI),'og');
            plot(movKins.push.startI,movKins.smJsVel(movKins.push.startI),'*b'); plot(movKins.push.stopI,movKins.smJsVel(movKins.push.stopI),'*r'); plot(movKins.push.maxVelI, movKins.smJsVel(movKins.push.maxVelI),'og'); hold off;
            title('smJsVel (mm/s)')
            subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pull.startI,movKins.smJsAcl(movKins.pull.startI),'*b'); plot(movKins.pull.stopI,movKins.smJsAcl(movKins.pull.stopI),'*r'); plot(movKins.pull.maxVelI, movKins.smJsAcl(movKins.pull.maxVelI),'og');
            plot(movKins.push.startI,movKins.smJsAcl(movKins.push.startI),'*b'); plot(movKins.push.stopI,movKins.smJsAcl(movKins.push.stopI),'*r'); plot(movKins.push.maxVelI, movKins.smJsAcl(movKins.push.maxVelI),'og'); hold off;
            title('smJsAcl (mm/s^2)')
            subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pull.startI,movKins.periodicAbsVelSum(movKins.pull.startI),'*b'); plot(movKins.pull.stopI,movKins.periodicAbsVelSum(movKins.pull.stopI),'*r'); plot(movKins.pull.maxVelI, movKins.periodicAbsVelSum(movKins.pull.maxVelI),'og');
            plot(movKins.push.startI,movKins.periodicAbsVelSum(movKins.push.startI),'*b'); plot(movKins.push.stopI,movKins.periodicAbsVelSum(movKins.push.stopI),'*r'); plot(movKins.push.maxVelI, movKins.periodicAbsVelSum(movKins.push.maxVelI),'og'); hold off;
            title('absVelSum10msCum')
        elseif ~isnan(movKins.pull.startI)&&isnan(movKins.push.startI)
            subplot(2,2,1); plot(movKins.smJsTraj); hold on; plot(movKins.pull.startI,movKins.smJsTraj(movKins.pull.startI),'*b'); plot(movKins.pull.stopI,movKins.smJsTraj(movKins.pull.stopI),'*r');  plot(movKins.pull.maxVelI, movKins.smJsTraj(movKins.pull.maxVelI),'og');
            title('jsTrajmm (mm)')
            subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pull.startI,movKins.smJsVel(movKins.pull.startI),'*b'); plot(movKins.pull.stopI,movKins.smJsVel(movKins.pull.stopI),'*r'); plot(movKins.pull.maxVelI, movKins.smJsVel(movKins.pull.maxVelI),'og');
            title('smJsVel (mm/s)')
            subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pull.startI,movKins.smJsAcl(movKins.pull.startI),'*b'); plot(movKins.pull.stopI,movKins.smJsAcl(movKins.pull.stopI),'*r'); plot(movKins.pull.maxVelI, movKins.smJsAcl(movKins.pull.maxVelI),'og');
            title('smJsAcl (mm/s^2)')
            subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pull.startI,movKins.periodicAbsVelSum(movKins.pull.startI),'*b'); plot(movKins.pull.stopI,movKins.periodicAbsVelSum(movKins.pull.stopI),'*r'); plot(movKins.pull.maxVelI, movKins.periodicAbsVelSum(movKins.pull.maxVelI),'og');
            title('absVelSum10msCum')
        else
        end
            
    case 'to' % for a timeout trial
        subplot(2,2,1); plot(movKins.smJsTraj);  
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel);  
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl);  
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); 
        title('absVelSum10msCum')
end

end

