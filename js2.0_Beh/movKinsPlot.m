function movKinsPlot(movKins)
%This function takes the structure movKins and plots Js movement
% trajectories (pos, vel, accel with reach start/stop/maxVel points marked on the trajectories). 

switch movKins.trialType
    case 'sp' % for a successful pull
        subplot(2,2,1); plot(movKins.jsTrajmm); hold on; plot(movKins.pullStart,movKins.jsTrajmm(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.jsTrajmm(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.jsTrajmm(movKins.pullMaxVelI),'og'); hold off;
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pullStart,movKins.smJsVel(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsVel(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsVel(movKins.pullMaxVelI),'og'); hold off; 
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pullStart,movKins.smJsAcl(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsAcl(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsAcl(movKins.pullMaxVelI),'og'); hold off; 
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pullStart,movKins.periodicAbsVelSum(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.periodicAbsVelSum(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.periodicAbsVelSum(movKins.pullMaxVelI),'og'); hold off; 
        title('absVelSum10msCum')
    case 'ps' % for a push
        subplot(2,2,1); plot(movKins.jsTrajmm); hold on; plot(movKins.pushStart,movKins.jsTrajmm(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.jsTrajmm(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.jsTrajmm(movKins.pushMaxVelI),'og'); hold off;
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pushStart,movKins.smJsVel(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.smJsVel(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.smJsVel(movKins.pushMaxVelI),'og'); hold off; 
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pushStart,movKins.smJsAcl(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.smJsAcl(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.smJsAcl(movKins.pushMaxVelI),'og'); hold off; 
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pushStart,movKins.periodicAbsVelSum(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.periodicAbsVelSum(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.periodicAbsVelSum(movKins.pushMaxVelI),'og'); hold off; 
        title('absVelSum10msCum')        
    case 'pm' % for a premature pull
        subplot(2,2,1); plot(movKins.jsTrajmm); hold on; plot(movKins.pullStart,movKins.jsTrajmm(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.jsTrajmm(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.jsTrajmm(movKins.pullMaxVelI),'og'); hold off;
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pullStart,movKins.smJsVel(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsVel(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsVel(movKins.pullMaxVelI),'og'); hold off; 
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pullStart,movKins.smJsAcl(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsAcl(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsAcl(movKins.pullMaxVelI),'og'); hold off; 
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pullStart,movKins.periodicAbsVelSum(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.periodicAbsVelSum(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.periodicAbsVelSum(movKins.pullMaxVelI),'og'); hold off; 
        title('absVelSum10msCum')

    case 'pmpp' % for a premature pull then push
        subplot(2,2,1); plot(movKins.jsTrajmm); hold on; plot(movKins.pullStart,movKins.jsTrajmm(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.jsTrajmm(movKins.pullStop),'*r');  plot(movKins.pullMaxVelI, movKins.jsTrajmm(movKins.pullMaxVelI),'og'); 
        plot(movKins.pushStart,movKins.jsTrajmm(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.jsTrajmm(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.jsTrajmm(movKins.pushMaxVelI),'og'); hold off;
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel); hold on; plot(movKins.pullStart,movKins.smJsVel(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsVel(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsVel(movKins.pullMaxVelI),'og'); 
        plot(movKins.pushStart,movKins.smJsVel(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.smJsVel(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.smJsVel(movKins.pushMaxVelI),'og'); hold off; 
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl); hold on; plot(movKins.pullStart,movKins.smJsAcl(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.smJsAcl(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.smJsAcl(movKins.pullMaxVelI),'og');
        plot(movKins.pushStart,movKins.smJsAcl(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.smJsAcl(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.smJsAcl(movKins.pushMaxVelI),'og'); hold off; 
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); hold on; plot(movKins.pullStart,movKins.periodicAbsVelSum(movKins.pullStart),'*b'); plot(movKins.pullStop,movKins.periodicAbsVelSum(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, movKins.periodicAbsVelSum(movKins.pullMaxVelI),'og'); 
        plot(movKins.pushStart,movKins.periodicAbsVelSum(movKins.pushStart),'*b'); plot(movKins.pushStop,movKins.periodicAbsVelSum(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, movKins.periodicAbsVelSum(movKins.pushMaxVelI),'og'); hold off; 
        title('absVelSum10msCum')
    case 'to' % for a timeout trial
        subplot(2,2,1); plot(movKins.jsTrajmm);  
        title('jsTrajmm (mm)')
        subplot(2,2,2); plot(movKins.smJsVel);  
        title('smJsVel (mm/s)')
        subplot(2,2,3); plot(movKins.smJsAcl);  
        title('smJsAcl (mm/s^2)')
        subplot(2,2,4); plot(movKins.periodicAbsVelSum); 
        title('absVelSum10msCum')
end

end

