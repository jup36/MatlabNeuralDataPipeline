% plot jsTraj
hold on; 
for t = 1:size(jsTime1k,2)
    if strcmpi(jsTime1k(t).trialType,'sp')
         plot(jsTime1k(t).trJsReady:jsTime1k(t).trJsReady+length(jsTime1k(t).movKins.sgJsTrajmm)-1, jsTime1k(t).movKins.sgJsTrajmm, 'b')
    end    
end

% plot licks
for t = 1:length(evtIdx1k.lickIdx)
    plot(evtIdx1k.lickIdx(t), .2, 'ro', 'MarkerSize', 3)
end
hold off;