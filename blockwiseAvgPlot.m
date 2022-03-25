function blockwiseAvgPlot(m_b_mat)
    cmat = {'b', 'b', 'r', 'r', 'b', 'b', 'r', 'r', 'b', 'b', 'r', 'r', 'b', 'b', 'r', 'r'}; 
    randX = [-.5 + rand(size(m_b_mat,1)*size(m_b_mat,2),1)].*0.25; 
    alph = [0.8, 0.4]; 
    x = cell(size(m_b_mat,1), size(m_b_mat,2)); 
    y = cell(size(m_b_mat,1), size(m_b_mat,2)); 
    
    figure; hold on
    val = 1; 
    for b = 1:size(m_b_mat,2)
        for m = 1:size(m_b_mat,1)
            x{m,b} = b+randX(val); 
            y{m,b} = m_b_mat(m,b); 
            if ~isnan(y{m,b})
                scatter(x{m,b}, y{m,b}, 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmat{b}, ...
                'MarkerFaceAlpha', alph(mod(b,2)+1))
                % connect dots across blocks
                if b > 1
                   plot([x{m,b-1}, x{m,b}], [y{m,b-1}, y{m,b}], 'k:')    
                end
                val = val + 1; 
            end
        end
        avg_y_b = nanmean([y{:,b}]); 
        % mark the block average
        plot([b-.3, b+.3], [avg_y_b, avg_y_b], cmat{b}, 'LineWidth', 2) 
    end
    hold off
    set(gca, 'TickDir', 'out')
    
    xlim([.5 size(m_b_mat,2)+.5])
    ypad = floor((max(m_b_mat(:))-min(m_b_mat(:)))*0.08);
    ylim([min(m_b_mat(:))-ypad, max(m_b_mat(:))+ypad])
end
