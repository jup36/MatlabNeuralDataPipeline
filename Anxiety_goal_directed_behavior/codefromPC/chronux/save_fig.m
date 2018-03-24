function save_fig(name)
saveas(gcf, [name, '.eps'],'psc2'); 
saveas(gcf, [name, '.fig']); 
saveas(gcf, [name, '.png']); 