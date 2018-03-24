filenm=input('input the file name:  ','s');
saveas(gcf,sprintf('%s.fig',filenm));
saveas(gcf,sprintf('%s.png',filenm));
saveas(gcf,sprintf('%s.eps',filenm),'psc2');