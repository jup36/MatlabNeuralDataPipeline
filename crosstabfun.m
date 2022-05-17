function [tbl, chi2stat, pval] = crosstabfun(tab1, tab2)
         chi_label = [repmat('tab1',length(tab1),1); repmat('tab2',length(tab2),1)];
         chi_dat = [ones(sum(tab1==1),1); ones(sum(tab1==0),1)+1; 
                    ones(sum(tab2==1),1); ones(sum(tab2==0),1)+1];
         [tbl,chi2stat,pval] = crosstab(chi_label,chi_dat);                      
end
