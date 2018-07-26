function [indices] = TNC_FindEventsInContData(the_cont_data,min_space)

dlogic =0;

% determine the threshold
threshold       = mean(the_cont_data) + 4.*std(the_cont_data);

% find events
thr_cont_data   = the_cont_data>threshold;
tmp_indices     = find([0 diff(thr_cont_data)]==1);
cleans          = find(diff([0 tmp_indices])>min_space);
indices         = tmp_indices(cleans);

if dlogic==1
    figure(1); plot(the_cont_data,'k'); hold on; plot(indices,threshold.*ones(1,numel(indices)),'ro','MarkerSize',10);
end