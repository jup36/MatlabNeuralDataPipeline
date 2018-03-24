function spk = extract_spk(events, t0t1)
t0 = t0t1(1);
t1 = t0t1(2);
[M,N] = size(events);
spk = cell(M,N);
for m=1:M
    for n=1:N
        tsp = events{m,n};
        tsp(tsp<=t0) = [];
        tsp(tsp>t1) = [];
        spk{m,n} = tsp-t0;
    end
end