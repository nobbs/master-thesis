Nspan = 1:50;

conds = zeros(length(Nspan), 1);
times = zeros(length(Nspan), 1);
density = zeros(length(Nspan), 1);

for ii = 1:length(Nspan)
    num_M = Nspan(ii);
    num_Q = Nspan(ii);
    
    ii
    tic;
    t_plot = 0;
    B = start(num_M, num_Q);
    
    density(ii) = nnz(B) / numel(B);
    s = svds(B, length(B));
    t = toc
    
    times(ii) = t;
    conds(ii) = max(s) / min(s);
end

