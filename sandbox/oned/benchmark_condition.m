Mspan = 1:10;
Qspan = 1:10;

conds = zeros(length(Mspan), length(Qspan));
times = zeros(length(Mspan), length(Qspan));
density = zeros(length(Mspan), length(Qspan));

for ii = 1:length(Mspan)
    for jj = 1:length(Qspan)
        num_M = Mspan(ii);
        num_Q = Qspan(jj);
        
        ii
        jj
        tic;
        t_plot = 0;
        start;
        
        density(ii, jj) = nnz(B) / numel(B);
        s = svds(B, length(B));
        t = toc
        
        times(ii, jj) = t;
        conds(ii, jj) = max(s) / min(s); 
    end
end

