Mspan = 1:15;
Qspan = 1:15;

conds    = zeros(length(Mspan), length(Qspan));
times    = zeros(length(Mspan), length(Qspan));
density  = zeros(length(Mspan), length(Qspan));
relconds = zeros(length(Mspan), length(Qspan));

for ii = 1:length(Mspan)
    for jj = 1:length(Qspan)
        num_M = Mspan(ii);
        num_Q = Qspan(jj);

        % Daten und Einstellungen laden
        [data, opt, omega] = setup(num_M, num_Q);

        ii
        tic;
        t_plot = 0;

        B = steifigkeitsmatrix(data, omega, opt);
        [NX, NY] = normierungsmatrizen(opt);
        B = NY * B * NX;

        density(ii, jj) = nnz(B) / numel(B);
        s = svds(B, length(B));
        t = toc

        times(ii, jj) = t;
        conds(ii, jj) = max(s) / min(s);
        relconds(ii, jj) = conds(ii, jj) / (num_M * (num_Q + 1))^2;
    end
end
