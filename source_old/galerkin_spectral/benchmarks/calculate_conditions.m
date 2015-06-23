Nspan = 1:20;

conds    = zeros(length(Nspan), 1);
times    = zeros(length(Nspan), 1);
density  = zeros(length(Nspan), 1);
relconds = zeros(length(Nspan), 1);

for ii = 1:length(Nspan)
	ii

    num_M = Nspan(ii);
    num_Q = Nspan(ii);

    % Daten und Einstellungen laden
    [data, opt, omega] = setup(num_M, num_Q);

    tic;
    t_plot = 0;

    B = steifigkeitsmatrix(data, omega, opt);
    [NX, NY] = normierungsmatrizen(opt);
    B = NY * B * NX;

    density(ii) = nnz(B) / numel(B);
    s = svds(B, length(B));
    t = toc

    times(ii) = t;
    conds(ii) = max(s) / min(s);
    relconds(ii) = conds(ii) / (num_M * num_Q)^2;
end
