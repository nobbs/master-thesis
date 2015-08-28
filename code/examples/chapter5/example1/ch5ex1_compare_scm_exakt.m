%% setup
% number of parameters:
nptest = 100;
% setup the paremeter discretization
ptest  = linspace(pspan(1), pspan(2), nptest);

%% computation
% compute the exact scm bounds
tic;
textprogressbar('Computing exact values: ');
ex_values = zeros(nptest, 1);
Ynorm = solver.TeNorm;
Xnorm = solver.TrNorm;
for pdx = 1:nptest
    textprogressbar(pdx / nptest * 100);
    Lhs            = solver.spacetimeSystemMatrix(ptest(pdx));
    supNorm        = Lhs.' * (Ynorm \ Lhs);
    ex_values(pdx) = sqrt(eigs(supNorm, Xnorm, 1, 'sm'));
end
textprogressbar(' done!');
toc

%% plotting
for idx = [1,2,3,10]
  figure(idx);
  plot(ptest, ex_values, '-', ...
       ptest, scm.scm_lower_bounds{idx}, '--', ...
       ptest, scm.scm_upper_bounds{idx}, '-.', ...
       scm.offMuCk(1:idx), sqrt(scm.offAlphaCk(1:idx)), '*');
  ylim([0.0920, 0.1020]);
  drawnow;
  matlab2tikz(['ch5ex1_scm_plot_', num2str(idx), '.tikz']);
end

