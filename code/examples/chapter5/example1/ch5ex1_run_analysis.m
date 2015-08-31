%% computation
% compute the exact scm bounds
t = tic;
textprogressbar('Computing exact values: ');
ex_values = zeros(paramNum, 1);
Ynorm = solver.TeNorm;
Xnorm = solver.TrNorm;
for pdx = 1:paramNum
    textprogressbar(pdx / paramNum * 100);
    Lhs            = solver.spacetimeSystemMatrix(ptrain(pdx));
    supNorm        = Lhs.' * (Ynorm \ Lhs);
    ex_values(pdx) = sqrt(eigs(supNorm, Xnorm, 1, 'sm'));
end
textprogressbar(' done!');
toc(t)

