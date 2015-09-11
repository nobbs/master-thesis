% Run a simple online query (use the rb-method to solve for a given parameter)

with_plots = false;

% choose a random parameter if none given
parameter = rand(1) * (abs(pspan(1)) + abs(pspan(2))) - (abs(pspan(1)) + abs(pspan(2))) / 2;

% compute the exact solution
tex = tic;
exsol = solver.solve(parameter);
Tex = toc(tex);

% compute the rb solution and the a posteriori error bound
trb = tic;
[rbsol, aperr] = rbm.onlineQuery(parameter);
Trb = toc(trb);

% plot the difference between the solutions (optional)
if with_plots
    [tmesh, xmesh] = meshgrid(pd.tgrid, pd.xgrid);

    figure(1);
    mesh(tmesh, xmesh, abs(solver.evaluateSolution(exsol) - rbm.evaluateSolutionRb(rbsol)));
    title('Differenz zwischen Truth- und RB-Lösung');
    xlabel('Zeit, t');
    ylabel('Ort, x');

    figure(2);
    mesh(tmesh, xmesh, solver.evaluateSolution(exsol));
    title('Truth-Lösung');
    xlabel('Zeit, t');
    ylabel('Ort, x');

    figure(3);
    mesh(tmesh, xmesh, rbm.evaluateSolutionRb(rbsol));
    title('RB-Lösung');
    xlabel('Zeit, t');
    ylabel('Ort, x');
end

% compute the exact error and effectivity
exerr = sqrt(full((rbm.trialSnapshots * rbsol - exsol).' * solver.TrNorm * (rbm.trialSnapshots * rbsol - exsol)));
effec = aperr / exerr;

% output
disp(sprintf('Solved with truth and rb solver for parameter %f', parameter));
disp(sprintf(' - time for truth:   %f', Tex));
disp(sprintf(' - time for rbm:     %f', Trb));
disp(sprintf(' - a posteriori err: %e', aperr));
disp(sprintf(' - real err:         %e', exerr));
disp(sprintf(' - effectivity:      %f', effec));
