% load model data
ch5ex1_problem_data;

runTimes = struct();

%% create the assembly objects
temporal = TemporalAssemblyNodal(pd, 0);
spatial  = SpatialAssemblySine(pd, spatialDim);
solver   = SolverNodal(pd, spatial, temporal);

%% create the parameter training set.
% two possibilities:
% first, equidistant grid
ptrain = linspace(pspan(1), pspan(2), paramNum);
% second, randomly chosen
% ptrain = rand(1, paramNum) * (pspan(2) - pspan(1)) - (pspan(2) - pspan(1)) / 2;

%% and now the real stuff: first let's execute the scm offline part
tic;
scm = SCM(solver, ptrain);
scm.save_scm_bounds = true;
scm.startOfflineStage();
runTimes.scm = toc;

%% and now the reduced basis offline part
tic;
rbm = RBM(pd, solver, scm);
rbm.resetTraining;
rbm.offlineStage(ptrain);
runTimes.rbm = toc;

% now we are ready to use the rbm solver in the online stage!
