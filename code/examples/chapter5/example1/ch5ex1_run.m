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
t = tic;
scm = SCM(solver, ptrain);
scm.save_scm_bounds = false;
scm.startOfflineStage();
runTimes.scm = toc(t);

%% and now the reduced basis offline part
t = tic;
rbm = RBM(pd, solver, scm);
rbm.rbm_save_error_bounds = false;
rbm.resetTraining;
rbm.offlineStage(ptrain, 1e-6);
runTimes.rbm = toc(t);

% now we are ready to use the rbm solver in the online stage!
% ch5ex1_online_query
