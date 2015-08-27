% Chapter 5 - First example.
%
% This example corresponds to the first example given in chapter 4. We use
% homogeneous boundary conditions in space and the remaining model data is given
% by:
%   * `\Omega = [0, 1]`, `I = [0, 1]`,
%   * the parametrich differential operator `A(\sigma)` is given by
%       `A \eta = - \Delta \eta + \eta + \sigma \sin(\pi \blank) \eta.`,
%   * where the parameter domain is given by `\mathcal P = [-1, 1]`,
%   * the source term is `g = 0` and the initial conditions are given by
%       `u_0 = sin(\pi \blank)`.

%% setup
% dimensions of the discrete subspaces
temporalDim = 100;
spatialDim  = 20;
% number of training parameters to use
paramNum    = 10000;

% model data
pd = ProblemData;
pd.laplacian = 1;
pd.offset    = 1;
pd.xspan     = [0 1];
pd.tspan     = [0 1];
pd.xgrid     = linspace(pd.xspan(1), pd.xspan(2), 100);
pd.tgrid     = linspace(pd.tspan(1), pd.tspan(2), temporalDim);
pd.f         = [];
pd.nF        = 1;
pd.nC        = 1;
pd.seriesIdx = [1];
pd.sourcefun = 0;
pd.icfun     = @(x) sin(pi * x);

% parameter domain span
pspan = [-1, 1];

%% create the assembly objects
temporal = TemporalAssemblyNodal(pd, 0);
spatial  = SpatialAssemblyFourier(pd, spatialDim);
solver   = SolverNodal(pd, spatial, temporal);

%% create the parameter training set.
% two possibilities:
% first, equidistant grid
ptrain = linspace(pspan(1), pspan(2), paramNum);
% second, randomly chosen
% ptrain = rand(1, paramNum) * (pspan(2) - pspan(1)) - (pspan(2) - pspan(1)) / 2;

%% and now the real stuff: first let's execute the scm offline part
scm = SCM(solver, ptrain);
scm.startOfflineStage();

%% and now the reduced basis offline part
rbm = RBM(pd, solver, scm);
rbm.resetTraining;
rbm.offlineStage(ptrain);

