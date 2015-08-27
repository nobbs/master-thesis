% Chapter 4 - First example, homogeneous boundary conditions in space.
%
% This example calculates the inf-sup-constant for the petrov-galerkin-
% discretization with different temporal and spatial dimensions. Additionally
% the inf-sup-bound given by Andreev is computed.
%
% This model uses the following data:
%   * `\Omega = [0, 1]`, `I = [0, 1]`,
%   * the differential operator `A` is given by
%       `A \eta = - \Delta \eta + \eta + \sin(\pi \blank) \eta.`.

%% setup
% set the values for the temporal and spatial dimensions
tvalues = [11,26:25:301];
svalues = [5:5:25];

% setup the problem data
pd = ProblemData;
pd.laplacian = 1;
pd.offset    = 1;
pd.xspan     = [0 1];
pd.tspan     = [0 1];
pd.xgrid     = linspace(pd.xspan(1), pd.xspan(2), 10);
pd.f         = [];
pd.nF        = 1;
pd.nC        = 1;
pd.seriesIdx = [1];

% testing parameters (sets the used field)
useSineExpansion = true;
fieldcoeffs      = [1];

%% execution
% get the number of values
nt = length(tvalues);
ns = length(svalues);

% create the shared assembly objects
SpatialAssembly  = @SpatialAssemblySine;

% and now run!
ch4_compute_stability;
