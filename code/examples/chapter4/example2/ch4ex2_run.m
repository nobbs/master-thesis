% Chapter 4 - Second example, homogeneous boundary conditions in space.
%
% This example calculates the inf-sup-constant for the petrov-galerkin-
% discretization with different temporal and spatial dimensions. Additionally
% the inf-sup-bound given by Andreev is computed.
%
% This model uses the following data:
%   * `\Omega = [0, 1]`, `I = [0, 1]`,
%   * the differential operator `A(t)` is given by
%       `A(t)\eta = - \Delta \eta + 3 \eta + \chi_{[0, 0.5)}(t) \left[ \sin(\pi \blank) -
%       \sin(6 \pi \blank) \right] \eta + \chi_{[0.5, 1]}(t) \left[ \sin(3 \pi \blank ) +
%       \sin(6 \pi \blank) \right] \eta.`.

%% setup
% set the values for the temporal and spatial dimensions
% (odd temporal dimensions to guarantee that pd.f is a gridpoint)
tvalues = [11 25 51 75 101 125 151 175 201 225 251 275 301];
svalues = [5:5:25];

% setup the problem data
pd = ProblemData;
pd.laplacian = 1;
pd.offset    = 3;
pd.xspan     = [0 1];
pd.tspan     = [0 1];
pd.xgrid     = linspace(pd.xspan(1), pd.xspan(2), 10);
pd.f         = [1/2];
pd.nF        = 2;
pd.nC        = 3;
pd.seriesIdx = [1,3,6];

% testing parameters
useSineExpansion = true;
fieldcoeffs = [1, 0, -1, 0, 1, 1];

%% execution
% number of elements
nt = length(tvalues);
ns = length(svalues);

% create the shared assembly objects
SpatialAssembly  = @SpatialAssemblySine;

% and now run!
ch4_compute_stability;
