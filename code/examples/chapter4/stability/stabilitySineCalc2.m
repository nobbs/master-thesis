% settings
tvalues = [11 25 51 75 101 125 151 175 201 225 251 275 301];
svalues = [5:5:25];

% number of elements
nt = length(tvalues);
ns = length(svalues);

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
fieldcoeffs = [1, 0, -1, 0, 1, 1];

SpatialAssembly = @SpatialAssemblySine;
useSineExpansion = true;

%% run
stability;

%% save the run
save('sineset_normal.mat');
clear
