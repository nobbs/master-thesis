% settings
tvalues = 5:5:105;
svalues = 1:2:20;

% number of elements
nt = length(tvalues);
ns = length(svalues);

% setup the problem data
pd = ProblemData;
pd.laplacian = 1;
pd.offset    = 1;
pd.xspan     = [0 1];
pd.tspan     = [0 1];
pd.xgrid     = linspace(pd.xspan(1), pd.xspan(2), 100);
pd.f         = [];
pd.nF        = 1;
pd.nC        = 1;
pd.seriesIdx = [1];

% testing parameters
fieldcoeffs = [1];

stability
