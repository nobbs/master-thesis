% Load data that was generated for thesis.
% The used problem data is the following:
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

% further we calculated the inf-sup-constant `\beta_{\mathcal N}` and the
% correspoding bound for the discretizations with the following values:
% number of temporal grid points
tvalues  = [11, 26, 51, 76, 101, 126, 151, 176, 201, 226, 251, 276, 301];
% number of spatial basis functions to use
svalues  = [5, 10, 15, 20, 25];

% the used spatial basis functions are fourier basis functions given through
% SpatialAssemblyFourier

% now, let's load the precomputed data
load('stabilityFourier1.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure(1);
semilogy(tvalues, exact, '-o')
ylim([1e-3, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

figure(2);
semilogy(tvalues, bounds, '-o')
ylim([1e-3, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');
