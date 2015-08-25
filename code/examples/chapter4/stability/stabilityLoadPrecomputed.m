%% Load data that was generated for thesis.
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

% the used spatial basis functions are sine basis functions given through
% SpatialAssemblySine

% now, let's load the precomputed data
load('stabilitySineData1.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure(1);
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_sine_dataset1_fig_1.tikz');

figure(2);
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_sine_dataset1_fig_2.tikz');

%% And now the second example
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

% now, let's load the precomputed data
load('stabilitySineData2.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure(3);
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_sine_dataset2_fig_1.tikz');

figure(4);
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_sine_dataset2_fig_2.tikz');

%% Load data that was generated for thesis.
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
load('stabilityFourierData.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure(5);
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_fourier_dataset1_fig_1.tikz');

figure(6);
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_fourier_dataset1_fig_2.tikz');

