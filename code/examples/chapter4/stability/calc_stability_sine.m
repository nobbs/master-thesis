% settings
tvalues = [11, 51:50:501];
svalues = [10:10:50];

% number of elements
nt = length(tvalues);
ns = length(svalues);

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

% testing parameters
fieldcoeffs = [1];

stability;

%% prepare for plotting
% [tgrid, sgrid] = meshgrid(tvalues, svalues);
% figure();
% mesh(sgrid, tgrid, exact.');
% hold on;
% % figure();
% mesh(sgrid, tgrid, bounds.');

%% series plots

figure();
semilogy(tvalues, exact, '-o')
ylim([1e-4, 1]);
figure();
semilogy(tvalues, bounds, '-o')
ylim([1e-4, 1]);
% figure();
% plot(cfls)

%% series plots - alternative
figure()
semilogy(svalues, exact.')
figure()
semilogy(svalues, bounds.')