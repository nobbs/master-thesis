% Chapter 5 - Second example.
%
% We use homogeneous boundary conditions in space and the remaining model data
% is given by:
%   * `\Omega = [0, 1]`, `I = [0, 1]`,
%   * the parametrical differential operator `A(\sigma)` is given by
%       `A \eta = - \frac{1}{10} \Delta \eta + 2\eta + \chi_{[0, 0.5)} [\sigma_1
%       \sin(\pi \blank) + \sigma_2 \sin(2 \pi \blank)] \eta + \chi_{[0.5, 1]}
%       [\sigma_3 \sin(\pi \blank) + \sigma_4 \sin(2 \pi \blank)] \eta`,
%   * where the parameter domain is given by `\mathcal P = [-1, 1]^4`,
%   * the source term is `g = 0` and the initial conditions are given by
%       `u_0 = 0.5 sin(\pi \blank)`.

% dimensions of the discrete subspaces
temporalDim = 200;
spatialDim  = 20;
% number of training parameters to use
paramNum    = 5000;
% parameter domain span
pspan = [-1, 1];

% model data
pd = ProblemData;
pd.laplacian        = 0.1;
pd.offset           = 2;
pd.xspan            = [0 1];
pd.tspan            = [0 1];
pd.xgrid            = linspace(pd.xspan(1), pd.xspan(2), 100);
pd.tgrid            = linspace(pd.tspan(1), pd.tspan(2), temporalDim);
pd.f                = [1/2];
pd.nF               = 2;
pd.nC               = 2;
pd.seriesIdx        = [1, 2];
pd.useSineExpansion = true;
pd.sourcefun        = 0;
pd.icfun            = @(x) 0.5 * sin(pi * x);
