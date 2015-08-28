% Chapter 5 - First example.
%
% This example corresponds to the first example given in chapter 4. We use
% homogeneous boundary conditions in space and the remaining model data is given
% by:
%   * `\Omega = [0, 1]`, `I = [0, 1]`,
%   * the parametrich differential operator `A(\sigma)` is given by
%       `A \eta = - \frac{1}{10} \Delta \eta + \eta + \frac{1}{2} \sigma \sin(\pi \blank) \eta.`,
%   * where the parameter domain is given by `\mathcal P = [-1, 1]`,
%   * the source term is `g = 0` and the initial conditions are given by
%       `u_0 = sin(\pi \blank)`.

% dimensions of the discrete subspaces
temporalDim = 50;
spatialDim  = 10;
% number of training parameters to use
paramNum    = 100;
% parameter domain span
pspan = [-1, 1];

% model data
pd = ProblemData;
pd.laplacian        = 0.1;
pd.offset           = 1;
pd.xspan            = [0 1];
pd.tspan            = [0 1];
pd.xgrid            = linspace(pd.xspan(1), pd.xspan(2), 100);
pd.tgrid            = linspace(pd.tspan(1), pd.tspan(2), temporalDim);
pd.f                = [];
pd.nF               = 1;
pd.nC               = 1;
pd.seriesIdx        = [1];
pd.useSineExpansion = true;
pd.sourcefun        = 0;
pd.icfun            = @(x) 0.5 * sin(pi * x);
